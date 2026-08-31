/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang.nan@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

// cuOpt backend: ships the model as JSON to the cuOpt HTTP server on the LAN
// GPU box (see docs/ops/cuopt-server-deploy.md) and polls for the solution.
// The server's "cuOpt_LP" action handles both LP and MILP: the request goes
// through MILP whenever "variable_types" is present. City3D's face-selection
// model is binary, so we always emit variable_types ("I" for BINARY/INTEGER,
// "C" for CONTINUOUS) and let cuOpt run its (beta) MILP solver.
//
// No HTTP/JSON dependencies: the model schema is trivial to serialize by hand
// (same approach as the MINDOPTPY backend) and a one-shot HTTP/1.1 exchange
// with "Connection: close" needs only POSIX sockets.

#include "linear_program_solver.h"
#include "../basic/logger.h"

#include <cctype>
#include <cfloat>
#include <cerrno>
#include <cstdlib>
#include <cstring>

#include <arpa/inet.h>
#include <fcntl.h>
#include <netdb.h>
#include <netinet/in.h>
#include <poll.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <unistd.h>

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <thread>
#include <vector>


namespace {

// Server address; override with CITY3D_CUOPT_URL (e.g. http://10.0.0.5:8900).
const char* CUOPT_URL_DEFAULT = "http://192.168.2.156:8900";
// Sent to the server in solver_config.time_limit; matches the hard limit of
// the in-process backends. The client waits well past it because the server
// clock starts only after its own presolve (which alone can take a minute on
// the DSM models) and encoding the result vector takes tens of seconds more;
// giving up early would waste the whole solve and trigger a needless local
// re-solve.
const double SOLVE_TIME_LIMIT_SEC = 600.0;
const double CLIENT_DEADLINE_SEC = SOLVE_TIME_LIMIT_SEC + 300.0;
const double POLL_INTERVAL_SEC = 0.5;

bool get_server_host_port(std::string& host, std::string& port) {
    const char* url = std::getenv("CITY3D_CUOPT_URL");
    std::string u = url ? url : CUOPT_URL_DEFAULT;
    const std::string prefix = "http://";
    if (u.rfind(prefix, 0) != 0) {
        std::cerr << "CITY3D_CUOPT_URL must start with " << prefix << std::endl;
        return false;
    }
    u = u.substr(prefix.size());
    // drop any path component; the endpoints are fixed
    std::size_t slash = u.find('/');
    if (slash != std::string::npos) u = u.substr(0, slash);
    std::size_t colon = u.rfind(':');
    if (colon == std::string::npos) {
        host = u;
        port = "80";
    } else {
        host = u.substr(0, colon);
        port = u.substr(colon + 1);
    }
    return !host.empty();
}

// One-shot HTTP/1.1 request/response over a fresh TCP connection.
bool http_call(const std::string& host, const std::string& port,
               const char* method, const std::string& path,
               const std::string& body, std::string& response) {
    addrinfo hints{};
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;
    addrinfo* res = nullptr;
    if (getaddrinfo(host.c_str(), port.c_str(), &hints, &res) != 0 || !res) {
        std::cerr << "cuOpt: cannot resolve " << host << std::endl;
        return false;
    }

    bool connected = false;
    int fd = -1;
    for (addrinfo* p = res; p && !connected; p = p->ai_next) {
        fd = ::socket(p->ai_family, p->ai_socktype, p->ai_protocol);
        if (fd < 0) continue;
        int flags = fcntl(fd, F_GETFL, 0);
        // bounded connect so a down/unreachable server fails in seconds,
        // not at the TCP default of ~2 minutes
        fcntl(fd, F_SETFL, flags | O_NONBLOCK);
        if (::connect(fd, p->ai_addr, p->ai_addrlen) != 0 && errno == EINPROGRESS) {
            pollfd pfd{fd, POLLOUT, 0};
            int err = 0;
            socklen_t len = sizeof(err);
            if (poll(&pfd, 1, 3000) == 1 &&
                getsockopt(fd, SOL_SOCKET, SO_ERROR, &err, &len) == 0)
                connected = (err == 0);
        } else {
            connected = true;
        }
        if (connected) fcntl(fd, F_SETFL, flags); // back to blocking
        else { ::close(fd); fd = -1; }
    }
    freeaddrinfo(res);
    if (!connected) {
        std::cerr << "cuOpt: cannot connect to " << host << ":" << port << std::endl;
        return false;
    }

    std::ostringstream req;
    req << method << " " << path << " HTTP/1.1\r\n"
        << "Host: " << host << ":" << port << "\r\n"
        << "CLIENT-VERSION: custom\r\n"
        << "Connection: close\r\n";
    if (!body.empty()) {
        req << "Content-Type: application/json\r\n"
            << "Content-Length: " << body.size() << "\r\n";
    }
    req << "\r\n" << body;

    const std::string req_str = req.str();
    std::size_t sent = 0;
    while (sent < req_str.size()) {
        ssize_t n = ::send(fd, req_str.data() + sent, req_str.size() - sent, 0);
        if (n <= 0) {
            ::close(fd);
            std::cerr << "cuOpt: send failed" << std::endl;
            return false;
        }
        sent += static_cast<std::size_t>(n);
    }

    timeval tv{30, 0}; // per-recv; the server answers each request promptly
    setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
    char buf[65536];
    response.clear();
    while (true) {
        ssize_t n = ::recv(fd, buf, sizeof(buf), 0);
        if (n < 0) {
            if (errno == EINTR) continue;
            ::close(fd);
            std::cerr << "cuOpt: receive failed or timed out" << std::endl;
            return false;
        }
        if (n == 0) break;
        response.append(buf, static_cast<std::size_t>(n));
    }
    ::close(fd);

    // strip headers
    std::size_t hdr_end = response.find("\r\n\r\n");
    if (hdr_end == std::string::npos) {
        std::cerr << "cuOpt: malformed HTTP response" << std::endl;
        return false;
    }
    response = response.substr(hdr_end + 4);
    return true;
}

void append_num(std::ostringstream& os, double v) {
    // JSON has no Infinity; cuOpt accepts the "inf"/"ninf" strings instead
    if (v >= DBL_MAX) os << "\"inf\"";
    else if (v <= -DBL_MAX) os << "\"ninf\"";
    else os << v;
}

// Locates "key": in the JSON and positions just past the colon (the server
// writes "key": value with a space, so whitespace after ':' is expected).
bool find_key_pos(const std::string& json, const char* key, std::size_t& pos) {
    std::string pat = std::string("\"") + key + "\":";
    std::size_t p = json.find(pat);
    if (p == std::string::npos) return false;
    pos = p + pat.size();
    while (pos < json.size() && std::isspace(static_cast<unsigned char>(json[pos]))) ++pos;
    return pos < json.size();
}

// Extracts the string value of "key":"..." (first occurrence).
bool find_string_value(const std::string& json, const char* key, std::string& value) {
    std::size_t p;
    if (!find_key_pos(json, key, p) || json[p] != '"') return false;
    ++p;
    std::size_t e = json.find('"', p);
    if (e == std::string::npos) return false;
    value = json.substr(p, e - p);
    return true;
}

// Parses the "[a, b, ...]" array following "key": into doubles.
// Returns false while the value is not a materialized numeric array (e.g. a
// poll that arrives before the result store finishes encoding, or the
// "primal_solution":null placeholder).
bool find_number_array(const std::string& json, const char* key, std::vector<double>& out) {
    std::size_t p;
    if (!find_key_pos(json, key, p) || json[p] != '[') return false;
    ++p; // past '['

    out.clear();
    const char* s = json.c_str() + p;
    char* end = nullptr;
    while (true) {
        while (*s == ' ' || *s == ',') ++s;
        if (*s == ']') break;
        double v = std::strtod(s, &end);
        if (end == s) return false; // null/NaN/nested — not a ready array
        out.push_back(v);
        s = end;
        while (*s == ' ' || *s == ',') ++s;
        if (*s == ']') break;
        if (*s == '\0') return false; // truncated response
    }
    return true;
}

bool find_number(const std::string& json, const char* key, double& value) {
    std::size_t p;
    if (!find_key_pos(json, key, p)) return false;
    const char* s = json.c_str() + p;
    char* end = nullptr;
    double v = std::strtod(s, &end);
    if (end == s) return false;
    value = v;
    return true;
}

} // namespace


bool LinearProgramSolver::_solve_CUOPT(const LinearProgram* program) {
    typedef Variable<double>			Variable;
    typedef LinearExpression<double>	Objective;
    typedef LinearConstraint<double>	Constraint;

    const std::vector<Variable>& variables = program->variables();
    if (variables.empty()) {
        std::cerr << "variable set is empty" << std::endl;
        return false;
    }
    const std::vector<Constraint>& constraints = program->constraints();
    const std::size_t n = variables.size();
    const std::size_t m = constraints.size();

    std::string host, port;
    if (!get_server_host_port(host, port)) return false;

    // ---- serialize the model (cuOpt server LP/MILP schema) ----
    std::ostringstream jos;
    jos << std::setprecision(17);
    jos << "{\"action\":\"cuOpt_LP\",\"client_version\":\"custom\",\"data\":{";

    // constraints as a CSR matrix; sorted indices so the request (and thus the
    // solver's floating-point summation) doesn't depend on hash-map order
    jos << "\"csr_constraint_matrix\":{\"offsets\":[";
    std::vector<std::size_t> offsets;
    offsets.reserve(m + 1);
    offsets.push_back(0);
    std::ostringstream jidx, jval;
    jidx << std::setprecision(17);
    jval << std::setprecision(17);
    for (std::size_t i = 0; i < m; ++i) {
        std::vector<std::pair<std::size_t, double>> coeffs(
            constraints[i].coefficients().begin(), constraints[i].coefficients().end());
        std::sort(coeffs.begin(), coeffs.end());
        for (const auto& c : coeffs) {
            jidx << c.first << ",";
            jval << c.second << ",";
        }
        offsets.push_back(offsets.back() + coeffs.size());
    }
    for (std::size_t i = 0; i < offsets.size(); ++i) {
        if (i) jos << ",";
        jos << offsets[i];
    }
    jos << "],\"indices\":[";
    {
        std::string s = jidx.str();
        if (!s.empty()) s.pop_back(); // trailing comma
        jos << s;
    }
    jos << "],\"values\":[";
    {
        std::string s = jval.str();
        if (!s.empty()) s.pop_back();
        jos << s;
    }
    jos << "]}";

    // row bounds
    jos << ",\"constraint_bounds\":{\"upper_bounds\":[";
    for (std::size_t i = 0; i < m; ++i) {
        if (i) jos << ",";
        double lb = -DBL_MAX, ub = DBL_MAX;
        switch (constraints[i].bound_type()) {
        case Constraint::FIXED: {
            double v = constraints[i].get_single_bound();
            lb = ub = v;
            break;
        }
        case Constraint::LOWER: lb = constraints[i].get_single_bound(); break;
        case Constraint::UPPER: ub = constraints[i].get_single_bound(); break;
        case Constraint::DOUBLE: constraints[i].get_double_bounds(lb, ub); break;
        default: break;
        }
        append_num(jos, ub);
    }
    jos << "],\"lower_bounds\":[";
    for (std::size_t i = 0; i < m; ++i) {
        if (i) jos << ",";
        double lb = -DBL_MAX, ub = DBL_MAX;
        switch (constraints[i].bound_type()) {
        case Constraint::FIXED: {
            double v = constraints[i].get_single_bound();
            lb = ub = v;
            break;
        }
        case Constraint::LOWER: lb = constraints[i].get_single_bound(); break;
        case Constraint::UPPER: ub = constraints[i].get_single_bound(); break;
        case Constraint::DOUBLE: constraints[i].get_double_bounds(lb, ub); break;
        default: break;
        }
        append_num(jos, lb);
    }
    jos << "]}";

    // objective as a dense coefficient vector
    jos << ",\"objective_data\":{\"coefficients\":[";
    {
        std::vector<double> obj(n, 0.0);
        for (const auto& c : program->objective().coefficients())
            obj[c.first] = c.second;
        for (std::size_t i = 0; i < n; ++i) {
            if (i) jos << ",";
            jos << obj[i];
        }
    }
    jos << "]}";

    // variable bounds + types (types present => the server runs MILP)
    jos << ",\"variable_bounds\":{\"upper_bounds\":[";
    for (std::size_t i = 0; i < n; ++i) {
        if (i) jos << ",";
        double lb, ub;
        variables[i].get_double_bounds(lb, ub);
        append_num(jos, ub);
    }
    jos << "],\"lower_bounds\":[";
    for (std::size_t i = 0; i < n; ++i) {
        if (i) jos << ",";
        double lb, ub;
        variables[i].get_double_bounds(lb, ub);
        append_num(jos, lb);
    }
    jos << "]}";
    jos << ",\"variable_types\":[";
    for (std::size_t i = 0; i < n; ++i) {
        if (i) jos << ",";
        jos << "\"" << (variables[i].variable_type() == Variable::CONTINUOUS ? "C" : "I") << "\"";
    }
    jos << "]";

    jos << ",\"maximize\":" << (program->objective_sense() == LinearProgram::MAXIMIZE ? "true" : "false");
    jos << ",\"solver_config\":{\"time_limit\":" << SOLVE_TIME_LIMIT_SEC << "}";
    jos << "}}";

    // ---- submit ----
    std::string response;
    if (!http_call(host, port, "POST", "/cuopt/request", jos.str(), response)) return false;
    std::string req_id;
    if (!find_string_value(response, "reqId", req_id)) {
        std::cerr << "cuOpt: no reqId in server reply" << std::endl;
        return false;
    }

    // ---- poll for the solution ----
    const std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::now() +
        std::chrono::duration_cast<std::chrono::steady_clock::duration>(
            std::chrono::duration<double>(CLIENT_DEADLINE_SEC));
    const std::string sol_path = "/cuopt/solution/" + req_id;
    int consecutive_failures = 0;
    while (true) {
        std::this_thread::sleep_for(
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::duration<double>(POLL_INTERVAL_SEC)));

        if (!http_call(host, port, "GET", sol_path, "", response)) {
            // a crashed solver worker respawns under the same server and can
            // refuse connections for a moment; tolerate ~60s of that before
            // declaring the backend down
            if (++consecutive_failures * POLL_INTERVAL_SEC > 60.0) return false;
            continue;
        }
        consecutive_failures = 0;

        // accept any status once a full primal solution vector is present
        // (MILP "FeasibleFound"/time-limit incumbents are usable, matching
        // the HiGHS backend's kTimeLimit handling)
        std::vector<double> primal;
        if (find_number_array(response, "primal_solution", primal)) {
            if (primal.size() != n) {
                std::cerr << "cuOpt: solution size mismatch: expected " << n
                          << ", got " << primal.size() << std::endl;
                return false;
            }
            result_ = primal;
            if (!find_number(response, "primal_objective", objective_value_)) {
                // fall back to recomputing from the objective coefficients
                std::vector<double> obj(n, 0.0);
                for (const auto& c : program->objective().coefficients())
                    obj[c.first] = c.second;
                objective_value_ = 0.0;
                for (std::size_t i = 0; i < n; ++i)
                    objective_value_ += obj[i] * primal[i];
            }
            return true;
        }

        std::string status;
        if (find_string_value(response, "status", status)) {
            if (status == "Infeasible" || status == "Unbounded") {
                std::cerr << "cuOpt: problem is " << status << std::endl;
                return false;
            }
            // "Optimal"/"FeasibleFound" can show up before the result store
            // finishes encoding the vectors — keep polling in that case
        }

        if (std::chrono::steady_clock::now() > deadline) {
            std::cerr << "cuOpt: no solution within " << CLIENT_DEADLINE_SEC
                      << "s (reqId " << req_id << ")" << std::endl;
            return false;
        }
    }
}
