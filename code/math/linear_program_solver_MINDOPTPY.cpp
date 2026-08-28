/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/ - liangliang.nan@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include "linear_program_solver.h"
#include "../basic/basic_types.h"
#include "../basic/logger.h"

#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

// hand-rolled serialization below; the model is small and the format trivial,
// so we don't pull a JSON dependency into the math module

namespace {

// The local MindOpt 2.3.0 SDK (pip's mindoptpy embeds an expired 2023-2024
// license, so the stand-alone SDK + community license is the only working
// route). LIBMILP_PATH locates the MILP component; the license file is picked
// up from ~/mindopt/ce_license.ini by default.
// NOTE: the shipped libmindopt.so.2.3.0 carries a RWE GNU_STACK flag that
// Ubuntu 26.04 refuses to load ("cannot enable executable stack"); both the
// lib/ and the lib/python/mindoptpy/ copies were patched in place (PT_GNU_STACK
// PF_X cleared, .bak kept). An upstream reinstall will need the same patch.
const char* MINDOPT_SDK_PY = "/home/nolan/mindopt/2.3.0/linux64-x86/lib/python";
const char* MINDOPT_LIB_DIR = "/home/nolan/mindopt/2.3.0/linux64-x86/lib";

} // namespace

// A thin backend that ships the model to a mindoptpy child process:
//   C++ writes model.json  ->  `uv run --with-editable <sdk> python <script> <model.json>`
//   -> mindoptpy (community license: no size limit, no time limit)
//   -> writes solution to <model.json>.sol: first line "obj <value>", then one
//      "<index> <value>" line per variable.
// Unlike the gurobipy free license there is no 2000x2000 size limit, so no
// presolve-to-fit step is needed here.
bool LinearProgramSolver::_solve_MINDOPTPY(const LinearProgram* program) {
	typedef Variable<double>			Variable;
	typedef LinearExpression<double>	Objective;
	typedef LinearConstraint<double>	Constraint;

	const std::vector<Variable>& variables = program->variables();
	if (variables.empty()) {
		std::cerr << "variable set is empty" << std::endl;
		return false;
	}

	// serialize the model as JSON
	std::ostringstream jos;
	jos << "{\"sense\":\"" << (program->objective_sense() == LinearProgram::MINIMIZE ? "min" : "max") << "\"";
	jos << ",\"variables\":[";
	for (std::size_t i = 0; i < variables.size(); ++i) {
		const Variable& var = variables[i];
		const char* vtype = var.variable_type() == Variable::BINARY ? "B"
			: var.variable_type() == Variable::INTEGER ? "I" : "C";
		double lb, ub;
		var.get_double_bounds(lb, ub);
		if (lb <= -DBL_MAX) lb = -1e100; // JSON-safe sentinels; mindoptpy treats them as infinite
		if (ub >= DBL_MAX) ub = 1e100;
		if (i) jos << ",";
		jos << "{\"t\":\"" << vtype << "\",\"lb\":" << lb << ",\"ub\":" << ub << "}";
	}
	jos << "],\"objective\":[";
	{
		const std::unordered_map<std::size_t, double>& coeffs = program->objective().coefficients();
		bool first = true;
		for (const auto& c : coeffs) {
			if (!first) jos << ",";
			first = false;
			jos << "[" << c.first << "," << c.second << "]";
		}
	}
	jos << "],\"constraints\":[";
	{
		const std::vector<Constraint>& constraints = program->constraints();
		bool first = true;
		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const Constraint& cstr = constraints[i];
			if (first) first = false;
			else jos << ",";
			jos << "{\"c\":[";
			const auto& coeffs = cstr.coefficients();
			bool f = true;
			for (const auto& c : coeffs) {
				if (!f) jos << ",";
				f = false;
				jos << "[" << c.first << "," << c.second << "]";
			}
			jos << "],";
			double lb = -1e100, ub = 1e100;
			switch (cstr.bound_type()) {
			case Constraint::FIXED: {
				double v = cstr.get_single_bound();
				lb = ub = v;
				break;
			}
			case Constraint::LOWER: lb = cstr.get_single_bound(); break;
			case Constraint::UPPER: ub = cstr.get_single_bound(); break;
			case Constraint::DOUBLE: cstr.get_double_bounds(lb, ub); break;
			default: break;
			}
			jos << "\"lb\":" << lb << ",\"ub\":" << ub << "}";
		}
	}
	jos << "]}";

	// embedded solver script: reads model file path on argv[1], writes the
	// solution to argv[1] + ".sol"
	static const char* SOLVER_SCRIPT = R"PY(
import json, os, sys

os.environ.setdefault("LIBMILP_PATH", "/home/nolan/mindopt/2.3.0/linux64-x86/lib")

import mindoptpy as mdo
from mindoptpy import Env, Model

path = sys.argv[1]
with open(path) as f:
    model = json.load(f)

sense = mdo.MDO.MINIMIZE if model["sense"] == "min" else mdo.MDO.MAXIMIZE
vtypes = {"B": "B", "I": "I", "C": "C"}

env = Env()
env.start()
m = Model("city3d", env)
m.setParam(mdo.MDO.Param.OutputFlag, 0)
m.setParam(mdo.MDO.Param.MaxTime, 600.0)

vars = model["variables"]
xs = []
for v in vars:
    lb = v["lb"] if v["lb"] > -1e99 else -mdo.MDO.INFINITY
    ub = v["ub"] if v["ub"] < 1e99 else mdo.MDO.INFINITY
    xs.append(m.addVar(lb=lb, ub=ub, vtype=vtypes[v["t"]]))

# mindoptpy's setObjective requires a LinExpr, not a bare Var/number
obj = mdo.quicksum(c * xs[int(i)] for i, c in model["objective"])
m.setObjective(obj, sense)

for con in model["constraints"]:
    expr = mdo.quicksum(coef * xs[int(v)] for v, coef in con["c"])
    lb = con["lb"] if con["lb"] > -1e99 else -mdo.MDO.INFINITY
    ub = con["ub"] if con["ub"] < 1e99 else mdo.MDO.INFINITY
    if lb == ub:
        m.addConstr(expr == lb)
    else:
        if lb > -mdo.MDO.INFINITY:
            m.addConstr(expr >= lb)
        if ub < mdo.MDO.INFINITY:
            m.addConstr(expr <= ub)

m.optimize()
if m.status != mdo.MDO.OPTIMAL:
    sys.exit(f"mindopt status {m.status} (1=optimal; see MDO status codes)")

with open(path + ".sol", "w") as f:
    # float(): m.objval and .X return np.float64 and numpy 2.x reprs
    # ("np.float64(...)") would break the C++ parser
    f.write(f"obj {float(m.objval)!r}\n")
    for i, x in enumerate(xs):
        f.write(f"{i} {float(x.X)!r}\n")
)PY";

	// model + script temp files
	char model_tmpl[] = "/tmp/city3d_lp_model_XXXXXX.json";
	char script_tmpl[] = "/tmp/city3d_lp_solver_XXXXXX.py";
	int model_fd = mkstemps(model_tmpl, 5);   // suffix ".json" (5 chars)
	int script_fd = mkstemps(script_tmpl, 3); // suffix ".py" (3 chars)
	if (model_fd < 0 || script_fd < 0) {
		if (model_fd >= 0) { close(model_fd); unlink(model_tmpl); }
		if (script_fd >= 0) { close(script_fd); unlink(script_tmpl); }
		std::cerr << "failed to create temp files for the mindoptpy backend" << std::endl;
		return false;
	}
	const std::string json_str = jos.str();
	const std::string script_str(SOLVER_SCRIPT);
	bool write_ok = write(model_fd, json_str.c_str(), json_str.size()) == (ssize_t)json_str.size()
		&& write(script_fd, script_str.c_str(), script_str.size()) == (ssize_t)script_str.size();
	close(model_fd);
	close(script_fd);
	if (!write_ok) {
		unlink(model_tmpl);
		unlink(script_tmpl);
		std::cerr << "failed to write the mindoptpy model/script" << std::endl;
		return false;
	}

	// run: uv run --with-editable <sdk> python <script> <model.json>
	// (the SDK is the stand-alone install, not pip's expired-license wheel)
	std::string sol_path = std::string(model_tmpl) + ".sol";
	std::string err_path = std::string(model_tmpl) + ".err";
	std::string cmd = "uv run --no-project --with-editable " + std::string(MINDOPT_SDK_PY)
		+ " python " + std::string(script_tmpl) + " " + model_tmpl + " 2>" + err_path + " >/dev/null";
	int rc = system(cmd.c_str());
	unlink(script_tmpl);

	if (rc != 0) {
		std::cerr << "mindoptpy failed (exit " << (rc >> 8) << ")" << std::endl;
		std::ifstream err(err_path.c_str());
		std::stringstream buf;
		buf << err.rdbuf();
		std::string msg = buf.str();
		if (!msg.empty())
			std::cerr << msg << std::endl;
		unlink(model_tmpl);
		unlink(err_path.c_str());
		return false;
	}

	// read back the solution
	std::ifstream sol(sol_path.c_str());
	if (!sol) {
		std::cerr << "mindoptpy did not produce a solution file" << std::endl;
		unlink(model_tmpl);
		unlink(err_path.c_str());
		return false;
	}
	std::string tag;
	sol >> tag >> objective_value_;
	if (!sol || tag != "obj") {
		// a malformed header means every subsequent extraction fails silently
		// and result_ stays all-zero — which reads as "select nothing" upstream
		std::cerr << "malformed solution header from mindoptpy" << std::endl;
		sol.close();
		unlink(model_tmpl);
		unlink(sol_path.c_str());
		unlink(err_path.c_str());
		return false;
	}
	result_.assign(variables.size(), 0.0);
	std::size_t idx;
	double val;
	std::size_t read = 0;
	while (sol >> idx >> val) {
		if (idx < result_.size())
			result_[idx] = val;
		++read;
	}
	if (read != variables.size()) {
		// partial reads (truncated file, stream error) would silently drop
		// variables to 0 and delete the corresponding faces
		std::cerr << "mindoptpy solution truncated: expected " << variables.size()
			<< " values, got " << read << std::endl;
		sol.close();
		unlink(model_tmpl);
		unlink(sol_path.c_str());
		unlink(err_path.c_str());
		return false;
	}
	sol.close();
	unlink(model_tmpl);
	unlink(sol_path.c_str());
	unlink(err_path.c_str());

	return true;
}
