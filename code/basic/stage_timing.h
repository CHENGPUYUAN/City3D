/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/ - liangliang.nan@gmail.com

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef _STAGE_TIMING_H_
#define _STAGE_TIMING_H_

#include "stop_watch.h"

#include <map>
#include <mutex>
#include <sstream>
#include <string>

// Lightweight per-stage wall-time accumulator for the reconstruction pipeline.
// The pipeline itself is serial; the mutex only keeps stage scopes safe if
// individual stages ever become internally parallel.
class StageTiming
{
public:
	static void record(const std::string& stage, double seconds) {
		std::lock_guard<std::mutex> guard(mutex());
		table()[stage] += seconds;
	}

	static void clear() {
		std::lock_guard<std::mutex> guard(mutex());
		table().clear();
	}

	// one line, stages in alphabetical order, e.g. "lp_solve=35.2s lp_build=1.1s ..."
	static std::string summary() {
		std::lock_guard<std::mutex> guard(mutex());
		std::ostringstream ss;
		for (const auto& entry : table())
			ss << entry.first << "=" << entry.second << "s ";
		return ss.str();
	}

private:
	static std::mutex& mutex() { static std::mutex m; return m; }
	static std::map<std::string, double>& table() { static std::map<std::string, double> t; return t; }
};

// RAII helper: measures the wall time of the enclosing scope.
//   { StageScope scope("lp_solve"); ... }
class StageScope
{
public:
	StageScope(const std::string& name) : name_(name) {}
	~StageScope() { StageTiming::record(name_, watch_.seconds()); }

private:
	std::string name_;
	StopWatch   watch_;
};

#endif // _STAGE_TIMING_H_
