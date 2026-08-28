/*
Copyright (C) 2017  Liangliang Nan
https://3d.bk.tudelft.nl/liangliang/ - liangliang.nan@gmail.com

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

#include "linear_program_solver.h"
#include "../basic/basic_types.h"

#ifdef HAS_HIGHS

#include "Highs.h"

#include <cstdlib>
#include <iostream>


bool LinearProgramSolver::_solve_HIGHS(const LinearProgram* program) {
	try {
		typedef Variable<double>			Variable;
		typedef LinearExpression<double>	Objective;
		typedef LinearConstraint<double>	Constraint;

		const std::vector<Variable>& variables = program->variables();
		if (variables.empty()) {
			std::cerr << "variable set is empty" << std::endl;
			return false;
		}

		HighsLp lp;
		lp.num_col_ = static_cast<HighsInt>(variables.size());
		lp.sense_ = program->objective_sense() == LinearProgram::MINIMIZE
			? ObjSense::kMinimize : ObjSense::kMaximize;
		lp.offset_ = 0.0;

		lp.col_cost_.assign(variables.size(), 0.0);
		lp.col_lower_.resize(variables.size());
		lp.col_upper_.resize(variables.size());
		lp.integrality_.resize(variables.size());
		bool has_integral = false;
		for (std::size_t i = 0; i < variables.size(); ++i) {
			const Variable& var = variables[i];
			double lb, ub;
			var.get_double_bounds(lb, ub);
			lp.col_lower_[i] = lb <= -DBL_MAX ? -kHighsInf : lb;
			lp.col_upper_[i] = ub >= DBL_MAX ? kHighsInf : ub;
			switch (var.variable_type())
			{
			case Variable::CONTINUOUS:
				lp.integrality_[i] = HighsVarType::kContinuous;
				break;
			case Variable::INTEGER:
				lp.integrality_[i] = HighsVarType::kInteger;
				has_integral = true;
				break;
			case Variable::BINARY:
				lp.integrality_[i] = HighsVarType::kInteger;
				has_integral = true;
				break;
			}
		}
		// an all-continuous model must not carry integrality data
		if (!has_integral)
			lp.integrality_.clear();

		// objective coefficients
		const std::unordered_map<std::size_t, double>& obj_coeffs = program->objective().coefficients();
		for (const auto& c : obj_coeffs)
			lp.col_cost_[c.first] = c.second;

		// constraints: build the column-wise sparse matrix
		const std::vector<Constraint>& constraints = program->constraints();
		lp.num_row_ = static_cast<HighsInt>(constraints.size());
		lp.row_lower_.resize(constraints.size());
		lp.row_upper_.resize(constraints.size());

		std::vector<HighsInt> col_count(variables.size(), 0);
		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const Constraint& cstr = constraints[i];
			for (const auto& c : cstr.coefficients())
				++col_count[c.first];

			double lb, ub;
			switch (cstr.bound_type()) {
			case Constraint::FIXED: {
				double v = cstr.get_single_bound();
				lb = ub = v;
				break;
			}
			case Constraint::LOWER: {
				lb = cstr.get_single_bound();
				ub = kHighsInf;
				break;
			}
			case Constraint::UPPER: {
				ub = cstr.get_single_bound();
				lb = -kHighsInf;
				break;
			}
			case Constraint::DOUBLE:
			default:
				cstr.get_double_bounds(lb, ub);
				break;
			}
			lp.row_lower_[i] = lb <= -DBL_MAX ? -kHighsInf : lb;
			lp.row_upper_[i] = ub >= DBL_MAX ? kHighsInf : ub;
		}

		HighsSparseMatrix& a = lp.a_matrix_;
		a.format_ = MatrixFormat::kColwise;
		std::size_t num_nz = 0;
		for (std::size_t j = 0; j < col_count.size(); ++j)
			num_nz += col_count[j];
		a.start_.resize(variables.size() + 1);
		a.start_[0] = 0;
		for (std::size_t j = 0; j < col_count.size(); ++j)
			a.start_[j + 1] = a.start_[j] + col_count[j];
		// current fill position per column
		std::vector<HighsInt> fill = a.start_;
		a.index_.resize(num_nz);
		a.value_.resize(num_nz);
		for (std::size_t i = 0; i < constraints.size(); ++i) {
			const auto& coeffs = constraints[i].coefficients();
			for (const auto& c : coeffs) {
				HighsInt pos = fill[c.first]++;
				a.index_[pos] = static_cast<HighsInt>(i);
				a.value_[pos] = c.second;
			}
		}

		Highs highs;
		// quiet: the CLI/GUI print their own progress lines
		highs.setOptionValue("output_flag", false);
		// same settings as the SCIP backend: 600 s wall limit, 1e-4 MIP gap
		highs.setOptionValue("time_limit", 600.0);
		highs.setOptionValue("mip_rel_gap", 1e-4);
		// experiment knob: looser MIP gap trades optimality proof for speed
		if (const char* gap_env = std::getenv("CITY3D_HIGHS_GAP"))
			highs.setOptionValue("mip_rel_gap", std::atof(gap_env));
		// optional thread-count override for experiments (unset/0 = HiGHS default).
		// Measured on case7: HiGHS MIP branch & bound is serial — 4/8 threads make
		// no difference, 16 threads oversubscribe and hurt 3x. Keep default.
		if (const char* threads_env = std::getenv("CITY3D_HIGHS_THREADS"))
			highs.setOptionValue("threads", std::atoi(threads_env));

		if (highs.passModel(lp) != HighsStatus::kOk) {
			std::cerr << "HiGHS failed to accept the model" << std::endl;
			return false;
		}

		if (highs.run() != HighsStatus::kOk) {
			std::cerr << "HiGHS failed to solve the model: "
				<< highs.modelStatusToString(highs.getModelStatus()) << std::endl;
			return false;
		}

		const HighsModelStatus status = highs.getModelStatus();
		if (status != HighsModelStatus::kOptimal && status != HighsModelStatus::kTimeLimit) {
			std::cerr << "HiGHS model status: " << highs.modelStatusToString(status) << std::endl;
			return false;
		}

		const HighsSolution& sol = highs.getSolution();
		if (!sol.value_valid || sol.col_value.size() != variables.size()) {
			if (status == HighsModelStatus::kTimeLimit)
				std::cerr << "time limit reached" << std::endl;
			return false;
		}

		objective_value_ = highs.getInfo().objective_function_value;
		result_.assign(sol.col_value.begin(), sol.col_value.end());
		return true;
	}
	catch (std::exception e) {
		std::cerr << "Error code = " << e.what() << std::endl;
	}
	catch (...) {
		std::cerr << "Exception during optimization" << std::endl;
	}
	return false;
}

#endif // HAS_HIGHS
