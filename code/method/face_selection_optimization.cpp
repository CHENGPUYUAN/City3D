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

#include "face_selection.h"
#include "polyfit_info.h"
#include "method_global.h"
#include "../model/point_set.h"
#include "../model/map_editor.h"
#include "../model/map_geometry.h"
#include "../model/map_circulators.h"
#include "../basic/logger.h"
#include "../basic/stage_timing.h"
#include "../model/map_builder.h"
#include "../model/map_io.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <limits>

std::vector<std::vector<Map::Facet*> > FaceSelection::find_multi_roofs(Map* mesh,
	Map::Facet* footprint,
	std::vector<Plane3d*>& v)
{
	facet_attrib_supporting_plane_.bind_if_defined(model_, "FacetSupportingPlane");
	std::vector<Map::Facet*> faces;
	FOR_EACH_FACET(Map, mesh, it)
	{
		int count = 0; //count the testing times
		for (std::size_t k = 0; k < v.size(); k++)
		{
			const Plane3d* ver = v[k];
			if (facet_attrib_supporting_plane_[it] == ver)
				break;
			count++;
		}
		if (count == v.size())
			faces.push_back(it);
	}
	if (faces.size() < 2)
		return {};

	const Plane3d& plane = Geom::facet_plane(footprint);

	// Project each face onto the footprint plane once and cache its 2D polygon,
	// bbox, and centroid (the original code rebuilt every polygon for every
	// query point, i.e. O(F^2) polygon constructions).
	struct FaceInfo
	{
		Map::Facet* face;
		Polygon2d  plg;
		vec2       centroid;
		double     xmin, ymin, xmax, ymax;
	};
	std::vector<FaceInfo> infos;
	infos.reserve(faces.size());
	double gmin_x = std::numeric_limits<double>::max(), gmin_y = std::numeric_limits<double>::max();
	double gmax_x = -std::numeric_limits<double>::max(), gmax_y = -std::numeric_limits<double>::max();
	for (std::size_t i = 0; i < faces.size(); ++i)
	{
		FaceInfo info;
		info.face = faces[i];
		vec3 c(0, 0, 0);
		int degree = 0;
		FacetHalfedgeCirculator fcir(info.face);
		for (; !fcir->end(); ++fcir)
		{
			const vec3& q = fcir->halfedge()->vertex()->point();
			c += q;
			++degree;
			vec3 proj = plane.projection(q);
			info.plg.push_back(plane.to_2d(proj));
		}
		c /= degree;
		info.centroid = plane.to_2d(plane.projection(c));

		info.xmin = info.ymin = std::numeric_limits<double>::max();
		info.xmax = info.ymax = -std::numeric_limits<double>::max();
		for (std::size_t j = 0; j < info.plg.size(); ++j)
		{
			const vec2& p = info.plg[j];
			info.xmin = std::min(info.xmin, p.x);
			info.ymin = std::min(info.ymin, p.y);
			info.xmax = std::max(info.xmax, p.x);
			info.ymax = std::max(info.ymax, p.y);
		}
		gmin_x = std::min(gmin_x, info.xmin);
		gmin_y = std::min(gmin_y, info.ymin);
		gmax_x = std::max(gmax_x, info.xmax);
		gmax_y = std::max(gmax_y, info.ymax);
		infos.push_back(info);
	}

	// Uniform grid: a face is registered in every cell its bbox overlaps. A
	// point outside a face's closed bbox can never be inside its polygon, so
	// querying only the cell candidates yields exactly the same containment
	// results as testing all faces.
	const int grid_dim = std::max(1, std::min<int>(
		static_cast<int>(std::sqrt(static_cast<double>(infos.size()))), 512));
	const double cell_w = std::max((gmax_x - gmin_x) / grid_dim, 1e-12);
	const double cell_h = std::max((gmax_y - gmin_y) / grid_dim, 1e-12);
	auto cell_coord = [&](double x, double y, int& cx, int& cy) {
		cx = std::min(std::max(static_cast<int>((x - gmin_x) / cell_w), 0), grid_dim - 1);
		cy = std::min(std::max(static_cast<int>((y - gmin_y) / cell_h), 0), grid_dim - 1);
	};

	std::vector<std::vector<int> > cells(static_cast<std::size_t>(grid_dim) * grid_dim);
	for (std::size_t i = 0; i < infos.size(); ++i)
	{
		const FaceInfo& info = infos[i];
		int x0, y0, x1, y1;
		cell_coord(info.xmin, info.ymin, x0, y0);
		cell_coord(info.xmax, info.ymax, x1, y1);
		for (int cy = y0; cy <= y1; ++cy)
			for (int cx = x0; cx <= x1; ++cx)
				cells[static_cast<std::size_t>(cy) * grid_dim + cx].push_back(static_cast<int>(i));
	}

	std::vector<std::vector<Map::Facet*> > multiple_roofs;
	for (std::size_t i = 0; i < infos.size(); ++i)
	{
		const FaceInfo& info = infos[i];
		std::vector<Map::Facet*> roofs({ info.face });

		int cx, cy;
		cell_coord(info.centroid.x, info.centroid.y, cx, cy);
		const std::vector<int>& candidates = cells[static_cast<std::size_t>(cy) * grid_dim + cx];
		for (std::size_t k = 0; k < candidates.size(); ++k)
		{
			int j = candidates[k];
			if (j == static_cast<int>(i))
				continue;
			const FaceInfo& other = infos[j];
			if (info.centroid.x < other.xmin || info.centroid.x > other.xmax ||
				info.centroid.y < other.ymin || info.centroid.y > other.ymax)
				continue;
			if (Geom::point_is_in_polygon(other.plg, info.centroid))
				roofs.push_back(other.face);
		}

		if (roofs.size() > 1)
			multiple_roofs.push_back(roofs);
	}
	return multiple_roofs;
}

bool FaceSelection::optimize(PolyFitInfo* polyfit_info,
	Map::Facet* footprint,
	std::vector<Plane3d*>& v, LinearProgramSolver::SolverName solver_name)
{
	if (pset_ == 0 || model_ == 0)
		return false;

	facet_attrib_supporting_vertex_group_.bind_if_defined(model_, Method::facet_attrib_supporting_vertex_group);
	if (!facet_attrib_supporting_vertex_group_.is_bound())
	{
		Logger::err("-") << "attribute " << Method::facet_attrib_supporting_vertex_group << " doesn't exist" << std::endl;
        return false;
	}
	facet_attrib_supporting_point_num_.bind_if_defined(model_, Method::facet_attrib_supporting_point_num);
	if (!facet_attrib_supporting_point_num_.is_bound())
	{
		Logger::err("-") << "attribute " << Method::facet_attrib_supporting_point_num << " doesn't exist" << std::endl;
        return false;
	}
	facet_attrib_facet_area_.bind_if_defined(model_, Method::facet_attrib_facet_area);
	if (!facet_attrib_facet_area_.is_bound())
	{
		Logger::err("-") << "attribute " << Method::facet_attrib_facet_area << " doesn't exist" << std::endl;
        return false;
	}
	facet_attrib_covered_area_.bind_if_defined(model_, Method::facet_attrib_covered_area);
	if (!facet_attrib_covered_area_.is_bound())
	{
		Logger::err("-") << "attribute " << Method::facet_attrib_covered_area << " doesn't exist" << std::endl;
        return false;
	}

	facet_attrib_supporting_plane_.bind_if_defined(model_, "FacetSupportingPlane");

	MapFacetAttribute<bool> is_confident_facet(model_, "is_confident_facet");

	// Faces marked as having essentially no data support are excluded from
	// the program entirely (they get no variable). They stay in the mesh:
	// after solving, the closure repair pulls back the excluded faces that
	// border the selected surface — they are the structural filler the
	// "2 or 0" rule would otherwise have forced into the solution.
	MapFacetAttribute<bool> is_unsupported_facet;
	is_unsupported_facet.bind_if_defined(model_, "is_unsupported_facet");

	//////////////////////////////////////////////////////////////////////////

	double total_points = double(pset_->points().size());
	std::size_t idx = 0;
	MapFacetAttribute<std::size_t> facet_indices(model_);
	FOR_EACH_FACET(Map, model_, it)
	{
		Map::Facet* f = it;
		if (is_unsupported_facet.is_bound() && is_unsupported_facet[f])
			continue; // excluded: no variable
		facet_indices[f] = idx;
		++idx;
	}
    //StopWatch w;
	const std::vector<FaceStar>& fans = adjacency_.extract(model_, polyfit_info->planes);
	std::map<const FaceStar*, std::size_t> edge_sharp_status;    // the edge is sharp or not (needed after solving)

	MapFacetAttribute<bool> is_z_outside(model_, "is_z_outside_facet");
	// Faces whose centroid lies far outside the data's z range are
	// extrapolation artifacts (the arrangement projects every plane across
	// the whole footprint, so tilted planes grow corners far above/below
	// the real geometry), never part of a correct model: pin them to 0.
	// They stay in the arrangement and keep every fan intact.
	{
		const double data_zmax = pset_->bbox().z_max();
		const double data_zmin = pset_->bbox().z_min();
		const double z_tol = std::max(1.0, 0.2 * (data_zmax - data_zmin));
		// experiment switch: CITY3D_ZPIN=0 disables the pinning entirely
		// (diagnoses whether z-pinning — not the support-point marking —
		// starves regions of candidates and leaves holes)
		bool zpin_off = (std::getenv("CITY3D_ZPIN") != nil && std::getenv("CITY3D_ZPIN")[0] == '0');
		FOR_EACH_FACET(Map, model_, it)
		{
			Map::Facet* f = it;
			if (zpin_off)
			{
				is_z_outside[f] = false;
				continue;
			}
			// any vertex outside the band pins the face: long thin faces on
			// extrapolating planes can keep their centroid well inside the
			// band while a far corner sticks tens of meters out
			double fz_min = 1e30, fz_max = -1e30;
			FacetHalfedgeCirculator cir(f);
			for (; !cir->end(); ++cir)
			{
				double vz = cir->halfedge()->vertex()->point().z;
				fz_min = std::min(fz_min, vz);
				fz_max = std::max(fz_max, vz);
			}
			is_z_outside[f] = (fz_max > data_zmax + z_tol || fz_min < data_zmin - z_tol);
		}
	}
	auto is_excluded = [&](Map::Facet* f) {
		return (is_unsupported_facet.is_bound() && is_unsupported_facet[f]) ||
		       (is_z_outside.is_bound() && is_z_outside[f]);
	};
	// z-outside faces are pinned by geometry, not by missing support: the
	// closure repair must never pull them back in
	auto repair_eligible = [&](Map::Facet* f) {
		return is_excluded(f) && !(is_z_outside.is_bound() && is_z_outside[f]);
	};
	// number of selectable (non-excluded) members per fan
	std::vector<std::size_t> fan_active(fans.size(), 0);
	for (std::size_t i = 0; i < fans.size(); ++i)
	{
		const FaceStar& fan = fans[i];
		for (std::size_t j = 0; j < fan.size(); ++j)
		{
			if (!is_excluded(fan[j]->facet()))
				++fan_active[i];
		}
	}

	{
	StageScope stage("09a_lp_build");
	MapHalfedgeAttribute<bool> is_bound(model_);

	FOR_EACH_HALFEDGE(Map, model_, it)
	{
		auto edge = *it;
		is_bound[edge] = false;
	}
	for (std::size_t i = 0; i < fans.size(); ++i)
	{
		const FaceStar& star = fans[i];
		if (star.size() == 1)
		{
			Map::Halfedge* h = *star.begin();
			is_bound[h] = true;
		}
	}

	std::size_t num_faces = idx; // selectable faces only (excluded ones got no index)
	std::size_t num_usage_vars = 0;
	// Height reference from the DATA, not the mesh: the arrangement projects
	// every plane across the whole footprint, so extrapolating planes grow
	// corners far outside the real z range and blow up the mesh bbox — a
	// mesh-bbox reference made tall extrapolation faces (near bbox_zmax)
	// almost free for the height term, and the solver started picking them.
	double bbox_zmax = pset_->bbox().z_max();
	double bbox_zheight = bbox_zmax - pset_->bbox().z_min();
	// one usage variable per edge with more than two selectable incident
	// faces (the "2 or 0" constraint below needs it); edges degraded by
	// excluded faces are left unconstrained — the closure repair after
	// solving re-adds the filler faces such edges need
	std::map<const FaceStar*, std::size_t> edge_usage_status;    // keep or remove an intersecting edges
	for (std::size_t i = 0; i < fans.size(); ++i)
	{
		const FaceStar& fan = fans[i];
		if (fan_active[i] > 2)
		{
			std::size_t var_idx = num_faces + num_usage_vars;
			edge_usage_status[&fan] = var_idx;
			++num_usage_vars;
		}
	}
    double coeff_data_fitting = Method::lambda_data_fitting;
	double coeff_hight = total_points * Method::lambda_model_height / bbox_zheight;
	coeff_hight /= num_faces;
	double coeff_complexity = total_points * Method::lambda_model_complexity / double(fans.size());

	typedef Variable<double> Variable;
	typedef LinearExpression<double> Objective;
	typedef LinearConstraint<double> Constraint;
	typedef LinearProgram<double> LinearProgram;
	Objective obj;

	std::size_t num_sharp_edges = 0;
	for (std::size_t i = 0; i < fans.size(); ++i)
	{
		const FaceStar& fan = fans[i];
		if (fan.size() == 4 && fan_active[i] == 4) // fully intact edge only
		{
			std::size_t var_idx = num_faces + num_usage_vars + num_sharp_edges;
			edge_sharp_status[&fan] = var_idx;
			// accumulate model complexity term
			obj.add_coefficient(var_idx, coeff_complexity);
			++num_sharp_edges;
		}
	}

	FOR_EACH_FACET(Map, model_, it)
	{
		Map::Facet* f = it;
		if (is_excluded(f))
			continue;
		std::size_t var_idx = facet_indices[f];
		vec3 c(0, 0, 0);
		int degree = 0;
		FacetHalfedgeCirculator cir(f);
		for (; !cir->end(); ++cir)
		{
			c += cir->halfedge()->vertex()->point();
			++degree;
		}
		c /= degree;
		double avg_z = bbox_zmax - c.z;
		// accumulate data fitting term
		double num = facet_attrib_supporting_point_num_[f];

		obj.add_coefficient(var_idx, -coeff_data_fitting * num);
		// accumulate face height  term
		obj.add_coefficient(var_idx,  coeff_hight * avg_z);

		auto v_plane = facet_attrib_supporting_plane_[f];
		bool v_face = false;
		for (int i = 0; i < v.size(); ++i)
		{
			if (v_plane == v[i])
			{
				v_face = true;
				break;
			}
		}
		double uncovered_area = (facet_attrib_facet_area_[f]);
		// accumulate vertical face area term,
        double coeff_vertical_coverage = 0.1*total_points * Method::lambda_model_complexity / model_->bbox().area();
        if (v_face)
		{
			obj.add_coefficient(var_idx, coeff_vertical_coverage * uncovered_area);
		}
	}
	program_.set_objective(obj, LinearProgram::MINIMIZE);

	std::size_t total_variables = num_faces + num_usage_vars + num_sharp_edges;
	typedef LinearProgram::Variable Variable;
	for (std::size_t i = 0; i < total_variables; ++i)
	{
		program_.add_variable(Variable(Variable::BINARY));
	}

	Logger::out("    -") << "num total variables: " << total_variables << std::endl;

	//////////////////////////////////////////////////////////////////////////

	typedef LinearProgram::Constraint Constraint;

	// Add constraints: the number of faces associated with an edge must be either 2 or 0.
	// Edges degraded by excluded faces (fewer than three selectable members
	// left) stay unconstrained — the LP may select one side alone, and the
	// closure repair after solving pulls the excluded filler face back in.
	for (std::size_t i = 0; i < fans.size(); ++i)
	{
		const FaceStar& fan = fans[i];

		if (fan_active[i] > 2)
		{
			Constraint constraint;
			for (std::size_t j = 0; j < fan.size(); ++j)
			{
				MapTypes::Facet* f = fan[j]->facet();
				if (is_excluded(f))
					continue;
				constraint.add_coefficient(facet_indices[f], 1.0);
			}

			constraint.add_coefficient(edge_usage_status[&fan], -2.0);  //
			constraint.set_bounds(Constraint::FIXED, 0.0, 0.0);
			program_.add_constraint(constraint);
		}
	}

	FOR_EACH_FACET(Map, model_, it)
	{
		Map::Facet* f = it;
		if (is_excluded(f))
			continue;
		auto v_plane = facet_attrib_supporting_plane_[f];
		bool v_face = false;
		for (int i = 0; i < v.size(); ++i)
		{
			if (v_plane == v[i])
			{
				v_face = true;
				break;
			}
		}

		//exclude the added vertical planar segments
		if (!v_face && is_confident_facet[f] )
		{
			Constraint constraint;
			std::size_t var_idx1 = facet_indices[f];
			constraint.add_coefficient(var_idx1, 1.0);
			constraint.set_bounds(Constraint::FIXED, 1.0, 1.0);
			program_.add_constraint(constraint);

		}
	}

	// Add constraints: for the sharp edges (fully intact degree-4 edges only)
	double M = 1.0;
	for (std::size_t i = 0; i < fans.size(); ++i)
	{
		const FaceStar& fan = fans[i];
		if (fan.size() != 4 || fan_active[i] != 4)
			continue;

		Constraint edge_used_constraint;
		std::size_t var_edge_usage_idx = edge_usage_status[&fan];
		edge_used_constraint.add_coefficient(var_edge_usage_idx, 1.0);

		std::size_t var_edge_sharp_idx = edge_sharp_status[&fan];
		edge_used_constraint.add_coefficient(var_edge_sharp_idx, -1.0);

		edge_used_constraint.set_bounds(Constraint::LOWER, 0.0, 0.0);
		program_.add_constraint(edge_used_constraint);

		for (std::size_t j = 0; j < fan.size(); ++j)
		{
			MapTypes::Facet* f1 = fan[j]->facet();
			Plane3d* plane1 = facet_attrib_supporting_plane_[f1];
			std::size_t fid1 = facet_indices[f1];
			for (std::size_t k = j + 1; k < fan.size(); ++k)
			{
				MapTypes::Facet* f2 = fan[k]->facet();
				Plane3d* plane2 = facet_attrib_supporting_plane_[f2];
				std::size_t fid2 = facet_indices[f2];

				if (plane1 != plane2)
				{
					Constraint edge_sharp_constraint;
					edge_sharp_constraint.add_coefficient(var_edge_sharp_idx, 1.0);
					edge_sharp_constraint.add_coefficient(fid1, -M);
					edge_sharp_constraint.add_coefficient(fid2, -M);
					edge_sharp_constraint.add_coefficient(var_edge_usage_idx, -M);
					edge_sharp_constraint.set_bounds(Constraint::LOWER, 1.0 - 3.0 * M, 0.0);
					program_.add_constraint(edge_sharp_constraint);
				}
			}
		}
	}

	// Add constraints: single-roof
	std::vector<std::vector<Map::Facet*> > multiple_roofs;
	{
		StageScope stage("09a2_multi_roofs");
		multiple_roofs = find_multi_roofs(model_, footprint, v);
	}
	std::set<std::vector<int>> multi_roof_set;
	for (std::size_t i = 0; i < multiple_roofs.size(); ++i)
	{
		const std::vector<Map::Facet*>& roofs = multiple_roofs[i];
		if (is_excluded(roofs[0]))
			continue; // group owned by an excluded face: wouldn't exist without it
		std::set<int> temp;
		std::vector<int> face_index;
		for (std::size_t j = 0; j < roofs.size(); ++j)
		{
			Map::Facet* f = roofs[j];
			if (is_excluded(f))
				continue;
			std::size_t fid = facet_indices[f];
			temp.insert(fid);
		}
		if (temp.size() < 2)
			continue; // would pin a single remaining face
		face_index.insert(face_index.end(), temp.begin(), temp.end());
		multi_roof_set.insert(face_index);

	}

	for (std::set<std::vector<int>>::iterator iter = multi_roof_set.begin(); iter != multi_roof_set.end(); ++iter)
	{
		const std::vector<int>& roofs = *iter;
		Constraint constraint;
		for (std::size_t j = 0; j < roofs.size(); ++j)
		{
			std::size_t fid = roofs[j];
			constraint.add_coefficient(fid, 1.0);
		}
		// "at most one" instead of "exactly one": a stack whose faces also
		// meet along one arrangement edge must satisfy the even "2 or 0"
		// fan rule, and forcing exactly one (odd) there makes the program
		// infeasible. This configuration is data-dependent and shows up
		// with the unsupported-face exclusion, so defuse the whole class.
		constraint.set_bounds(Constraint::UPPER, 0.0, 1.0);
		program_.add_constraint(constraint);
	}


	std::vector<std::size_t> vertical_group;
	FOR_EACH_FACET(Map, model_, it)
	{
		Map::Facet* f = it;
		if (is_excluded(f))
			continue;
		bool vertical_facet = false;
		for (std::size_t k = 0; k < v.size(); k++)
		{
			if (facet_attrib_supporting_plane_[it] == v[k])
			{
				vertical_facet = true;
				break;
			}
		}
		if (vertical_facet)
		{
			int count_boundary = 0;
			FacetHalfedgeCirculator circulator(f);
			for (; !circulator->end(); ++circulator)
			{
				Map::Halfedge* f_edge = circulator->halfedge();
				if (is_bound[f_edge])
					count_boundary++;
			}
			if (count_boundary > 2)
			{
				std::size_t v_index = facet_indices[f];
				auto fplg = f->to_polygon();
				vertical_group.push_back(v_index);
			}
		}
	}
 // the boundary vertical facets cannot be selected
	for (int l = 0; l < vertical_group.size(); ++l)
	{
		Constraint v_constraint;
		std::size_t fid = vertical_group[l];
		v_constraint.add_coefficient(fid, 1.0);
		v_constraint.set_bounds(Constraint::FIXED, 0.0, 0.0);
		program_.add_constraint(v_constraint);
	}
	//////////////////////////////////////////////////////////////////////////
	} // stage 09a_lp_build

	// Near-optimal guidance for the selection program: MIP time on the
	// marked hard cases swings from seconds to minutes between runs, and
	// proving the last fraction of optimality never changed the geometry.
	// A gap plus a time cap keeps the whole reconstruction predictable;
	// the defaults are loose enough that easy models still solve exactly.
	{
		double lp_time = 20.0, lp_gap = 1e-3;
		if (const char* env = std::getenv("CITY3D_LP_TIME"))
			lp_time = std::atof(env);
		if (const char* env = std::getenv("CITY3D_LP_GAP"))
			lp_gap = std::atof(env);
		program_.set_solver_guidance(lp_time, lp_gap, false);
	}

	LinearProgramSolver solver;
	bool success;
	{
		StageScope stage("09b_lp_solve");
		success = solver.solve(&program_, solver_name);
	}
	if (success) {
		StageScope stage("09c_lp_postprocess");
		// mark results
		const std::vector<double>& X = solver.get_result();
		std::set<Map::Facet*> keep;
		FOR_EACH_FACET(Map, model_, it)
		{
			Map::Facet* f = it;
			if (is_excluded(f))
				continue;
			if (static_cast<int>(std::round(X[facet_indices[f]])) == 1)
				keep.insert(f);
		}

		// Closure repair: fans degraded by excluded faces were left
		// unconstrained, so the selection can end up with an odd number of
		// kept faces around such an edge. Restore the even "2 or 0" parity
		// with a minimum-area closure program over the repair-eligible
		// excluded faces instead of a local greedy walk: pulling back the
		// same-plane neighbour fan by fan floods whole extrapolated sheets
		// onto the selection (case7: roof area x10, 80 z-levels instead of
		// 27). The closure program decides globally which faces complete
		// the surface with the least added area. z-outside faces stay
		// pinned to 0 (geometrically absurd), and fans that no filler can
		// close pay a large slack and may remain odd.
		const std::size_t lp_kept_count = keep.size();
		auto greedy_flood = [&]() {
			bool changed = true;
			while (changed)
			{
				changed = false;
				for (std::size_t i = 0; i < fans.size(); ++i)
				{
					const FaceStar& fan = fans[i];
					if (fan.size() < 2)
						continue;
					std::size_t kept_count = 0;
					std::set<Plane3d*> kept_planes;
					for (std::size_t j = 0; j < fan.size(); ++j)
					{
						Map::Facet* f = fan[j]->facet();
						if (f != nil && keep.count(f) != 0)
						{
							++kept_count;
							kept_planes.insert(facet_attrib_supporting_plane_[f]);
						}
					}
					if (kept_count % 2 == 0)
						continue;
					Map::Facet* add = nil;
					for (std::size_t j = 0; j < fan.size() && add == nil; ++j)
					{
						Map::Facet* f = fan[j]->facet();
						if (f == nil || keep.count(f) == 0)
							continue;
						Map::Facet* g = fan[j]->opposite()->facet();
						if (g != nil && keep.count(g) == 0 && repair_eligible(g) &&
							kept_planes.count(facet_attrib_supporting_plane_[g]) != 0)
							add = g; // mesh-edge partner on an endorsed plane
					}
					if (add == nil)
					{
						for (std::size_t j = 0; j < fan.size() && add == nil; ++j)
						{
							Map::Facet* g = fan[j]->facet();
							if (g != nil && keep.count(g) == 0 && repair_eligible(g) &&
								kept_planes.count(facet_attrib_supporting_plane_[g]) != 0)
								add = g;
						}
					}
					if (add != nil)
					{
						keep.insert(add);
						changed = true;
					}
				}
			}
		};
		{
			// reverse index: face -> indices of the fans containing it
			std::map<Map::Facet*, std::vector<std::size_t> > face_fans;
			for (std::size_t i = 0; i < fans.size(); ++i)
			{
				const FaceStar& fan = fans[i];
				for (std::size_t j = 0; j < fan.size(); ++j)
				{
					Map::Facet* f = fan[j]->facet();
					if (f != nil)
						face_fans[f].push_back(i);
				}
			}

			// relevant fans and filler candidates, breadth-limited around the
			// selection: the minimal closure never strays far from the
			// surface boundary, and unbounded candidate sets make the
			// parity program itself a hard MIP (case1: 25k candidates ran
			// into the solver time limit). Fans on the last ring stay in the
			// constraints and may keep an odd remainder via their slack.
			std::size_t closure_rings = 3;
			if (const char* env = std::getenv("CITY3D_CLOSURE_RINGS"))
				closure_rings = std::atoi(env);
			std::set<std::size_t> relevant_fans;
			std::set<Map::Facet*> fillers;
			std::vector<std::pair<Map::Facet*, std::size_t> > worklist;
			for (std::set<Map::Facet*>::const_iterator it = keep.begin(); it != keep.end(); ++it)
				worklist.push_back(std::make_pair(*it, std::size_t(0)));
			while (!worklist.empty())
			{
				Map::Facet* f = worklist.back().first;
				std::size_t depth = worklist.back().second;
				worklist.pop_back();
				const std::vector<std::size_t>& f_fans = face_fans[f];
				for (std::size_t k = 0; k < f_fans.size(); ++k)
				{
					std::size_t i = f_fans[k];
					relevant_fans.insert(i);
					if (depth >= closure_rings)
						continue;
					const FaceStar& fan = fans[i];
					for (std::size_t j = 0; j < fan.size(); ++j)
					{
						Map::Facet* g = fan[j]->facet();
						if (g == nil || keep.count(g) != 0 || !repair_eligible(g) || fillers.count(g) != 0)
							continue;
						fillers.insert(g);
						worklist.push_back(std::make_pair(g, depth + 1));
					}
				}
			}

			bool closed = false;
			if (!fillers.empty())
			{
				MapFacetAttribute<double> facet_area(model_, Method::facet_attrib_facet_area);
				typedef LinearProgram<double> ClosureProgram;
				ClosureProgram closure;
				std::map<Map::Facet*, std::size_t> filler_idx;
				ClosureProgram::Objective obj;
				for (std::set<Map::Facet*>::const_iterator it = fillers.begin(); it != fillers.end(); ++it)
				{
					std::size_t idx = filler_idx.size();
					filler_idx[*it] = idx;
					// area plus a small per-face term: among closures of
					// equal area prefer the one with fewer faces
					obj.add_coefficient(idx, facet_area[*it] + 0.01);
				}
				const std::size_t num_filler_vars = filler_idx.size();

				// one u/s pair per fan keeps the variable indices simple;
				// fans outside the relevant set get free, unused slots
				for (std::size_t i = 0; i < fans.size(); ++i)
				{
					if (relevant_fans.count(i) != 0)
						obj.add_coefficient(num_filler_vars + 2 * i + 1, 1e4); // penalize odd remainders
				}
				closure.set_objective(obj, ClosureProgram::MINIMIZE);
				// near-optimal guidance: proving optimality on this parity
				// model is what stalled hard cases for minutes; a small gap
				// plus a time cap returns the incumbent quickly
				double closure_time = 20.0, closure_gap = 0.02;
				if (const char* env = std::getenv("CITY3D_CLOSURE_TIME"))
					closure_time = std::atof(env);
				if (const char* env = std::getenv("CITY3D_CLOSURE_GAP"))
					closure_gap = std::atof(env);
				closure.set_solver_guidance(closure_time, closure_gap, true);

				for (std::set<std::size_t>::const_iterator it = relevant_fans.begin(); it != relevant_fans.end(); ++it)
				{
					std::size_t i = *it;
					const FaceStar& fan = fans[i];
					std::size_t kept_count = 0;
					for (std::size_t j = 0; j < fan.size(); ++j)
					{
						Map::Facet* f = fan[j]->facet();
						if (f != nil && keep.count(f) != 0)
							++kept_count;
					}
					// kept_count + (fillers) + s - 2u = 0: u counts the face
					// pairs meeting at this edge, s the odd remainder when
					// no filler set can close the fan
					ClosureProgram::Constraint constraint;
					for (std::size_t j = 0; j < fan.size(); ++j)
					{
						Map::Facet* g = fan[j]->facet();
						if (g != nil && filler_idx.count(g) != 0)
							constraint.add_coefficient(filler_idx[g], 1.0);
					}
					constraint.add_coefficient(num_filler_vars + 2 * i + 1, 1.0);
					constraint.add_coefficient(num_filler_vars + 2 * i, -2.0);
					constraint.set_bounds(ClosureProgram::Constraint::FIXED,
						-static_cast<double>(kept_count), -static_cast<double>(kept_count));
					closure.add_constraint(constraint);
				}

				for (std::size_t i = 0; i < num_filler_vars; ++i)
					closure.add_variable(ClosureProgram::Variable(ClosureProgram::Variable::BINARY));
				for (std::size_t i = 0; i < fans.size(); ++i)
				{
					ClosureProgram::Variable u(ClosureProgram::Variable::INTEGER);
					u.set_bounds(ClosureProgram::Variable::DOUBLE, 0.0, 4.0);
					closure.add_variable(u);
					ClosureProgram::Variable s(ClosureProgram::Variable::CONTINUOUS);
					s.set_bounds(ClosureProgram::Variable::DOUBLE, 0.0, 1.0);
					closure.add_variable(s);
				}

				LinearProgramSolver closure_solver;
				StageScope stage("09c2_closure_solve");
				if (closure_solver.solve(&closure, solver_name))
				{
					const std::vector<double>& Y = closure_solver.get_result();
					std::size_t added = 0;
					for (std::map<Map::Facet*, std::size_t>::const_iterator it = filler_idx.begin();
						it != filler_idx.end(); ++it)
					{
						if (static_cast<int>(std::round(Y[it->second])) == 1)
						{
							keep.insert(it->first);
							++added;
						}
					}
					closed = true;
					Logger::out("-") << "closure program: " << fillers.size() << " candidates, "
						<< added << " fillers added" << std::endl;
				}
				else
					Logger::warn("-") << "closure program failed, falling back to greedy repair" << std::endl;
			}
			if (!closed)
				greedy_flood();
		}
		Logger::out("-") << keep.size() << " faces kept after closure repair (LP selected " << lp_kept_count << ")" << std::endl;

		std::vector<Map::Facet*> to_delete;
		FOR_EACH_FACET(Map, model_, it)
		{
			Map::Facet* f = it;
			if (keep.count(f) == 0)
				to_delete.push_back(f);
		}

		MapEditor editor(model_);
		for (std::size_t i = 0; i < to_delete.size(); ++i)
		{
			Map::Facet* f = to_delete[i];
			editor.erase_facet(f->halfedge());
		}

		//////////////////////////////////////////////////////////////////////////

		// mark the sharp edges
		MapHalfedgeAttribute<bool> edge_is_sharp(model_, "SharpEdge");
		FOR_EACH_EDGE(Map, model_, it) edge_is_sharp[it] = false;

		for (std::size_t i = 0; i < fans.size(); ++i)
		{
			const FaceStar& fan = fans[i];
			if (fan.size() != 4 || fan_active[i] != 4)
				continue;

			std::size_t idx_sharp_var = edge_sharp_status[&fan];
			if (static_cast<int>(X[idx_sharp_var]) == 1)
			{
				for (std::size_t j = 0; j < fan.size(); ++j)
				{
					Map::Halfedge* e = fan[j];
					Map::Facet* f = e->facet();
					if (f)
					{ // some faces may be deleted
						std::size_t fid = facet_indices[f];
						if (static_cast<int>(X[fid]) == 1)
						{
							edge_is_sharp[e] = true;
							break;
						}
					}
				}
			}
		}

		MapFacetNormal normal(model_);
		FOR_EACH_FACET(Map, model_, it)
		{
			normal[it] = facet_attrib_supporting_plane_[it]->normal();
		}

		// Cosmetic, last: glue adjacent facets that share one supporting
		// plane back into single polygons. The selection can keep many
		// coplanar pieces (closure-repair filler); merging leaves the
		// geometry identical but drops the face count back to a clean
		// level. Runs after the sharp-edge marking so its halfedge
		// iteration never touches merged-away edges; the merged edge itself
		// joins two same-plane faces and can never have been sharp.
		{
			// disabled by default: the splice loses area on large merge
			// chains (case7: roof area halved) — root cause not yet found;
			// enable with CITY3D_MERGE_OUTPUT=1 to experiment
			const char* env = std::getenv("CITY3D_MERGE_OUTPUT");
			if (env == nil || std::atoi(env) == 0) {
				// skip
			} else {
			std::size_t num_merged = 0;
			bool progress = true;
			while (progress)
			{
				progress = false;
				std::vector<Map::Halfedge*> candidates;
				FOR_EACH_HALFEDGE(Map, model_, it)
				{
					Map::Halfedge* h = it;
					if (h > h->opposite())
						continue; // collect one halfedge per edge only
					Map::Facet* f1 = h->facet();
					Map::Facet* f2 = h->opposite()->facet();
					if (f1 == nil || f2 == nil || f1 == f2)
						continue;
					if (facet_attrib_supporting_plane_[f1] != facet_attrib_supporting_plane_[f2])
						continue;
					candidates.push_back(h);
				}
				for (std::size_t k = 0; k < candidates.size(); ++k)
				{
					Map::Halfedge* h = candidates[k];
					Map::Facet* f1 = h->facet();
					Map::Facet* f2 = h->opposite()->facet();
					if (f1 == nil || f2 == nil || f1 == f2)
						continue; // consumed by an earlier merge in this round
					if (facet_attrib_supporting_plane_[f1] != facet_attrib_supporting_plane_[f2])
						continue;
					if (editor.merge_facets_along_edge(h))
					{
						++num_merged;
						progress = true;
					}
				}
			}
			if (num_merged > 0)
				Logger::out("-") << num_merged << " coplanar facet pairs merged for output" << std::endl;
			}
		}
	}

	facet_attrib_supporting_vertex_group_.unbind();
	facet_attrib_supporting_point_num_.unbind();
	facet_attrib_facet_area_.unbind();
	facet_attrib_covered_area_.unbind();
	facet_attrib_supporting_plane_.unbind();

    if (success)
        return true;
    else {
        Logger::warn("-") << "solving the binary program failed." << std::endl;
        return false;
    }
}
