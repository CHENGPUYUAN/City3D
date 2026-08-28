//
// Created by Jin on 19/07/2022.
//
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

#include "../model/point_set.h"
#include "../model/map.h"
#include "../model/map_io.h"
#include "../model/point_set_io.h"
#include "../method/method_global.h"
#include "../method/reconstruction.h"

#include "../basic/logger.h"

#include <iostream>


int main(int argc, char **argv)
{
    // Usage: City3D_cli <input_cloud.(ply|las|laz)> <input_footprint.(obj|geojson)> <output.obj> [pixel_size] [min_points] [solver]
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_cloud.(ply|las|laz)> <input_footprint.(obj|geojson)> <output.obj> [pixel_size] [min_points] [solver]"
                  << std::endl;
        std::cerr << "  solver: gurobi | mindoptpy | highs | scip10 | scip (default: best available)"
                  << std::endl;
        return EXIT_FAILURE;
    }
    const std::string input_cloud_file = argv[1];
    const std::string input_footprint_file = argv[2];
    const std::string output_file = argv[3];

    ///ToDo: user may need to tune these parameters for their datasets
    Method::min_points = argc > 5 ? std::atoi(argv[5]) : 40;
    Method::pixel_size = argc > 4 ? std::atof(argv[4]) : 0.15;

    // the Logger only prints 'out' messages for registered features (the GUI
    // registers "*"); do the same so progress lines show up in the CLI
    Logger::initialize();
    Logger::instance()->set_value(Logger::LOG_REGISTER_FEATURES, "*");

    // load input point cloud
    std::cout << "loading input point cloud data from file: " << input_cloud_file << std::endl;
    PointSet *pset = PointSetIO::read(input_cloud_file);
    if (!pset) {
        std::cerr << "failed to load point cloud data from file: " << input_cloud_file << std::endl;
        return EXIT_FAILURE;
    }

    // load input footprint data
    std::cout << "loading input footprint data from file: " << input_footprint_file << std::endl;
    const vec3& offset = pset->offset();
    /// ToDo: in this demo the Z coordinate of the footprint/ground is set to the min_Z of the point cloud.
    ///       This is not optimal (at least noise, outliers, and incompleteness not considered).
    ///       In practice, the Z coordinate of each building should be determined by extracting local ground planes,
    ///       or directly from available DTM or DSM data.
    Map *footprint = MapIO::read(input_footprint_file, vec3(offset.x, offset.y, -pset->bbox().z_min()));
    if (!footprint) {
        std::cerr << "failed to load footprint data from file: " << input_footprint_file << std::endl;
        return EXIT_FAILURE;
    }

    Reconstruction recon;

    // Step 1: segmentation to obtain point clouds of individual buildings
    std::cout << "segmenting individual buildings..." << std::endl;
    recon.segmentation(pset, footprint);

    // Step 2: extract planes from the point cloud of each building (for all buildings)
    std::cout << "extracting roof planes..." << std::endl;
    if (!recon.extract_roofs(pset, footprint)) {
        std::cerr << "no roofs could be extracted from the point cloud" << std::endl;
        return EXIT_FAILURE;
    }

    // Step 3: reconstruction of all the buildings in the scene
    Map *result = new Map;
    std::string solver_arg = argc > 6 ? argv[6] : "auto";
    LinearProgramSolver::SolverName solver;
    bool solver_ok = true;
    if (solver_arg == "gurobi") {
#ifdef HAS_GUROBI
        solver = LinearProgramSolver::GUROBI;
#else
        solver_ok = false;
#endif
    } else if (solver_arg == "mindoptpy") {
        solver = LinearProgramSolver::MINDOPTPY;
    } else if (solver_arg == "highs") {
#ifdef HAS_HIGHS
        solver = LinearProgramSolver::HIGHS;
#else
        solver_ok = false;
#endif
    } else if (solver_arg == "scip10") {
#ifdef HAS_SCIP10
        solver = LinearProgramSolver::SCIP10;
#else
        solver_ok = false;
#endif
    } else if (solver_arg == "scip") {
        solver = LinearProgramSolver::SCIP;
    } else {
        // auto: best available — gurobi > highs > mindoptpy > scip
        solver_arg = "auto";
#ifdef HAS_GUROBI
        solver = LinearProgramSolver::GUROBI;
#elif defined(HAS_HIGHS)
        solver = LinearProgramSolver::HIGHS;
#elif defined(HAS_MINDOPTPY)
        solver = LinearProgramSolver::MINDOPTPY;
#else
        solver = LinearProgramSolver::SCIP;
#endif
    }
    if (!solver_ok) {
        std::cerr << "solver '" << solver_arg << "' is not available in this build" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "reconstructing the buildings (using the " << solver_arg << " solver)..." << std::endl;
    bool status = recon.reconstruct(pset, footprint, result, solver);

    if (status && result->size_of_facets() > 0)
    {
        if (MapIO::save(output_file, result))
        {
            std::cout << "reconstruction result saved to file: " << output_file << std::endl;
            return EXIT_SUCCESS;
        } else
            std::cerr << "failed to save reconstruction result to file: " << output_file << std::endl;
    } else
        std::cerr << "reconstruction failed" << std::endl;

    delete pset;
    delete footprint;
    delete result;

    return EXIT_FAILURE;
}
