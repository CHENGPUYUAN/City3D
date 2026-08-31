get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(math_HEADERS
        box.h
        line.h
        math_types.h
        matrix.h
        plane.h
        polygon2d.h
        principal_axes.h
        semi_definite_symmetric_eigen.h
        vecg.h
        linear_program.h
        linear_program_solver.h
        )

set(math_SOURCES
        math_types.cpp
        polygon2d.cpp
        principal_axes.cpp
        semi_definite_symmetric_eigen.cpp
        linear_program_solver.cpp
        linear_program_solver_SCIP.cpp
        linear_program_solver_GUROBI.cpp
        linear_program_solver_MINDOPTPY.cpp
        linear_program_solver_HIGHS.cpp
        linear_program_solver_CUOPT.cpp
        )

add_library(${PROJECT_NAME} STATIC ${math_SOURCES} ${math_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "City3D")

include(../cmake/FindGUROBI.cmake)
if (GUROBI_FOUND)
    # the native Gurobi libraries need a license file; without one GRBEnv()
    # throws at runtime, so fall back to the mindoptpy backend
    set(_gurobi_lic "$ENV{GRB_LICENSE_FILE}")
    if (NOT EXISTS "${_gurobi_lic}")
        set(_gurobi_lic "$ENV{HOME}/gurobi.lic")
    endif()
    if (NOT EXISTS "${_gurobi_lic}")
        set(_gurobi_lic "$ENV{GUROBI_HOME}/gurobi.lic")
    endif()
    if (EXISTS "${_gurobi_lic}")
        message(STATUS "GUROBI_INCLUDE_DIRS: ${GUROBI_INCLUDE_DIRS}")
        message(STATUS "GUROBI_LIBRARIES: ${GUROBI_LIBRARIES}")
        target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_GUROBI)
        target_include_directories(${PROJECT_NAME} PUBLIC ${GUROBI_INCLUDE_DIRS})
        target_link_libraries(${PROJECT_NAME} basic ${GUROBI_LIBRARIES})
    else()
        message(STATUS "GUROBI found but no license file; using the mindoptpy backend instead")
        set(GUROBI_FOUND FALSE)
    endif()
endif()

# System SCIP (apt scipoptsuite 10.x) replaces the bundled 8.0.3 when present:
# same C API, THREADSAFE build with TBB-backed parallel presolve (PaPILO) and
# SoPlex. Falling back to the bundled static build keeps the tree self-contained.
find_package(SCIP QUIET CONFIG)
if (TARGET libscip)
    message(STATUS "SCIP backend: system ${SCIP_VERSION} (replaces bundled 8.0.3)")
    set_target_properties(libscip PROPERTIES IMPORTED_GLOBAL TRUE)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_SYSTEM_SCIP)
    target_link_libraries(${PROJECT_NAME} libscip)
else()
    message(STATUS "SCIP backend: bundled 8.0.3")
    target_include_directories(${PROJECT_NAME} PUBLIC
            ${POLYFIT_scip_DIR}
            ${POLYFIT_scip_DIR}/scip
            )
    target_link_libraries(${PROJECT_NAME} 3rd_scip 3rd_soplex)
endif()

# HiGHS: local install in ~/.local (built from source, no license needed).
# Preferred over mindoptpy — in-process library, no subprocess/license-server hop.
find_package(HIGHS QUIET)
if (NOT HIGHS_FOUND)
    find_package(highs QUIET)
endif()
if (HIGHS_FOUND OR TARGET highs::highs)
    message(STATUS "HIGHS backend enabled (${HIGHS_VERSION})")
    target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_HIGHS)
    target_link_libraries(${PROJECT_NAME} basic highs::highs)
endif()

if (NOT GUROBI_FOUND)
    # mindoptpy runs off the local MindOpt stand-alone SDK (pip's wheel embeds
    # an expired license); the community license has no size limit. The backend
    # shells out to `uv run --with-editable <sdk>`.
    find_program(UV_EXECUTABLE uv)
    if (UV_EXECUTABLE AND EXISTS "$ENV{HOME}/mindopt/2.3.0/linux64-x86/lib/python")
        message(STATUS "mindoptpy backend enabled (via ${UV_EXECUTABLE} and ~/mindopt/2.3.0)")
        target_compile_definitions(${PROJECT_NAME} PUBLIC HAS_MINDOPTPY)
    endif()
endif()

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _CRT_SECURE_NO_WARNINGS _CRT_SECURE_NO_DEPRECATE)
endif()