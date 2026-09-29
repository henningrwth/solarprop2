
set(HEADERS arrows.h
            binning.h
            density.h
            filehelper.h
	    green.h
            model.h
            modelfactory.h
            parkerspiral.h
            parkerspiralhcs.h
            modifiedspiral.h
            numdiff.h
	    particle.h
	    properties.h
	    random.h
	    sde.h
            sde2d.h
            sde3d.h
            spectrum.h
	    stringtools.h
	    units.h
            solarprop.h
	    )

set(SOURCES arrows.cc
            binning.cc
	    density.cc
            filehelper.cc
	    green.cc
	    model.cc
            parkerspiral.cc
            parkerspiralhcs.cc
            modifiedspiral.cc
            numdiff.cc
	    particle.cc
            properties.cc
	    random.cc
            sde.cc
            sde2d.cc
            sde3d.cc
	    spectrum.cc
	    stringtools.cc
            units.cc
            solarprop.cc
	    )

include_directories(${GSL_INCLUDE_DIRS})

if(OpenMP_FOUND)
  include_directories(${OpenMP_CXX_INCLUDE_DIRS})
endif()

# Prepare the python module library.
python_add_library(solarprop MODULE ${SOURCES} ${HEADERS} WITH_SOABI)
target_link_libraries(solarprop PRIVATE pybind11::headers pybind11::module ${GSL_LIBRARIES})
target_include_directories(solarprop PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# When not building with pip, but cmake is invoked manually,
# create a standalone executable for debugging and profiling purposes.
if(NOT DEFINED SKBUILD)
    add_executable(sprop_test sproptest.cc ${SOURCES} ${HEADERS})
    target_link_libraries(sprop_test PRIVATE pybind11::headers pybind11::embed Python::Python ${GSL_LIBRARIES})
    target_include_directories(sprop_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    install(TARGETS sprop_test DESTINATION ${CMAKE_SOURCE_DIR})
endif()

add_subdirectory(models)

# Avoid pybind stuff when compiling standalone executable.
target_compile_definitions(solarprop PRIVATE PYTHON_BINDINGS)

# Pass the version as a define.
target_compile_definitions(solarprop PRIVATE VERSION_INFO=${PROJECT_VERSION})

if(DEFINED SKBUILD)
    # The install directory is the output (wheel) directory
    install(TARGETS solarprop DESTINATION .)
endif()
