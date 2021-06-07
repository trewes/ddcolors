# This module finds the gurobi header and libraries.
#
# User can give GUROBI_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
#  GUROBI_FOUND              - Set to false, or undefined, if gurobi isn't found.
#  GUROBI_INCLUDE_DIRS       - include directory
#  GUROBI_LIBRARIES          - library files

# use given hint directory or look in parent/gurobi folder
IF(NOT DEFINED GUROBI_ROOT_DIR)
    GET_FILENAME_COMPONENT(DIR_ONE_ABOVE ../ ABSOLUTE)
    SET(GUROBI_ROOT_DIR ${DIR_ONE_ABOVE}/gurobi)
ENDIF()


FIND_PATH(GUROBI_INCLUDE_DIR gurobi_c.h gurobi_c++.h
        PATH_SUFFIXES include
        PATHS ${GUROBI_ROOT_DIR}
        )
MESSAGE(STATUS "GUROBI Include Dir: ${GUROBI_INCLUDE_DIR}")

FIND_LIBRARY(GUROBI_LIBRARY
        NAMES gurobi libgurobi90.so
        PATH_SUFFIXES lib
        PATHS ${GUROBI_ROOT_DIR}
        )
MESSAGE(STATUS "GUROBI Library: ${GUROBI_LIBRARY}")

FIND_LIBRARY(GUROBI_LIBRARY_LIGHT
        NAMES gurobi libgurobi90_light.so
        PATH_SUFFIXES lib
        PATHS ${GUROBI_ROOT_DIR}
        )
MESSAGE(STATUS "GUROBI light Library: ${GUROBI_LIBRARY_LIGHT}")



# check whether required things have been found and set GUROBI_FOUND accordingly
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GUROBI DEFAULT_MSG GUROBI_INCLUDE_DIR GUROBI_LIBRARY)

# if package has been found, collect include dirs and libraries
IF(GUROBI_FOUND)
    SET(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})
    SET(GUROBI_LIBRARIES ${GUROBI_LIBRARY} ${GUROBI_LIBRARY_LIGHT} ${GUROBI_LIBRARY_5_2})

    IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        SET(GUROBI_LIBRARIES "${GUROBI_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
ENDIF(GUROBI_FOUND)

# Advance the variables not needed anymore
MARK_AS_ADVANCED(GUROBI_INCLUDE_DIR GUROBI_LIBRARY)



