# This module finds the scip header and libraries.
#
# User can give SCIP_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
#  SCIP_FOUND              - Set to false, or undefined, if scip isn't found.
#  SCIP_INCLUDE_DIRS       - include directory
#  SCIP_LIBRARIES          - library files

# use given hint directory or look in parent/scip folder
IF(NOT DEFINED SCIP_ROOT_DIR)
    GET_FILENAME_COMPONENT(DIR_ONE_ABOVE ../ ABSOLUTE)
    SET(SCIP_ROOT_DIR ${DIR_ONE_ABOVE}/scip)
ENDIF()


FIND_PATH(SCIP_INCLUDE_DIR scip.h scip
        PATH_SUFFIXES include
        PATHS ${SCIP_ROOT_DIR}
        )
MESSAGE(STATUS "SCIP Include Dir: ${SCIP_INCLUDE_DIR}")

FIND_LIBRARY(SCIP_LIBRARY
        NAMES scip libscip.a
        PATH_SUFFIXES lib
        PATHS ${SCIP_ROOT_DIR}
        )
MESSAGE(STATUS "SCIP Library: ${SCIP_LIBRARY}")

# check whether required things have been found and set SCIP_FOUND accordingly
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCIP DEFAULT_MSG SCIP_INCLUDE_DIR SCIP_LIBRARY)

# if package has been found, collect include dirs and libraries
IF(SCIP_FOUND)
    SET(SCIP_INCLUDE_DIRS ${SCIP_INCLUDE_DIR})
    SET(SCIP_LIBRARIES ${SCIP_LIBRARY})

    IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        SET(SCIP_LIBRARIES "${SCIP_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
ENDIF(SCIP_FOUND)

# Advance the variables not needed anymore
MARK_AS_ADVANCED(SCIP_INCLUDE_DIR SCIP_LIBRARY)



