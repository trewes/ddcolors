# This module finds the exactcolors header and libraries.
#
# User can give EXACTCOLORS_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
#  EXACTCOLORS_FOUND              - Set to false, or undefined, if exactcolors isn't found.
#  EXACTCOLORS_INCLUDE_DIRS       - include directory
#  EXACTCOLORS_LIBRARIES          - library files

# use given hint directory or look in parent/exactcolors folder
IF(NOT DEFINED EXACTCOLORS_ROOT_DIR)
    SET(EXACTCOLORS_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../exactcolors)
ENDIF()

FIND_PATH(EXACTCOLORS_INCLUDE_DIR lp.h
          PATH_SUFFIXES include
          PATHS ${EXACTCOLORS_ROOT_DIR}
          )
MESSAGE(STATUS "EXACTCOLORS Include Dir: ${EXACTCOLORS_INCLUDE_DIR}")

FIND_PATH(EXACTCOLORS_MWIS_SEWELL_DIR mwss.h
          PATH_SUFFIXES include
          HINTS ${EXACTCOLORS_ROOT_DIR}/mwis_sewell
          PATHS ${EXACTCOLORS_ROOT_DIR}
          )
MESSAGE(STATUS "MWIS SEWELL Include Dir: ${EXACTCOLORS_MWIS_SEWELL_DIR}")

FIND_LIBRARY(EXACTCOLORS_LIBRARY
             NAMES exactcolor libexactcolor.a
             PATH_SUFFIXES lib
             PATHS ${EXACTCOLORS_ROOT_DIR}
             )
MESSAGE(STATUS "EXACTCOLORS Library: ${EXACTCOLORS_LIBRARY}")

FIND_LIBRARY(EXACTCOLORS_SEWELL_LIBRARY
             NAMES sewell libsewell.a
             PATH_SUFFIXES lib
             HINTS ${EXACTCOLORS_ROOT_DIR}/mwis_sewell
             PATHS ${EXACTCOLORS_ROOT_DIR}
             )
MESSAGE(STATUS "SEWELL Library: ${EXACTCOLORS_SEWELL_LIBRARY}")

# check whether required things have been found and set EXACTCOLORS_FOUND accordingly
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(EXACTCOLORS DEFAULT_MSG
                                  EXACTCOLORS_INCLUDE_DIR EXACTCOLORS_MWIS_SEWELL_DIR EXACTCOLORS_LIBRARY EXACTCOLORS_SEWELL_LIBRARY)

# if package has been found, collect include dirs and libraries
IF(EXACTCOLORS_FOUND)
    SET(EXACTCOLORS_INCLUDE_DIRS ${EXACTCOLORS_INCLUDE_DIR} ${EXACTCOLORS_MWIS_SEWELL_DIR})
    SET(EXACTCOLORS_LIBRARIES ${EXACTCOLORS_LIBRARY} ${EXACTCOLORS_SEWELL_LIBRARY})

    IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        SET(EXACTCOLORS_LIBRARIES "${EXACTCOLORS_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
ENDIF(EXACTCOLORS_FOUND)

# Advance the variables not needed anymore
MARK_AS_ADVANCED(EXACTCOLORS_INCLUDE_DIR EXACTCOLORS_MWIS_SEWELL_DIR EXACTCOLORS_LIBRARY EXACTCOLORS_SEWELL_LIBRARY)

