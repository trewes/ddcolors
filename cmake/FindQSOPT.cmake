# This module finds the qsopt header and libraries.
#
# User can give QSOPT_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
#  QSOPT_FOUND              - Set to false, or undefined, if qsopt isn't found.
#  QSOPT_INCLUDE_DIRS       - include directory
#  QSOPT_LIBRARIES          - library files

# use given hint directory or look in parent/qsopt folder
IF(NOT DEFINED QSOPT_ROOT_DIR)
    GET_FILENAME_COMPONENT(DIR_ONE_ABOVE ../ ABSOLUTE)
    SET(QSOPT_ROOT_DIR ${DIR_ONE_ABOVE}/qsopt)
ENDIF()


FIND_PATH(QSOPT_INCLUDE_DIR qsopt.h
        PATH_SUFFIXES include
        PATHS ${QSOPT_ROOT_DIR}
        )
MESSAGE(STATUS "QSOPT Include Dir: ${QSOPT_INCLUDE_DIR}")

FIND_LIBRARY(QSOPT_LIBRARY
        NAMES qsopt qsopt.a
        PATH_SUFFIXES lib
        PATHS ${QSOPT_ROOT_DIR}
        )



# check whether required things have been found and set QSOPT_FOUND accordingly
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(QSOPT DEFAULT_MSG QSOPT_INCLUDE_DIR QSOPT_LIBRARY)

# if package has been found, collect include dirs and libraries
IF(QSOPT_FOUND)
    SET(QSOPT_INCLUDE_DIRS ${QSOPT_INCLUDE_DIR})
    SET(QSOPT_LIBRARIES ${QSOPT_LIBRARY})

    IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        SET(QSOPT_LIBRARIES "${QSOPT_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
ENDIF(QSOPT_FOUND)

# Advance the variables not needed anymore
MARK_AS_ADVANCED(QSOPT_INCLUDE_DIR QSOPT_LIBRARY)




## This module finds qsopt.
#
#get_filename_component(DIR_ONE_ABOVE ../ ABSOLUTE)
##message(STATUS ${DIR_ONE_ABOVE})
#set(QSOPT_DIR ${DIR_ONE_ABOVE}/qsopt)
##message(STATUS ${QSOPT_DIR})
## The places to look for the tinyxml2 folders
#set(FIND_Qsopt_PATHS ${QSOPT_DIR} /opt/qsopt)
#
## The location of the include folder (and thus the header files)
## find_path uses the paths we defined above as places to look
## Saves the location of the header files in a variable called TINYXML2_INCLUDE_DIR
#find_path(Qsopt_INCLUDE_DIR qsopt.h   # The variable to store the path in and the name of the header files
#        PATH_SUFFIXES include               # The folder name containing the header files
#        PATHS ${FIND_Qsopt_PATHS})       # Where to look (defined above)
#
## The location of the lib folder (and thus the .a file)
## find_library uses the paths we defined above as places to look
## Saves the location of the .a file in a variable called TINYXML2_LIBRARY
#find_library(Qsopt_LIBRARY               # The variable to store where it found the .a files
#        NAMES qsopt qsopt.a                     # The name of the .a file (without the extension and without the 'lib')
#        PATH_SUFFIXES lib                   # The folder the .a file is in
#        PATHS ${FIND_Qsopt_PATHS})               # Where to look (defined above)
