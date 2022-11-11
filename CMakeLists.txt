######################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(ddcolors)

SET(CMAKE_CXX_STANDARD 20)

SET(SOURCE_FILES main.cpp Graph.cpp Graph.h DecisionDiagram.cpp DecisionDiagram.h DDColors.cpp DDColors.h)


####################
SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR arm)
####################


ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE_FILES})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")



LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")




FIND_PACKAGE(EXACTCOLORS)
IF(EXACTCOLORS_FOUND)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${EXACTCOLORS_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${EXACTCOLORS_LIBRARIES})
ELSE()
    MESSAGE(STATUS "No precompiled version of exactcolors found, cloning and building from github")
    IF(USE_SOLVER STREQUAL "CPLEX")
        SET(SOLVER_EXPORT ${CPLEX_HOME_EXPORT})

    ENDIF()
    INCLUDE(ExternalProject)
    SET(ENV{CPLEX_HOME} ${CPLEX_ROOT_DIR}/cplex)
    ExternalProject_Add(EXACTCOLORS
                        GIT_REPOSITORY https://github.com/heldstephan/exactcolors
                        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/exactcolors
                        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/exactcolors
                        PATCH_COMMAND cmake -E echo "Skipping patch step."
                        UPDATE_COMMAND cmake -E echo "Skipping update step."
                        CONFIGURE_COMMAND cmake -E echo "Skipping configure step."
                        BUILD_COMMAND make clean && export "CPLEX_HOME=${CPLEX_ROOT_DIR}/cplex" && make "libexactcolor.a" && cd mwis_sewell && make "libsewell.a"
                        BUILD_IN_SOURCE 1
                        #            EXCLUDE_FROM_ALL TRUE
                        INSTALL_COMMAND cmake -E echo "Skipping install step."
                        )
    ADD_DEPENDENCIES(${PROJECT_NAME} EXACTCOLORS)
    ExternalProject_Get_Property(EXACTCOLORS SOURCE_DIR)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${SOURCE_DIR})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${SOURCE_DIR}/libexactcolor.a)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${SOURCE_DIR}/mwis_sewell)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${SOURCE_DIR}/mwis_sewell/libsewell.a)
ENDIF()

SET(USE_SOLVER "CPLEX")

IF(USE_SOLVER STREQUAL "CPLEX")
    FIND_PACKAGE(CPLEX REQUIRED)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${CPLEX_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CPLEX_LIBRARIES})

    SET(CPLEX_HOME_EXPORT CPLEX_HOME=${CPLEX_ROOT_DIR}/cplex/)
    SET(ENV{CPLEX_HOME} ${CPLEX_ROOT_DIR}/cplex/)
ENDIF()


UNSET(CPLEX_ROOT_DIR CACHE)
