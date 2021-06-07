


SET(USE_SOLVER "CPLEX") # "QSOPT" "GUROBI", "CPLEX"



SET(QSOPT_HOME_EXPORT QSOPT_HOME=/opt/qsopt/)
#SET(QSOPT_HOME_EXPORT QSOPT_HOME=C:/Users/Timo/Documents/Uni/BA/Code/qsopt)

SET(GUROBI_HOME_EXPORT GUROBI_HOME=/opt/gurobi901/linux64)

SET(CPLEX_HOME_EXPORT CPLEX_HOME=/opt/cplex/CPLEX_Studio126/cplex/)

######################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(ddcolors)

SET(CMAKE_CXX_STANDARD 20)

SET(SOURCE_FILES main.cpp Graph.cpp Graph.h DecisionDiagram.cpp DecisionDiagram.h DDColors.cpp DDColors.h scip_interface.cpp scip_interface.h)
#set(SOURCE_FILES ${SOURCE_FILES} asc2bin.c asc2bin.c)


####################
SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR arm)
####################


ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE_FILES})


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")#-O3 oder -gp
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wshadow -W  -Wmissing-declarations -Wpointer-arith -Wundef -Wcast-qual -Wcast-align -Wwrite-strings")


LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


#SET(EXACTCOLORS_ROOT_DIR /cygdrive/c/Users/Timo/Documents/Uni/BA/clone_exactcolors/exactcolors)
FIND_PACKAGE(EXACTCOLORS)
IF(EXACTCOLORS_FOUND)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${EXACTCOLORS_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${EXACTCOLORS_LIBRARIES})
ELSE()
    MESSAGE(STATUS "No precompiled version of exactcolors found, cloning and building from github")
    IF(USE_SOLVER STREQUAL "QSOPT")
        SET(SOLVER_EXPORT ${QSOPT_HOME_EXPORT})
    ENDIF()
    IF(USE_SOLVER STREQUAL "GUROBI")
        SET(SOLVER_EXPORT ${GUROBI_HOME_EXPORT})
    ENDIF()
        IF(USE_SOLVER STREQUAL "CPLEX")
        SET(SOLVER_EXPORT ${CPLEX_HOME_EXPORT})
    ENDIF()
    INCLUDE(ExternalProject)
    MESSAGE(STATUS ${SOLVER_EXPORT})
    ExternalProject_Add(EXACTCOLORS
            GIT_REPOSITORY https://github.com/heldstephan/exactcolors
            DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/exactcolors
            SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/exactcolors
            PATCH_COMMAND cmake -E echo "Skipping patch step."
            UPDATE_COMMAND cmake -E echo "Skipping update step."
            CONFIGURE_COMMAND cmake -E echo "Skipping configure step."
            BUILD_COMMAND make clean && export ${SOLVER_EXPORT} && cmake -E echo "HI" && make "libexactcolor.a" && cd mwis_sewell && make "libsewell.a"
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


#INCLUDE(CHECKSYMBOLEXISTS)
#CHECK_SYMBOL_EXISTS(COLORlp_init_env ${EXACTCOLORS_LIBRARIES} FOUND_FUNCTION)
#MESSAGE(STATUS ${FOUND_FUNCTION})


IF(USE_SOLVER STREQUAL "QSOPT")
    SET(QSOPT_ROOT_DIR /opt/qsopt)
#    SET(QSOPT_ROOT_DIR C:/Users/Timo/Documents/Uni/BA/Code/qsopt)
    FIND_PACKAGE(QSOPT REQUIRED)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${QSOPT_INCLUDE_DIR})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${QSOPT_LIBRARY})
    MESSAGE(STATUS "Warning: QSOPT is only an LP solver and will treat the variables of the IP as continuous.")
    MESSAGE(STATUS "Warning: QSOPT This can find the chromatic number too but is not guaranteed to do so!")
ENDIF()

IF(USE_SOLVER STREQUAL "GUROBI")
    SET(GUROBI_ROOT_DIR /opt/gurobi901/linux64)
    FIND_PACKAGE(GUROBI REQUIRED)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${GUROBI_INCLUDE_DIR})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${GUROBI_LIBRARIES})
ENDIF()

IF(USE_SOLVER STREQUAL "CPLEX")
    SET(CPLEX_ROOT_DIR /opt/cplex/CPLEX_Studio126)
    FIND_PACKAGE(CPLEX REQUIRED)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${CPLEX_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CPLEX_LIBRARIES})
ENDIF()
#message make sure you also compiled exactcolors with the right solver

IF(TRUE)
    SET(SCIP_ROOT_DIR /home/brand/local_opt/scip-6.0.2/installdir/)
    FIND_PACKAGE(SCIP REQUIRED)
    TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${SCIP_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${SCIP_LIBRARIES})
ENDIF()


#export ILOG_LICENSE_FILE=/opt/ILOG/CPLEX_Studio_AcademicResearch122/licenses/access.ilm
#export CPLEX_HOME=/opt/cplex/CPLEX_Studio126/cplex/

#unset CPLEX_HOME
#export GUROBI_HOME=/opt/gurobi901/linux64

#export QSOPT_HOME=/opt/qsopt/


