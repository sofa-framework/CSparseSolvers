cmake_minimum_required(VERSION 3.24)
project(CSparseSolvers VERSION 1.0 LANGUAGES CXX)

find_package(Sofa.Config REQUIRED)
sofa_find_package(Sofa.Component.LinearSolver.Direct REQUIRED)

# List all files
set(CSPARSESOLVERS_SRC_DIR src/${PROJECT_NAME})
set(HEADER_FILES
    ${CSPARSESOLVERS_SRC_DIR}/config.h.in

    ${CSPARSESOLVERS_SRC_DIR}/SparseLUSolver.h
    ${CSPARSESOLVERS_SRC_DIR}/SparseLUSolver.inl
    ${CSPARSESOLVERS_SRC_DIR}/SparseCholeskySolver.h
    ${CSPARSESOLVERS_SRC_DIR}/SparseCholeskySolver.inl
)
set(SOURCE_FILES
    ${CSPARSESOLVERS_SRC_DIR}/init.cpp

    ${CSPARSESOLVERS_SRC_DIR}/SparseLUSolver.cpp
    ${CSPARSESOLVERS_SRC_DIR}/SparseCholeskySolver.cpp
)
set(README_FILES
    README.md
)

add_subdirectory(extlibs/csparse)

# Create the plugin library.
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Component.LinearSolver.Direct)
target_link_libraries(${PROJECT_NAME} PUBLIC csparse)

# Create package Config, Version & Target files.
# Deploy the headers, resources, scenes & examples.
# Set the plugin 'relocatable' if built within SOFA.
# --> see SofaMacros.cmake
sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
)
