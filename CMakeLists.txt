cmake_minimum_required(VERSION 3.9)

project(DOMAIN_DECOMPOSITION LANGUAGES CXX C)

# Create variables
###############################################################################################
set(DOMAIN_DECOMPOSITION_SOURCES ${PROJECT_NAME}_sources)
set(DOMAIN_DECOMPOSITION_HEADERS ${PROJECT_NAME}_headers)
set(DOMAIN_DECOMPOSITION_INCLUDE ${PROJECT_NAME}_include)
set(DOMAIN_DECOMPOSITION_LINKED_LIBRARIES ${PROJECT_NAME}_linked_libraries)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

# Add dependencies
###############################################################################################

## GeDiM
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gedim/GeDiM)
if (NOT TARGET GeDiM::GeDiM)
	find_package(GeDiM CONFIG REQUIRED)
	list(APPEND ${DOMAIN_DECOMPOSITION_LINKED_LIBRARIES} PRIVATE GeDiM::GeDiM)
else()
	get_target_property(GEDIM_LIBRARY_HEADERS GeDiM::GeDiM PUBLIC_HEADER)
	get_target_property(GEDIM_LIBRARY_NAME GeDiM::GeDiM EXPORT_NAME)
	get_target_property(GEDIM_LIBRARY_LINKED_LIBRARIES GeDiM::GeDiM LINK_LIBRARIES)
	get_target_property(GEDIM_LIBRARY_INCLUDE_DIRECTORIES GeDiM::GeDiM INCLUDE_DIRECTORIES)
	get_target_property(GEDIM_LIBRARY_SOURCES GeDiM::GeDiM SOURCES)

	list(APPEND ${DOMAIN_DECOMPOSITION_LINKED_LIBRARIES} PRIVATE ${GEDIM_LIBRARY_NAME} ${GEDIM_LIBRARY_LINKED_LIBRARIES})
	list(APPEND ${DOMAIN_DECOMPOSITION_HEADERS} ${GEDIM_LIBRARY_HEADERS})
	list(APPEND ${DOMAIN_DECOMPOSITION_INCLUDE} ${GEDIM_LIBRARY_INCLUDE_DIRECTORIES})
endif()

## Blas
find_package(BLAS REQUIRED)

## Lapack
find_package(LAPACK REQUIRED)

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)

## VTK
if(${ENABLE_VTK})
  find_package(VTK CONFIG REQUIRED)
endif()

## triangle
find_package(triangle CONFIG REQUIRED)

## tetgen
find_package(tetgen CONFIG REQUIRED)

### MPI
find_package(MPI REQUIRED) # Find the MPI library
list(APPEND ${DOMAIN_DECOMPOSITION_LINKED_LIBRARIES} PUBLIC MPI::MPI_CXX)

### Openmp
find_package(OpenMP) # Find the OpenMP library
list(APPEND ${DOMAIN_DECOMPOSITION_LINKED_LIBRARIES} PUBLIC OpenMP::OpenMP_CXX)

# Insert Sources
###############################################################################################
add_subdirectory(src)
list(APPEND ${DOMAIN_DECOMPOSITION_SOURCES} ${domain_decomposition_sources})
list(APPEND ${DOMAIN_DECOMPOSITION_HEADERS} ${domain_decomposition_headers})
list(APPEND ${DOMAIN_DECOMPOSITION_INCLUDE} ${domain_decomposition_includes})

# Create library test executable
###############################################################################################
add_executable(${PROJECT_NAME}
	main.cpp
	${${DOMAIN_DECOMPOSITION_SOURCES}}
    ${${DOMAIN_DECOMPOSITION_HEADERS}})

target_include_directories(${PROJECT_NAME}
    PRIVATE ${${DOMAIN_DECOMPOSITION_INCLUDE}})

target_link_libraries(${PROJECT_NAME} ${${DOMAIN_DECOMPOSITION_LINKED_LIBRARIES}})
