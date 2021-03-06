cmake_minimum_required(VERSION 3.9)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Aboria/cmake"
        ${CMAKE_MODULE_PATH})

# Boost
find_package(Boost 1.50.0 COMPONENTS python REQUIRED)
list(APPEND LIBRARIES ${Boost_LIBRARIES})
list(APPEND INCLUDES ${Boost_INCLUDE_DIRS})

# VTK
find_package(VTK REQUIRED)
if (VTK_FOUND)
    add_definitions(-DHAVE_VTK)
endif(VTK_FOUND)
list(APPEND LIBRARIES ${VTK_LIBRARIES})
list(APPEND INCLUDES ${VTK_INCLUDE_DIRS})

# Eigen
find_package(Eigen3 REQUIRED)
list(APPEND INCLUDES ${EIGEN3_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11")


# Aboria
set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})
list(APPEND INCLUDES Aboria/src)
list(APPEND INCLUDES Aboria/third-party)

include_directories(src ${INCLUDES})


#add_executable(fixed_gradient fixed_gradient.cpp)
#add_executable(fixed_gradient2 fixed_gradient2.cpp)
#add_executable(varying_chemo varying_chemo_new.cpp)
#add_executable(no_filopodia_varying_chemo no_filopodia_varying_chemo.cpp)
#add_executable(linear_chemo linear_chemo.cpp)
#add_executable(linear_chemo_random linear_chemo_random.cpp)
#add_executable(linear_chemo_random_copy linear_chemo_random_copy.cpp)
#add_executable(linear_move_up_the_gradient linear_move_up_the_gradient.cpp)
#add_executable(varying_chemo_random varying_chemo_random.cpp)
#add_executable(neat_fixed_gradient neat_fixed_gradient.cpp)
#add_executable(growing_domain growing_domain.cpp)
#add_executable(full_dynamics full_dynamics.cpp)
#add_executable(full_dynamics_fixed_matrix_size full_dynamics_fixed_matrix_size.cpp)
#add_executable(full_dynamics_multiple_direction full_dynamics_multiple_direction.cpp)
add_executable(leaders-followers leaders-followers.cpp)
#add_executable(old_multiple_direction old_multiple_direction.cpp)
#add_executable(chemical_trail_leaders_followers chemical_trail_leaders_followers.cpp)
#add_executable(new_version_fixed_gradient new_version_fixed_gradient.cpp)
#add_executable(no_cells no_cells.cpp)


#target_link_libraries(fixed_gradient ${LIBRARIES})
#target_link_libraries(fixed_gradient2 ${LIBRARIES})
#target_link_libraries(varying_chemo ${LIBRARIES})
#target_link_libraries(no_filopodia_varying_chemo ${LIBRARIES})
#target_link_libraries(linear_chemo ${LIBRARIES})
#target_link_libraries(linear_chemo_random ${LIBRARIES})
#target_link_libraries(linear_move_up_the_gradient ${LIBRARIES})
#target_link_libraries(varying_chemo_random ${LIBRARIES})
#target_link_libraries(neat_fixed_gradient ${LIBRARIES})
#target_link_libraries(growing_domain ${LIBRARIES})
#target_link_libraries(full_dynamics ${LIBRARIES})
#target_link_libraries(full_dynamics_fixed_matrix_size ${LIBRARIES})
#target_link_libraries(full_dynamics_multiple_direction ${LIBRARIES})
target_link_libraries(leaders-followers ${LIBRARIES})
#target_link_libraries(old_multiple_direction ${LIBRARIES})
#target_link_libraries(chemical_trail_leaders_followers ${LIBRARIES})
#target_link_libraries(new_version_fixed_gradient ${LIBRARIES})
#target_link_libraries(no_cells ${LIBRARIES})
