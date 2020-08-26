cmake_minimum_required(VERSION 3.5)
project(mock_se)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(PythonLibs 2.7)


include_directories(include
  include
  PRIVATE ${PYTHON_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
  src/Trackgraph.cpp
  src/Path_planning_utils.cpp
  src/cubic_spline.cpp
)

add_executable(talker src/mock_se_pub.cpp)
ament_target_dependencies(talker rclcpp std_msgs geometry_msgs)
add_executable(listener src/mock_se_sub.cpp)
ament_target_dependencies(listener rclcpp std_msgs geometry_msgs)
add_executable(path_planner src/path_planner.cpp)
ament_target_dependencies(path_planner rclcpp std_msgs geometry_msgs)

target_link_libraries(
  path_planner
  ${PROJECT_NAME}
  ${PYTHON_LIBRARIES}
)


install(TARGETS
  talker
  listener
  path_planner
  DESTINATION lib/${PROJECT_NAME}

)



# if(BUILD_TESTING)
#   find_package(ament_lint_auto REQUIRED)
#   # the following line skips the linter which checks for copyrights
#   # uncomment the line when a copyright and license is not present in all source files
#   #set(ament_cmake_copyright_FOUND TRUE)
#   # the following line skips cpplint (only works in a git repo)
#   # uncomment the line when this package is not in a git repo
#   #set(ament_cmake_cpplint_FOUND TRUE)
#   ament_lint_auto_find_test_dependencies()
# endif()

ament_package()
