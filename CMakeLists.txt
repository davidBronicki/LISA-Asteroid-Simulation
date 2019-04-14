cmake_minimum_required(VERSION 3.0)
project(LISA_Simulation CXX)

include_directories("${PROJECT_SOURCE_DIR}/Secular_Simulation_Code")
find_package(Threads)

file(GLOB Secular_Simulation_Code_Pile
	"Secular_Simulation_Code/*.cpp")

add_executable(sec_sim ${Secular_Simulation_Code_Pile})
set_property(TARGET sec_sim PROPERTY CXX_STANDARD 17)
target_link_libraries(sec_sim ${CMAKE_THREAD_LIBS_INIT})
