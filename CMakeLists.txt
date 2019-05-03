cmake_minimum_required(VERSION 3.0)
project(LISA_Simulation CXX)


find_package(Threads)


include_directories("${PROJECT_SOURCE_DIR}/Secular_Simulation_Code")
file(GLOB Secular_Simulation_Code_Pile
	"Secular_Simulation_Code/*.cpp")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Secular_Output_Directory")

add_executable(secular_simulation ${Secular_Simulation_Code_Pile})
set_property(TARGET secular_simulation PROPERTY CXX_STANDARD 17)
target_link_libraries(secular_simulation ${CMAKE_THREAD_LIBS_INIT})



include_directories("${PROJECT_SOURCE_DIR}/Bolen_Test")
file(GLOB Bolen_Test_Code_Pile
	"Bolen_Test/*.cpp")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Bolen_Test_Output_Directory")

add_executable(bolen_test ${Bolen_Test_Code_Pile})
set_property(TARGET bolen_test PROPERTY CXX_STANDARD 17)

