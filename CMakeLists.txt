cmake_minimum_required(VERSION 3.12)
project(GeneralModule LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# 编译优化参数
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    add_compile_options(-O3 -fPIC -Wall -Wextra)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options(-O3 -fPIC)
endif()

file(GLOB_RECURSE SRC_FILES "src/*.f90")

# 生成静态库与动态共享库
add_library(GeneralModule_static STATIC ${SRC_FILES})
set_target_properties(GeneralModule_static PROPERTIES OUTPUT_NAME "general_module")

add_library(GeneralModule_shared SHARED ${SRC_FILES})
set_target_properties(GeneralModule_shared PROPERTIES OUTPUT_NAME "general_module")

target_include_directories(GeneralModule_static PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)
target_include_directories(GeneralModule_shared PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

# 编译测试用例
enable_testing()

file(GLOB TEST_SOURCES "tests/*.f90")
foreach(TEST_SRC ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} PRIVATE GeneralModule_static)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
