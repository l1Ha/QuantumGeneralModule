cmake_minimum_required(VERSION 3.12)
project(GeneralModule LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# 编译优化与告警参数
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    add_compile_options(-O3 -fPIC -Wall -Wextra)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options(-O3 -fPIC)
endif()

# 统一模块输出目录与输出路径
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# 源文件发现
file(GLOB SRC_FILES "src/*.f90")

# 使用 OBJECT 库单次编译 Fortran 源代码，彻底杜绝并行编译时的 .mod 写入竞争冲突
add_library(GeneralModule_obj OBJECT ${SRC_FILES})
set_target_properties(GeneralModule_obj PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_include_directories(GeneralModule_obj PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
)

# 生成静态库与动态共享库
add_library(GeneralModule_static STATIC $<TARGET_OBJECTS:GeneralModule_obj>)
set_target_properties(GeneralModule_static PROPERTIES OUTPUT_NAME "general_module")

add_library(GeneralModule_shared SHARED $<TARGET_OBJECTS:GeneralModule_obj>)
set_target_properties(GeneralModule_shared PROPERTIES OUTPUT_NAME "general_module")

target_include_directories(GeneralModule_static PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>
)
target_include_directories(GeneralModule_shared PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>
)

# 单元测试构建配置 (41 个单元测试套件)
include(CTest)
if(BUILD_TESTING)
    file(GLOB TEST_SOURCES "tests/*.f90")
    foreach(TEST_SRC ${TEST_SOURCES})
        get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
        add_executable(${TEST_NAME} ${TEST_SRC})
        target_include_directories(${TEST_NAME} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${CMAKE_Fortran_MODULE_DIRECTORY}
        )
        target_link_libraries(${TEST_NAME} PRIVATE GeneralModule_static)
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    endforeach()
endif()

# 示例程序构建配置 (36 个前沿物理应用示例)
option(BUILD_EXAMPLES "Build physics examples" ON)
if(BUILD_EXAMPLES)
    file(GLOB EXAMPLE_SOURCES "examples/*.f90")
    foreach(EX_SRC ${EXAMPLE_SOURCES})
        get_filename_component(EX_NAME ${EX_SRC} NAME_WE)
        add_executable(${EX_NAME} ${EX_SRC})
        target_include_directories(${EX_NAME} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${CMAKE_Fortran_MODULE_DIRECTORY}
        )
        target_link_libraries(${EX_NAME} PRIVATE GeneralModule_static)
    endforeach()
endif()

# 安装目标规范
install(TARGETS GeneralModule_static GeneralModule_shared
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include
    FILES_MATCHING PATTERN "*.mod"
)
