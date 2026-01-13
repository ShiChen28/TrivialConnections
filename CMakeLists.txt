cmake_minimum_required(VERSION 3.15)
project(connection VERSION 1.0 LANGUAGES CXX)

# -----------------------------
# 编译模式设置
# -----------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# -----------------------------
# SuiteSparse 手动安装配置
# -----------------------------
set(SuiteSparse_ROOT_DIR "E:/Library/SR_ROOT/suitesparse-metis-for-windows/build/install")
set(SuiteSparse_INCLUDE_DIR "${SuiteSparse_ROOT_DIR}/include/suitesparse")
set(SuiteSparse_DIR "${SuiteSparse_ROOT_DIR}/lib/cmake/suitesparse-7.5.1")

# 查找 SuiteSparse
find_package(SuiteSparse CONFIG REQUIRED)

# -----------------------------
# 查找 OpenGL / GLUT / GLEW
# -----------------------------
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)

# -----------------------------
# 收集源文件
# -----------------------------
file(GLOB SOURCES
    "src/*.cpp"
    "src/*.cxx"
    "src/*.h"
    "src/*.hpp"
)

# -----------------------------
# 创建可执行文件
# -----------------------------
add_executable(${PROJECT_NAME} ${SOURCES})

# -----------------------------
# 配置头文件包含路径
# -----------------------------
target_include_directories(${PROJECT_NAME} PUBLIC 
    "include"
    "src"
    ${SuiteSparse_INCLUDE_DIR}
)

# -----------------------------
# 链接库 (修复 LNK2019 和之前的错误)
# -----------------------------
target_link_libraries(${PROJECT_NAME} PRIVATE 
    SuiteSparse::umfpack          # 必须链接这个来解决 umfpack_zl 错误
    SuiteSparse::cholmod
    SuiteSparse::spqr
    SuiteSparse::amd
    SuiteSparse::colamd
    SuiteSparse::suitesparseconfig
    OpenGL::GL
    GLUT::GLUT
    GLEW::GLEW
)

# -----------------------------
# 针对 Windows 的宏定义 (修复 C1189)
# -----------------------------
target_compile_definitions(${PROJECT_NAME} PUBLIC 
    _USE_MATH_DEFINES 
    NOMINMAX
    # 注意：千万不要在这里定义 UF_long=long long
)

# -----------------------------
# 编译选项
# -----------------------------
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE
        /W3 /MP
        $<$<CONFIG:RELEASE>:/O2>
        $<$<CONFIG:DEBUG>:/Zi /Od /RTC1>
    )
endif()

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)