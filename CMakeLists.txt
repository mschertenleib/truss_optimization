cmake_minimum_required(VERSION 3.19)

project(olc_code_jam LANGUAGES CXX)


set(TARGET_WEB OFF)
if(${CMAKE_SYSTEM_NAME} MATCHES "Emscripten")
    set(TARGET_WEB ON)
endif()


set(CPM_DOWNLOAD_VERSION 0.38.2)

if (CPM_SOURCE_CACHE)
    set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
elseif (DEFINED ENV{CPM_SOURCE_CACHE})
    set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
else ()
    set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
endif ()

if (NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
    message(STATUS "Downloading CPM.cmake to ${CPM_DOWNLOAD_LOCATION}")
    file(DOWNLOAD
        https://github.com/TheLartians/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
        ${CPM_DOWNLOAD_LOCATION}
    )
endif ()

include(${CPM_DOWNLOAD_LOCATION})


CPMAddPackage(
    GITLAB_REPOSITORY libeigen/eigen
    GIT_TAG c487a4fe
    DOWNLOAD_ONLY YES
)
add_compile_definitions(
    #EIGEN_NO_AUTOMATIC_RESIZING
    EIGEN_RUNTIME_NO_MALLOC
    EIGEN_INITIALIZE_MATRICES_BY_NAN
)

if (NOT TARGET_WEB)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    CPMAddPackage(
        GITHUB_REPOSITORY glfw/glfw
        GIT_TAG 3.4
    )
endif ()

CPMAddPackage(
    GITHUB_REPOSITORY ocornut/imgui
    GIT_TAG v1.91.9
    DOWNLOAD_ONLY YES
)
add_library(imgui STATIC)
target_sources(imgui PRIVATE
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui.h
    ${imgui_SOURCE_DIR}/imgui_internal.h
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui SYSTEM PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui PRIVATE glfw)


add_executable(${CMAKE_PROJECT_NAME})
target_sources(${CMAKE_PROJECT_NAME} PRIVATE
    src/main.cpp
    src/application.hpp src/application.cpp
    src/analysis.hpp src/analysis.cpp
    src/unique_resource.hpp
    src/vec.hpp
)
target_include_directories(${CMAKE_PROJECT_NAME} SYSTEM PRIVATE ${eigen_SOURCE_DIR})
target_compile_features(${CMAKE_PROJECT_NAME} PRIVATE cxx_std_23)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE imgui)


if (TARGET_WEB)
    set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/web)
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${EMSCRIPTEN_ROOT_PATH}/sysroot/include)
    target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE
        -Os
        -sDISABLE_EXCEPTION_CATCHING=0
    )
    target_link_options(${CMAKE_PROJECT_NAME} PRIVATE
        -Os
        -sWASM=1
        -sUSE_WEBGL2=1
        -sUSE_GLFW=3
        -sFULL_ES3=1
        -sASSERTIONS=1
        -sDISABLE_EXCEPTION_CATCHING=0
        -sFORCE_FILESYSTEM=1
        -sALLOW_MEMORY_GROWTH=1
        --emrun
        --shell-file ${CMAKE_SOURCE_DIR}/src/shell.html
    )
else ()
    target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE glfw)
endif ()


set(clang_warnings
    -Wall
    -Wextra
    -Wshadow
    -Wnon-virtual-dtor
    -Wold-style-cast
    -Wcast-align
    -Wunused
    -Woverloaded-virtual
    -Wpedantic
    -Wconversion
    -Wsign-conversion
    -Wnull-dereference
    -Wdouble-promotion
)
set(gcc_warnings
    ${clang_warnings}
    -Wmisleading-indentation
    -Wduplicated-cond
    -Wduplicated-branches
    -Wlogical-op
    -Wswitch-enum
)
if (CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE ${clang_warnings})
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE ${gcc_warnings})
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE /W4)
else ()
    message(WARNING "No warnings set for compiler ${CMAKE_CXX_COMPILER_ID}")
endif ()


if ((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
    target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE $<$<CONFIG:Debug>:-fsanitize=undefined>)
    #target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE $<$<CONFIG:Debug>:-fsanitize=address>)
    target_link_options(${CMAKE_PROJECT_NAME} PRIVATE $<$<CONFIG:Debug>:-fsanitize=undefined>)
    #target_link_options(${CMAKE_PROJECT_NAME} PRIVATE $<$<CONFIG:Debug>:-fsanitize=address>)
endif ()


message(CHECK_START "Checking for IPO support")
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported_result OUTPUT ipo_supported_output)
if (ipo_supported_result)
    message(CHECK_PASS "supported")
    set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION_MIN_SIZE_REL TRUE)
else ()
    message(CHECK_FAIL "not supported")
endif ()
