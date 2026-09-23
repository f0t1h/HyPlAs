cmake_minimum_required(VERSION 3.20)

# ---------------------------------------------------------------------------
# Environment-aware defaults (must precede project()/find_package()).
# ---------------------------------------------------------------------------
# Make the active conda/virtualenv discoverable for find_package() and default
# the install prefix to it, matching the behavior of the old configure script.
if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
    set(_hyplas_env_prefix "$ENV{CONDA_PREFIX}")
elseif(DEFINED ENV{VIRTUAL_ENV})
    set(_hyplas_env_prefix "$ENV{VIRTUAL_ENV}")
endif()

if(DEFINED _hyplas_env_prefix AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${_hyplas_env_prefix}" CACHE PATH
        "Install prefix (defaulted to active conda/virtual env)" FORCE)
endif()

project(hyplas VERSION 0.0.1 LANGUAGES CXX)

# Default to an optimized build when the caller hasn't chosen one.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Options (mirror ./configure --enable-*).
# ---------------------------------------------------------------------------
option(HYPLAS_SANITIZERS "Build with AddressSanitizer + UBSan" OFF)
option(HYPLAS_STATIC     "Prefer static libgcc/libstdc++ linking" OFF)

# ---------------------------------------------------------------------------
# Dependencies.
# ---------------------------------------------------------------------------
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

include(FetchContent)
FetchContent_Declare(gtl
    GIT_REPOSITORY https://github.com/greg7mdp/gtl.git
    GIT_TAG        v1.2.0)

# shrn: subprocess execution + output verification.
FetchContent_Declare(shrn
    GIT_REPOSITORY https://github.com/f0t1h/shrn.git
    GIT_TAG        v0.11.0)
FetchContent_MakeAvailable(gtl shrn)

# ---------------------------------------------------------------------------
# Core library: everything except the CLI entrypoint, shared by the binary and
# the test harness so both build the pipeline logic exactly once.
# ---------------------------------------------------------------------------
add_library(hyplas_core STATIC
    src/log.cpp
    src/gfa.cpp
    src/fastx.cpp
    src/pipeline.cpp
)
target_include_directories(hyplas_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(hyplas_core PUBLIC -Wall -Wextra -Wpedantic)
target_link_libraries(hyplas_core PUBLIC ZLIB::ZLIB Threads::Threads gtl shrn::shrn)

if(HYPLAS_SANITIZERS)
    target_compile_options(hyplas_core PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(hyplas_core PUBLIC -fsanitize=address,undefined)
endif()

if(HYPLAS_STATIC)
    target_link_options(hyplas_core PUBLIC -static-libgcc -static-libstdc++)
endif()

# ---------------------------------------------------------------------------
# Executable + install.
# ---------------------------------------------------------------------------
add_executable(hyplas src/hyplas_main.cpp)
target_link_libraries(hyplas PRIVATE hyplas_core)

install(TARGETS hyplas RUNTIME DESTINATION bin)

# ---------------------------------------------------------------------------
# Tests (ctest). BUILD_TESTING defaults ON via CTest module.
# ---------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTING)
    add_executable(hyplas_tests tests/test_transforms.cpp)
    target_link_libraries(hyplas_tests PRIVATE hyplas_core)
    add_test(NAME transforms COMMAND hyplas_tests)
endif()
