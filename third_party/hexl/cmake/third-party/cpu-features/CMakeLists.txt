# Copyright (C) 2020 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

# Add cpu-features directly to our build. This defines the cpu_features target.
# Cache needed to override option in cpu-features
hexl_cache_variable(BUILD_SHARED_LIBS)
hexl_cache_variable(BUILD_PIC)
hexl_cache_variable(BUILD_TESTING)

set(BUILD_PIC ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

add_subdirectory(${HEXL_ROOT_DIR}/../cpu_features
                 ${CMAKE_CURRENT_BINARY_DIR}/cpu-features-build
                 EXCLUDE_FROM_ALL)

# Prevent variables from leaking
unset(BUILD_PIC CACHE)
unset(BUILD_SHARED_LIBS CACHE)
unset(BUILD_TESTING CACHE)

hexl_uncache_variable(BUILD_SHARED_LIBS)
hexl_uncache_variable(BUILD_PIC)
hexl_uncache_variable(BUILD_TESTING)
