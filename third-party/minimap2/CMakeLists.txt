cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0048 NEW)
project(minimap2 VERSION 2.5.572 LANGUAGES C)

set(minimap2_RootDir    ${minimap2_SOURCE_DIR})
set(minimap2_SourceDir  ${minimap2_RootDir}/src)
set(minimap2_IncludeDir ${minimap2_RootDir}/src)

## DEPENDENCIES
# Threads
if (NOT Threads)
    find_package(Threads REQUIRED)
endif()
# ZLIB
if (NOT ZLIB_INCLUDE_DIRS OR NOT ZLIB_LIBRARIES)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(ZLIB zlib)
else()
    set(ZLIB_LDFLAGS ${ZLIB_LIBRARIES})
endif()

add_library(minimap2
    ${minimap2_SourceDir}/align.c
    ${minimap2_SourceDir}/bseq.c
    ${minimap2_SourceDir}/chain.c
    ${minimap2_SourceDir}/esterr.c
    ${minimap2_SourceDir}/format.c
    ${minimap2_SourceDir}/hit.c
    ${minimap2_SourceDir}/index.c
    ${minimap2_SourceDir}/pe.c
    ${minimap2_SourceDir}/ksw2_extd2_sse.c
    ${minimap2_SourceDir}/ksw2_exts2_sse.c
    ${minimap2_SourceDir}/ksw2_extz2_sse.c
    ${minimap2_SourceDir}/ksw2_ll_sse.c
    ${minimap2_SourceDir}/kalloc.c
    ${minimap2_SourceDir}/kthread.c
    ${minimap2_SourceDir}/map.c
    ${minimap2_SourceDir}/misc.c
    ${minimap2_SourceDir}/sdust.c
    ${minimap2_SourceDir}/sketch.c
)

# definitions
target_compile_definitions(minimap2
    PUBLIC
    HAVE_KALLOC=1
)
# includes
target_include_directories(minimap2
    PUBLIC
    ${minimap2_IncludeDir}
    ${ZLIB_INCLUDE_DIRS}
)
# link dependencies
target_link_libraries(minimap2
    PUBLIC
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

set(minimap2_INCLUDE_DIRS
    ${minimap2_IncludeDir}
    ${ZLIB_INCLUDE_DIRS}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)
set(minimap2_LIBRARIES
    $<TARGET_FILE:minimap2>
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    CACHE INTERNAL
    "${PROJECT_NAME}: Libraries"
    FORCE
)

# CFLAGS=		-g -Wall -O2 -Wc++-compat
# CPPFLAGS=	-DHAVE_KALLOC
# INCLUDES=
# OBJS=		kthread.o kalloc.o misc.o bseq.o sketch.o sdust.o index.o chain.o align.o hit.o map.o format.o pe.o ksw2_ll_sse.o
# PROG=		minimap2
# PROG_EXTRA=	sdust minimap2-lite
# LIBS=		-lm -lz -lpthread
