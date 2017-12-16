add_executable(pbmm2
    ${PBMM2_SourceDir}/main.cpp
    ${PBMM2_SourceDir}/Workflow.cpp
    ${PBMM2_SourceDir}/Settings.cpp
)
set_target_properties(pbmm2 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

set(PBMM2_INCLUDE_DIRS
    ${Boost_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${minimap2_INCLUDE_DIRS}
    ${PacBioBAM_INCLUDE_DIRS}
    ${pbcopper_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

set(PBMM2_LIBRARIES
    ${ZLIB_LDFLAGS}
    ${HTSLIB_LIBRARIES}
    ${minimap2_LIBRARIES}
    ${PacBioBAM_LIBRARIES}
    ${pbcopper_LIBRARIES}
    CACHE INTERNAL
    "${PROJECT_NAME}: Libraries"
    FORCE
)

target_include_directories(pbmm2
    PUBLIC
    ${PBMM2_INCLUDE_DIRS}
)
target_link_libraries(pbmm2
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
    ${PBMM2_LIBRARIES}
)

install(TARGETS pbmm2 RUNTIME DESTINATION bin)

add_executable(rscheck
    ${PBMM2_SourceDir}/rscheck.cpp
)
set_target_properties(rscheck PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
target_include_directories(rscheck
    PUBLIC
    ${PBMM2_INCLUDE_DIRS}
)
target_link_libraries(rscheck
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
    ${PBMM2_LIBRARIES}
)
install(TARGETS rscheck RUNTIME DESTINATION bin)
