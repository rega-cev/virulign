SUBDIRS(libseq)
SUBDIRS(mxml)
SUBDIRS(mxml-utils)
SUBDIRS(indicators)

SET(LIB_SOURCES
    Alignment.cpp
    CLIUtils.cpp
    Utils.cpp
    ReferenceSequence.cpp
    ResultsExporter.cpp
)

include_directories(libseq mxml mxml-utils indicators)

ADD_LIBRARY(virulignlib ${LIB_SOURCES})
ADD_EXECUTABLE(virulign Virulign.cpp)
TARGET_LINK_LIBRARIES(virulign virulignlib seq mxml mxml-utils indicators)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
   # Disable warnings for all the unknown OpenMP pragmas.
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

find_package(Git)
find_package(Hg)
add_custom_target(version
  ${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/src/virulign_version.h.in
                   -D DST=${CMAKE_SOURCE_DIR}/src/virulign_version.h
                   -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                   -D HG_EXECUTABLE=${HG_EXECUTABLE}
                   -P ${CMAKE_SOURCE_DIR}/src/GenerateVersionHeader.cmake
  )
add_dependencies(virulign version)

install(TARGETS virulign DESTINATION bin)
