SUBDIRS(libseq)
SUBDIRS(mxml)
SUBDIRS(mxml-utils)

SET(LIB_SOURCES
    Alignment.cpp
    CLIUtils.cpp
    Utils.cpp
    ReferenceSequence.cpp
    ResultsExporter.cpp
)

include_directories(libseq mxml mxml-utils)

ADD_LIBRARY(virulignlib ${LIB_SOURCES})
ADD_EXECUTABLE(virulign Virulign.cpp)
TARGET_LINK_LIBRARIES(virulign virulignlib seq mxml mxml-utils)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

install(TARGETS virulign DESTINATION bin)
