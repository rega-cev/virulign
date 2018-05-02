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
