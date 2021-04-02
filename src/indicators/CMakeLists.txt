SET(SOURCES
    progress_bar.cpp
    stream_helper.cpp
    setting.cpp
    terminal_size.cpp
    display_width.cpp
    progress_type.cpp
    cursor_control.cpp
)
    
ADD_LIBRARY(indicators ${SOURCES})
