SET(SOURCES
    progress_bar.cpp
    stream_helper.cpp
    color.cpp
    setting.cpp
    terminal_size.cpp
    display_width.cpp
    termcolor.cpp
    font_style.cpp
    progress_type.cpp
    cursor_control.cpp
)
    
ADD_LIBRARY(indicators ${SOURCES})
