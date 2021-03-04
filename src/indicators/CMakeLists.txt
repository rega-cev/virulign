SET(SOURCES
    block_progress_bar.cpp
    color.cpp
    cursor_control.cpp
    cursor_movement.cpp
    display_width.cpp
    dynamic_progress.cpp
    font_style.cpp
    indeterminate_progress_bar.cpp
    multi_progress.cpp
    progress_bar.cpp
    progress_spinner.cpp
    progress_type.cpp
    setting.cpp
    stream_helper.cpp
    termcolor.cpp
    terminal_size.cpp
)
    
ADD_LIBRARY(indicators ${SOURCES})
