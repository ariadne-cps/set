add_library(SET_SRC OBJECT
    curve.cpp
    function_set.cpp
    multifunction_image.cpp
    affine_set.cpp
)
target_include_directories(SET_SRC PRIVATE "${PROJECT_SOURCE_DIR}/include/set")
