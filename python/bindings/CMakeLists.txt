include_directories(${Python3_INCLUDE_DIRS})
get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-set-submodule-obj OBJECT set_submodule.cpp)
target_include_directories(pyariadne-set-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    $<TARGET_PROPERTY:pyariadne-solving,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-set-submodule-obj PRIVATE
    pyariadne-python-common
    pybind11::module
    set
)
set_target_properties(pyariadne-set-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
    target_compile_options(pyariadne-set-submodule-obj PRIVATE /bigobj /wd4702)
endif()

add_library(pyariadne-set INTERFACE)
target_sources(pyariadne-set INTERFACE $<TARGET_OBJECTS:pyariadne-set-submodule-obj>)
target_include_directories(pyariadne-set INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-solving,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-set INTERFACE pyariadne-solving)
ariadne_add_python_module(NAME set MODULE_SOURCE set_module.cpp)
