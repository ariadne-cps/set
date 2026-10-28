#include "pybind11.hpp"
void set_submodule(pybind11::module& module);
void set_module(pybind11::module& module) { set_submodule(module); }
