if(NOT TARGET linear_programming)
    ariadne_add_submodule(linear-programming)
endif()

if(NOT TARGET constraint_solving)
    ariadne_add_submodule(constraint-solving)
endif()
