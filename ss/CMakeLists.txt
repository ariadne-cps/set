if(NOT TARGET linear_programming)
    ensure_submodule_initialized(linear-programming)
    add_subdirectory(linear-programming)
endif()

if(NOT TARGET constraint_solving)
    ensure_submodule_initialized(constraint-solving)
    add_subdirectory(constraint-solving)
endif()
