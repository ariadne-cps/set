if(NOT TARGET solving)
    ensure_submodule_initialized(solving)
    add_subdirectory(solving)
endif()
