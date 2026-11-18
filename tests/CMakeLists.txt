include(CTest)
set(UNIT_TESTS
    test_affine_sets
    test_constraint_set
    test_validated_constrained_image_set
    test_function_set_sensitivity
    test_constrained_image_set
    test_curve
    test_list_set
    test_multifunction_image
    test_polytope
    test_support_polygon
    test_spre22_outer_projection
)
foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} PRIVATE ariadne-set set ariadne-testing)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}> --concise --logtonull)
    set_property(TEST ${TEST} PROPERTY LABELS "set")
    if(WIN32)
        set_tests_properties(${TEST} PROPERTIES ENVIRONMENT_MODIFICATION
            "PATH=path_list_prepend:$<TARGET_FILE_DIR:ariadne-set>")
    endif()
endforeach()
add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS})
# The real SPRE22 replay must never hold a CI runner indefinitely if a
# future simplex regression reintroduces cycling.
set_tests_properties(test_spre22_outer_projection PROPERTIES TIMEOUT 120)
