set(C_TOOLS
        teslib
        tesrecons)

set(CPP_TOOLS
        gennoisespec
        )

foreach(TOOL ${C_TOOLS})
    add_executable(${TOOL} "${TOOL}/${TOOL}.c")
endforeach()

foreach(TOOL ${CPP_TOOLS})
    add_executable(${TOOL} "${TOOL}/${TOOL}.cpp")
endforeach()

foreach(TOOL ${C_TOOLS};${CPP_TOOLS})
    target_include_directories(${TOOL} PUBLIC
            ${EXT_INCLUDE_DIRS}
            ${SIMPUT_INCLUDE_DIR}
            ${SIXTE_INCLUDE_DIR})
    target_link_directories(${TOOL} PUBLIC
            ${SIMPUT_LIBRARY_DIR}
            ${SIXTE_LIBRARY_DIR})
    target_link_libraries(${TOOL} PUBLIC
            compiler_flags
            ${EXT_LIBRARIES}
            ${SIMPUT_LIBRARIES}
            ${SIXTE_LIBRARIES}
	    sirena)

    install(TARGETS ${TOOL} DESTINATION bin)
    install(FILES ${TOOL}/${TOOL}.par DESTINATION share/sirena/pfiles)

    ### Apple Adjustments for missing boost headers ###
    if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        target_compile_definitions(${TOOL} PUBLIC unary_function=__unary_function)
        target_compile_definitions(${TOOL} PUBLIC binary_function=__binary_function)
    endif ()
endforeach()
