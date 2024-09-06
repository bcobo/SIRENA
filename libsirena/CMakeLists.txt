set(LIBSIRENA_SOURCES
        initSIRENA.c
        testriggerfile.c
        genutils.cpp
        inoututils.cpp
        integraSIRENA.cpp
        log.cpp
        pulseprocess.cpp
        scheduler.cpp
        tasksSIRENA.cpp)

set(LIBSIRENA_HEADERS
       genutils.h
       initSIRENA.h
       inoututils.h
       integraSIRENA.h
       log.h
       pulseprocess.h
       scheduler.h
       tasksSIRENA.h
       threadsafe_queue.h
       versionSIRENA.h)

#SET(LIB_NAME libsirena)
SET(LIB_NAME sirena)

add_library(${LIB_NAME} SHARED ${LIBSIRENA_SOURCES} ${LIBSIRENA_HEADERS})
set_target_properties(${LIB_NAME} PROPERTIES VERSION "9.1.5")
target_include_directories(${LIB_NAME} PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${EXT_INCLUDE_DIRS}
        ${SIMPUT_INCLUDE_DIR}
        ${SIXTE_INCLUDE_DIR})
target_link_directories(${LIB_NAME} PUBLIC
        ${SIMPUT_LIBRARY_DIR}
        ${SIXTE_LIBRARY_DIR})

target_link_libraries(${LIB_NAME} PUBLIC
        compiler_flags
        ${EXT_LIBRARIES}
        ${SIMPUT_LIBRARIES}
        ${SIXTE_LIBRARIES})

# INSTALL
set(installable_libs ${LIB_NAME})
install(TARGETS ${installable_libs})
install(FILES ${LIBSIRENA_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

### Apple Adjustments for missing boost headers ###
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        target_compile_definitions(${LIB_NAME} PUBLIC unary_function=__unary_function)
        target_compile_definitions(${LIB_NAME} PUBLIC binary_function=__binary_function)
        target_compile_definitions(${LIB_NAME} PUBLIC sincos=__sincos)
endif ()
