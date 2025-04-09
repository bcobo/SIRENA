cmake_minimum_required(VERSION 3.16)
project(sirena VERSION "10.0.0")

include(CTest)
include(GNUInstallDirs)

message(STATUS "CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")


# Throw an error if CMAKE_INSTALL_PREFIX is not set by the user.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX not set!\n"
            "Use -DCMAKE_INSTALL_PREFIX=sixtedir to specify the SIXTE installation directory.")
endif()


### Global compiler flags ###
# TODO: Add -pedantic
set(GLOBAL_COMPILER_FLAGS
        -g -O2 -D_GNU_SOURCE -Wall -Wextra -Wshadow -Wformat=2 -Wunused)
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE c_std_99 cxx_std_17)
target_compile_options(compiler_flags INTERFACE
        $<BUILD_INTERFACE:${GLOBAL_COMPILER_FLAGS}>)


### SIMPUT ###
# By default, the SIMPUT installation is expected to be located at
# CMAKE_INSTALL_PREFIX. If SIMPUT is installed to a different location, an
# alternative path can be set via the SIMPUT_ROOT variable.
# Usage: cmake -DSIMPUT_ROOT=simputdir).
set(SIMPUT_ROOT ${CMAKE_INSTALL_PREFIX} CACHE STRING "Directory that contains the SIMPUT installation")
set(SIMPUT_INCLUDE_DIR ${SIMPUT_ROOT}/${CMAKE_INSTALL_INCLUDEDIR})
set(SIMPUT_LIBRARY_DIR ${SIMPUT_ROOT}/${CMAKE_INSTALL_LIBDIR})
set(SIMPUT_LIBRARIES ape_simput atFunctions_simput cfitsio_simput fftw3_simput heainit_simput heaio_simput heasp_simput heautils_simput labnh posstring simput wcs_simput-7.7)
find_library(SIMPUT_LIBRARY
        NAMES ${SIMPUT_LIBRARIES}
        PATHS ${SIMPUT_LIBRARY_DIR})
if(NOT SIMPUT_LIBRARY)
    message(FATAL_ERROR "Could not find the SIMPUT libraries at \n${SIMPUT_LIBRARY_DIR}\n"
            "Use -DSIMPUT_ROOT=simputdir to specify SIMPUT installation directory.")
endif()
message(STATUS "Using SIMPUT at ${SIMPUT_ROOT}")

### SIXTE ###
# By default, the SIXTE installation is expected to be located at
# CMAKE_INSTALL_PREFIX. If SIXTE is installed to a different location, an
# alternative path can be set via the SIXTE_ROOT variable.
# Usage: cmake -DSIXTE_ROOT=sixtedir).
set(SIXTE_ROOT ${CMAKE_INSTALL_PREFIX} CACHE STRING "Directory that contains the SIXTE installation")
set(SIXTE_INCLUDE_DIR ${SIXTE_ROOT}/${CMAKE_INSTALL_INCLUDEDIR})
set(SIXTE_LIBRARY_DIR ${SIXTE_ROOT}/${CMAKE_INSTALL_LIBDIR})

set(SIXTE_LIBRARIES libsixt) # TODO may need to add libs
find_library(SIXTE_LIBRARY
        NAMES ${SIXTE_LIBRARIES}
        PATHS ${SIXTE_LIBRARY_DIR})
if(NOT SIXTE_LIBRARY)
    message(FATAL_ERROR "Could not find the SIXTE libraries at \n${SIXTE_LIBRARY_DIR}\n"
            "Use -DSIXTE_ROOT=simputdir to specify SIXTE installation directory.")
endif()
message(STATUS "Using SIXTE at ${SIXTE_ROOT}")


### External packages ###
set(EXT_INCLUDE_DIRS)
set(EXT_LIBRARIES)

# GSL
find_package(GSL REQUIRED)
list(APPEND EXT_INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
list(APPEND EXT_LIBRARIES ${GSL_LIBRARIES})

# Boost
find_package(Boost REQUIRED)
list(APPEND EXT_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
list(APPEND EXT_LIBRARIES ${Boost_LIBRARIES})

# Threads
find_package(Threads REQUIRED)
list(APPEND EXT_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

# CURSES
set(CURSES_NEED_NCURSES TRUE)
find_package(Curses REQUIRED)
list(APPEND EXT_INCLUDE_DIRS ${CURSES_INCLUDE_DIRS})
list(APPEND EXT_LIBRARIES ${CURSES_LIBRARIES})
message(STATUS "Using CURSES_INCLUDE_DIRS at ${CURSES_INCLUDE_DIRS}")
message(STATUS "Using CURSES_LIBRARIES at ${CURSES_LIBRARIES}")

# SIXTE extlibs
list(APPEND EXT_INCLUDE_DIRS
          ${SIXTE_INCLUDE_DIR}/progressbar)
list(APPEND EXT_LIBRARIES progressbar)

# Additional libraries
list(APPEND EXT_LIBRARIES expat m ncurses readline)


### Configured header file ###
#configure_file(sixteconfig.h.in sixteconfig.h)
#install(FILES "${PROJECT_BINARY_DIR}/sixteconfig.h" DESTINATION include)


### Install sirena-install.sh/csh scripts ###
configure_file(sirena-install.sh.in sirena-install.sh @ONLY)
configure_file(sirena-install.csh.in sirena-install.csh @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/sirena-install.sh"
              "${PROJECT_BINARY_DIR}/sirena-install.csh"
        DESTINATION ${CMAKE_INSTALL_BINDIR})


### Project subdirectories ###
add_subdirectory(libsirena)
add_subdirectory(tools)

### Packaging ###
set(CPACK_SOURCE_IGNORE_FILES
        /\\\\.gitignore # Regex to match all .gitignore files. The leading "/"
                        # indicates start of regex and the double backslashes
                        # are used to escape the dot.
        .idea .git .gitlab-ci.yml .gitmodules)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}")
set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)

#if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
#        install(
#                SCRIPT ${CMAKE_SOURCE_DIR}/mac_solve_sixte_dependencies.sh
#        )
#endif()

message(STATUS "CMAKE_INSTALL_LIBDIR is ${CMAKE_INSTALL_LIBDIR}") 

