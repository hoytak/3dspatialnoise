if(BUILD_PYTHON_EXTENSION)
  # get the library paths correct.
  get_directory_property(CMAKE_INCLUDE_DIRECTORIES DIRECTORY "../src" "INCLUDE_DIRECTORIES")
  get_directory_property(CMAKE_LINK_DIRECTORIES DIRECTORY "../src" "LINK_DIRECTORIES")

  function(get_argument_list DEST)
    string(REGEX REPLACE "${DEST}[; ]*" "" V "${ARGV}")
    string(REGEX REPLACE "[; :]+" ":" V "${V}")
    set(${DEST} "${V}" PARENT_SCOPE)
  endfunction()

  get_argument_list(PYTHON_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  message("PYTHON_CXX_FLAGS = ${PYTHON_CXX_FLAGS}")
  get_argument_list(PYTHON_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS})
  message("PYTHON_LINKER_FLAGS = ${CMAKE_SHARED_LINKER_FLAGS}")
  get_argument_list(PYTHON_INCLUDE_DIRECTORIES ${CMAKE_INCLUDE_DIRECTORIES})
  message("PYTHON_INCLUDE_DIRECTORIES = ${PYTHON_INCLUDE_DIRECTORIES}")
  get_argument_list(PYTHON_LIBRARY_DIRECTORIES ${CMAKE_LINK_DIRECTORIES})
  message("PYTHON_LIBRARY_DIRECTORIES = ${PYTHON_LIBRARY_DIRECTORIES}")
  get_argument_list(PYTHON_LIBRARIES ${STANDARD_LIBRARIES_TO_LINK_TO})
  message("PYTHON_LIBRARIES=${PYTHON_LIBRARIES}")

  SET(PYTHON_CALL_ARGUMENTS
    "CMAKE_LINKER_FLAGS=${PYTHON_LINKER_FLAGS}" 
    "CMAKE_CXX_FLAGS=${PYTHON_CXX_FLAGS}" 
    "CMAKE_INCLUDE_DIRECTORIES=${PYTHON_INCLUDE_DIRECTORIES}" 
    "CMAKE_LINK_DIRECTORIES=../src:${PYTHON_LIBRARY_DIRECTORIES}"
    "CMAKE_LIBRARIES=${PYTHON_LIBRARIES}")

  SET(PYTHON_TEMP_BUILD_DIR "${CMAKE_BINARY_DIR}/CMakeFiles/python-build")

  SET(PYTHON_INSTALL_COMMAND 
    python setupcmake.py install --prefix=${CMAKE_INSTALL_PREFIX} ${PYTHON_CALL_ARGUMENTS})

  SET(PYTHON_BUILD_COMMAND 
    python setupcmake.py build ${PYTHON_CALL_ARGUMENTS})

  add_custom_target(spatial_noise_python_build ALL ${PYTHON_BUILD_COMMAND}
    COMMENT "Building python extension.")

  install(CODE "execute_process(COMMAND ${PYTHON_INSTALL_COMMAND} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})")

  # install(TARGETS qbpmm_python_install)

endif()
