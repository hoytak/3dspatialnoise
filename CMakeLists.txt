project(3dspatialnoise)
cmake_minimum_required(VERSION 2.4)

if(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
endif()

option(USE_OPENMP "Enable mutithreaded processing using OpenMP, if available" Yes)

################################################################################
# Locate needed libraries

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(CheckCXXCompilerFlag)

SET(CMAKE_CXX_FLAGS "-g -fPIC -Wall -Wcast-qual -Wno-sign-compare -Wno-write-strings -Wno-deprecated -Wno-deprecated-declarations -Wno-cast-qual -Wno-reorder -Wno-unused-label")

set(STANDARD_LIBRARIES_TO_LINK_TO "")
set(STANDARD_LINKER_FLAGS "")

########################################
# Get the local include paths
set(CMAKE_INCLUDE_PATH "${CMAKE_INCLUDE_PATH}:$ENV{CMAKE_INCLUDE_PATH}:$ENV{INCLUDE_PATH}")
string(REGEX REPLACE ":+" ";" CMAKE_INCLUDE_PATH "${CMAKE_INCLUDE_PATH}")


########################################
# Get the local library include paths
set(CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH}:$ENV{CMAKE_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
string(REGEX REPLACE ":+" ";" CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH}")

link_directories("${CMAKE_INSTALL_PREFIX}/lib")
include_directories("${CMAKE_INSTALL_PREFIX}/include")

find_package(StandardMathLibrary)

if(NOT STANDARD_MATH_LIBRARY_FOUND)

  message(FATAL_ERROR
    "Can't link to the standard math library.")
else()
  
  set(STANDARD_LIBRARIES_TO_LINK_TO "${STANDARD_LIBRARIES_TO_LINK_TO};${STANDARD_MATH_LIBRARY}")
endif()

############################################################
# Find BLAS

find_package(BLAS)

if(BLAS_FOUND)
  set(STANDARD_LIBRARIES_TO_LINK_TO "${STANDARD_LIBRARIES_TO_LINK_TO};${BLAS_LIBRARIES}")
  set(STANDARD_LINKER_FLAGS "${STANDARD_LINKER_FLAGS};${BLAS_LINKER_FLAGS}")
  message("BLAS LIB: ${BLAS_LIBRARIES}")
  message("BLAS LIB_DIR: ${BLAS_LIBRARIES_DIR}")
  link_directories("${BLAS_LIBRARIES_DIR}")
else()
  message(FATAL_ERROR "BLAS libraries not found!")
endif()

################################################################################
# Check for compiler options; particularly OpenMP

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fipa-struct-reorg -fipa-cp -ftree-loop-linear -mfpmath=sse -ffast-math -fstrict-aliasing -floop-interchange -floop-block -floop-strip-mine -fno-common -Wno-long-long -Wwrite-strings")

  set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -UNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -DNDEBUG -march=native -fipa-struct-reorg -fipa-cp -ftree-loop-linear -mfpmath=sse -ffast-math -fstrict-aliasing -floop-interchange -floop-block -floop-strip-mine")

  if(CMAKE_BUILD_TYPE_TOLOWER STREQUAL "debug")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  elseif(CMAKE_BUILD_TYPE_TOLOWER STREQUAL "release")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  endif()
    
  check_cxx_compiler_flag("-Wno-variadic-macros" COMPILER_SUPPORT_WNOVARIADICMACRO)
  if(COMPILER_SUPPORT_WNOVARIADICMACRO)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-variadic-macros")
  endif()

  if(USE_OPENMP)

    check_cxx_compiler_flag("-fopenmp" COMPILER_SUPPORT_OPENMP)

    if(COMPILER_SUPPORT_OPENMP)
      set(Boost_USE_MULTITHREADED  ON)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DUSE_OPENMP")
      set(STANDARD_LIBRARIES_TO_LINK_TO "${STANDARD_LIBRARIES_TO_LINK_TO};gomp")
      message("Enabling OpenMP, pass -DUSE_OPENMP=No to cmake to enable it.")
    else()
      set(Boost_USE_MULTITHREADED  OFF)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -UUSE_OPENMP")
      message("Disabling OpenMP; compiler doesn't seem to support it.")
    endif()
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -UUSE_OPENMP")
    message("Disabling OpenMP, pass -DUSE_OPENMP=Yes to cmake to enable it.")
  endif()
  
endif(CMAKE_COMPILER_IS_GNUCXX)

if(MSVC)
  # C4127 - conditional expression is constant
  # C4714 - marked as __forceinline not inlined (I failed to deactivate it selectively)
  #         We can disable this warning in the unit tests since it is clear that it occurs
  #         because we are oftentimes returning objects that have a destructor or may
  #         throw exceptions - in particular in the unit tests we are throwing extra many
  #         exceptions to cover indexing errors.
  # C4505 - unreferenced local function has been removed (impossible to deactive selectively)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /wd4127 /wd4505 /wd4714")

  # replace all /Wx by /W4
  string(REGEX REPLACE "/W[0-9]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

  check_cxx_compiler_flag("/openmp" COMPILER_SUPPORT_OPENMP)
  if(COMPILER_SUPPORT_OPENMP)
    option(USE_OPENMP "Enable/Disable OpenMP." OFF)
    if(USE_OPENMP)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
      set(Boost_USE_MULTITHREADED  ON)
      message("Enabling OpenMP in tests/examples")
    else()
      set(Boost_USE_MULTITHREADED  OFF)
    endif()
  endif()

  option(USE_SSE2 "Enable/Disable SSE2 in tests/examples" OFF)
  if(USE_SSE2)
    if(NOT CMAKE_CL_64)
      # arch is not supported on 64 bit systems, SSE is enabled automatically.
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
    endif(NOT CMAKE_CL_64)
    message("Enabling SSE2 in tests/examples")
  endif(USE_SSE2)
endif(MSVC)

########################################
# Boost

find_package( Boost )

if(Boost_FOUND)
  message("BOOST version ${Boost_VERSION} found.")
  include_directories(${Boost_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "BOOST library not found!")
endif()

########################################
# GSL

find_package( GSL )

if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIR})
  link_directories(${GSL_LINK_DIRECTORIES})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_GSL_CXX_FLAGS}")
  set(STANDARD_LIBRARIES_TO_LINK_TO "${STANDARD_LIBRARIES_TO_LINK_TO};gsl")
else()
  message(FATAL_ERROR "GSL Libraries not found!")
endif()

add_definitions(${CMAKE_CXX_FLAGS})

################################################################################
# Install the C++ libraries
add_subdirectory(src)

install(
  DIRECTORY "src/"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/include/3dspatialnoise"
  FILES_MATCHING PATTERN "*.h*" 
)