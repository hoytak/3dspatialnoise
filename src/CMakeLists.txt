project(3dspatialnoise)

add_definitions(${CMAKE_CXX_FLAGS})

# For the python extension; note; this requires STATIC above
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
link_directories("${CMAKE_CURRENT_SOURCE_DIR}")

