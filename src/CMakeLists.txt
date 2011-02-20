add_library (spatialnoise SHARED
  physics.h
  cutting_planes.cpp
  physics.cpp
  ground.cpp
  ground.h
  ground_base.h
  ground_base.hpp
  typedefs.h
  boost_shared_ptr.h
  cutting_planes.h
  geometry.h
)

target_link_libraries(spatialnoise ${STANDARD_LIBRARIES_TO_LINK_TO})

install(TARGETS spatialnoise DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")