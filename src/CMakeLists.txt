list(APPEND SolidiPlatonici_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)
list(APPEND SolidiPlatonici_sources ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp)

list(APPEND SolidiPlatonici_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(SolidiPlatonici_sources ${SolidiPlatonici_sources} PARENT_SCOPE)
set(SolidiPlatonici_headers ${SolidiPlatonici_headers} PARENT_SCOPE)
set(SolidiPlatonici_includes ${SolidiPlatonici_includes} PARENT_SCOPE)