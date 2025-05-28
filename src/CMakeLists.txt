# CMakeLists.txt in src/

add_library(utils
    utils.cpp
)

# Dichiara il path degli header locali (per includere "Utils.hpp")
target_include_directories(utils
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${EIGEN3_INCLUDE_DIR}    # 👈 aggiunto per garantire visibilità a Eigen
)
