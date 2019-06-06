cmake_minimum_required (VERSION 2.9)

project (KT-3P)

set(CMAKE_CXX_STANDARD 11)
include_directories("include")
include_directories("src")
include_directories("scripts")

find_package(ROOT REQUIRED MathMore)
if (ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIRS})
    link_directories(${ROOT_LIBRARY_DIRS})
else()
    message(SEND_ERROR "Root not found!")
endif()

add_executable( omega_BW
scripts/omega_BW.cpp
src/amplitude.cpp
src/pipi.cpp
src/omnes.cpp
src/aux_math.cpp
src/breit_wigner.cpp
src/dalitz.cpp
src/dalitz_fit.cpp
src/poly_exp.cpp
)
target_link_libraries( omega_BW ${ROOT_LIBRARIES})

add_executable( omega_KT
scripts/omega_KT.cpp
src/amplitude.cpp
src/pipi.cpp
src/omnes.cpp
src/aux_math.cpp
src/breit_wigner.cpp
src/dalitz.cpp
src/dalitz_fit.cpp
src/poly_exp.cpp
src/kt_isobar.cpp
src/kt_iteration.cpp
src/kt.cpp
)
target_link_libraries( omega_KT ${ROOT_LIBRARIES})
