cmake_minimum_required (VERSION 2.9)

project (KT-3P)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

include_directories("include")
include_directories("src")
include_directories("executables")

find_package(ROOT REQUIRED MathMore)
if (ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIRS})
    link_directories(${ROOT_LIBRARY_DIRS})
else()
    message(SEND_ERROR "Root not found!")
endif()

find_package(GSL REQUIRED)
if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    link_directories(${GSL_LIBRARY_DIRS})
else()
    message(SEND_ERROR "GSL not found!")
endif()

add_executable(BW
executables/BW.cpp
src/decay_kinematics.cpp
src/pipi.cpp
src/omnes.cpp
src/utilities.cpp
src/breit_wigner.cpp
src/dalitz.cpp
src/dalitz_fit.cpp
src/poly_exp.cpp
)
target_link_libraries(BW ${ROOT_LIBRARIES})

add_executable(once_subtracted
executables/once_subtracted.cpp
src/decay_kinematics.cpp
src/pipi.cpp
src/omnes.cpp
src/utilities.cpp
src/kt_isobar.cpp
src/kt_equations.cpp
src/kt_amplitude.cpp
src/kt_ang_integral.cpp
src/kt_disp_integral.cpp
src/dalitz.cpp
src/dalitz_fit.cpp
src/poly_exp.cpp
src/breit_wigner.cpp
)
target_link_libraries(once_subtracted ${ROOT_LIBRARIES})
target_link_libraries(once_subtracted ${GSL_LIBRARIES})

add_executable(twice_subtracted
executables/twice_subtracted.cpp
src/decay_kinematics.cpp
src/pipi.cpp
src/omnes.cpp
src/utilities.cpp
src/kt_isobar.cpp
src/kt_equations.cpp
src/kt_amplitude.cpp
src/kt_ang_integral.cpp
src/kt_disp_integral.cpp
src/dalitz.cpp
src/dalitz_fit.cpp
src/poly_exp.cpp
src/breit_wigner.cpp
)
target_link_libraries(twice_subtracted ${ROOT_LIBRARIES})
target_link_libraries(twice_subtracted ${GSL_LIBRARIES})
