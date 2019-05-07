cmake_minimum_required (VERSION 2.6)

project (KT-3P)

set(CMAKE_CXX_STANDARD 11)
include_directories("include")
include_directories("src")
include_directories("scripts")

find_package(ROOT REQUIRED)
if (ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIRS})
    link_directories(${ROOT_LIBRARY_DIRS})
else()
    message(SEND_ERROR "Root not found!")
endif()

add_executable( omega_BW
scripts/omega_BW.cpp
src/amp.cpp
src/pipi.cpp
src/omnes.cpp
src/aux_math.cpp
src/dalitz.cpp
)
