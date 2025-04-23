cmake_minimum_required(VERSION 3.5.0)
project(fast_af VERSION 0.1.0 LANGUAGES C CXX)

add_executable(fast_af fast_af.cpp)

# Add EAGLE
include(FetchContent)

FetchContent_Declare(
  EAGLE
  GIT_REPOSITORY https://github.com/Hexive-UN-03/EAGLE.git
  GIT_TAG        main
)

FetchContent_MakeAvailable(EAGLE)

# Manually point to HTSlib headers and link the system library
target_include_directories(fast_af PRIVATE /usr/include/htslib)
target_link_libraries(fast_af PRIVATE 
    EAGLE::EAGLE
    -lhts
    -ldeflate
    -lz
    -llzma
    -lbz2
)
