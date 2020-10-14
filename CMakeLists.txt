cmake_minimum_required(VERSION 2.8)
project(ObliviousYannakakisDemo)

set(CMAKE_CXX_STANDARD 11)
ENABLE_LANGUAGE(C ASM)

add_definitions("-g -O2 -Wall -msse2 -msse -march=native -maes")

add_library(blake STATIC blake3/c/blake3.c blake3/c/blake3_dispatch.c blake3/c/blake3_portable.c blake3/c/blake3_sse41_x86-64_unix.S blake3/c/blake3_avx2_x86-64_unix.S blake3/c/blake3_avx512_x86-64_unix.S)

find_path (ntlpath NTL/ZZ.h)
include_directories(${ntlpath})

find_library(NTL_LIB ntl)
if(NOT NTL_LIB)
  message(FATAL_ERROR "ntl library not found.")
endif()

find_library(GMP_LIB gmp)
if(NOT GMP_LIB)
  message(FATAL_ERROR "gmp library not found.")
endif()

add_executable(demo demo.cpp src/OT.cpp src/aes.cpp src/OPRF.cpp src/PSI.cpp src/poly.cpp src/MurmurHash3.cpp src/PRNG.cpp src/OEP.cpp src/utils.cpp src/relation.cpp src/GC.cpp)
target_link_libraries(demo ${NTL_LIB} ${GMP_LIB} blake)