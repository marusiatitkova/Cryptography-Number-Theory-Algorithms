cmake_minimum_required(VERSION 3.12)
project(Crypography_NumberTheory_Algorithms)

set(CMAKE_CXX_STANDARD 14)

add_executable(Crypography_NumberTheory_Algorithms main.cpp BigInteger.cpp BigInteger.h Algorithms.cpp Algorithms.h EllipticCurve.h ECPoint.h)