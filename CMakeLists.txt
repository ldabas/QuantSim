cmake_minimum_required(VERSION 3.15)
project(quantsim)

set(CMAKE_CXX_STANDARD 14)

add_executable(quantsim main.cpp gates.cpp gates.h types_temp.h)