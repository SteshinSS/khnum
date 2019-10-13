cmake_minimum_required(VERSION 3.0)

project(catch)

add_library(catch INTERFACE)
target_include_directories(catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})