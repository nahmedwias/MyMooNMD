# CMakeLists file for the creation of a static and dynamic library from Triangle.

cmake_minimum_required(VERSION 3.0)

project(Triangle)

add_library(triangle_LINUX64 STATIC triangle.c)

add_library(triangle_LINUX64 SHARED triangle.c)