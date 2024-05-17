cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)

file(GLOB_RECURSE files src/*.cc)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")


add_executable (EBSDConograph ${files})
