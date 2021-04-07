cmake_minimum_required(VERSION 3.17)
project(TransitionMatrixRoute)

set(CMAKE_CXX_STANDARD 17)
include_directories(include)
find_package(Threads)
find_package(Boost COMPONENTS system filesystem REQUIRED)
file(GLOB source_files
        "include/*.h"
        "include/*.hpp"
        )

add_executable(TransitionMatrixRoute main.cpp)
target_link_libraries( TransitionMatrixRoute
        pthread
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_FILESYSTEM_LIBRARY}
        SQLiteCpp sqlite3
        )