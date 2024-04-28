cmake_minimum_required(VERSION 3.28)
project(GeneticAlgorithm)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(target "GeneticAlgorithm")

add_executable(${target} main.cpp
        GeneticAlgorithm.h
        GeneticAlgorithm.cpp
        DataParserQuick.h
        Individuals.h
        Population.h
        Population.cpp
)

set(-DCMAKE_BUILD_TYPE=Release)