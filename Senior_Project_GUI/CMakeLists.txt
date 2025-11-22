cmake_minimum_required(VERSION 3.16)
project(DNAVisualizer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 COMPONENTS Core Widgets REQUIRED)

set(SOURCES
    main.cpp
    core/DNASequence.cpp
    core/Feature.cpp
    features/Gene.cpp
    features/SearchMatch.cpp
    algorithms/Aligner.cpp
    algorithms/OrfFinder.cpp
    algorithms/TrieIndex.cpp
    algorithms/TrieNode.cpp
    algorithms/SearchAlgorithm.cpp
    io/FastaParser.cpp
    gui/Controller.cpp
    gui/ResultsPanel.cpp
    gui/mainwindow.cpp
    gui/Sequence.cpp
)

set(HEADERS
    core/DNASequence.h
    core/Feature.h
    features/Gene.h
    features/SearchMatch.h
    algorithms/Aligner.h
    algorithms/OrfFinder.h
    algorithms/TrieIndex.h
    algorithms/TrieNode.h
    algorithms/SearchAlgorithm.h
    io/FastaParser.h
    gui/Controller.h
    gui/ResultsPanel.h
    gui/mainwindow.h
    gui/Sequence.h
)

add_executable(DNAVisualizer ${SOURCES} ${HEADERS})

target_link_libraries(DNAVisualizer Qt6::Core Qt6::Widgets)

target_include_directories(DNAVisualizer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})