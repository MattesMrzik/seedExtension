project(seedExtensionProgram)

# sources as library to make them testable
add_library(seedExtensionLib STATIC AllTips.cpp AllTips.hpp
                                    VisualizeGraph.hpp VisualizeGraph.cpp
									Configuration.h
									ExtendSeed.cpp ExtendSeed.hpp
									PathBundleTip.cpp PathBundleTip.hpp)
target_include_directories(seedExtensionLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# set C++ standard
target_compile_features(seedExtensionLib PUBLIC cxx_std_17)

# link metagraph (important that this comes first)
target_link_libraries(seedExtensionLib PUBLIC metagraphInterface)

# link third party libraries
target_include_directories(seedExtensionLib SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(seedExtensionLib PUBLIC cxx-prettyprint)
target_link_libraries(seedExtensionLib PRIVATE Boost::program_options)

# the main program
add_executable(seedExtension main.cpp)

# link sources
target_link_libraries(seedExtension PRIVATE seedExtensionLib)
