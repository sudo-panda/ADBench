project("ADBench_Clad" CXX)

add_library("Clad" MODULE)
target_sources("Clad" PRIVATE "${CMAKE_SOURCE_DIR}/src/cpp/shared/utils.cpp")

target_sources("Clad" PRIVATE "ba/ba.c")
target_sources("Clad" PRIVATE "ba/ba_clad.cpp")
target_sources("Clad" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/CladBA.cpp")
