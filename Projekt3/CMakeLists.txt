cmake_minimum_required(VERSION 3.14)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/external/vcpkg/scripts/buildsystems/vcpkg.cmake"
        CACHE STRING "Vcpkg toolchain file")
endif()

set(VCPKG_INSTALLED_DIR "${CMAKE_SOURCE_DIR}/vcpkg_installed")

set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-windows/share")

set(PYTHON_EXECUTABLE "${CMAKE_SOURCE_DIR}/envs/sygnalytest_env/python.exe" CACHE FILEPATH "Python interpreter")

project(signals)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



find_package(pybind11 REQUIRED)
add_subdirectory(external/matplotplusplus)
#add_subdirectory(external/pybind11)

find_package(fftw3 REQUIRED)
find_package(libpng REQUIRED)
find_package(jpeg REQUIRED)
find_package(lapack REQUIRED)
find_package(liblzma REQUIRED)
find_package(nodesoup REQUIRED)
find_package(openblas REQUIRED)
find_package(tiff REQUIRED)
find_package(libjpeg-turbo REQUIRED)
find_package(zlib REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/src)

pybind11_add_module(signals
    src/binding.cpp
    src/main.cpp
)


set_target_properties(signals PROPERTIES
    PREFIX ""  # Removes the "lib" prefix on Unix
    SUFFIX ".pyd"  # Makes it .pyd instead of .so or .dll
)
#add_executable(sygnalytest src/main.cpp)

# Include directories
#target_include_directories(sygnalytest PRIVATE
    #"${CMAKE_SOURCE_DIR}/envs/sygnalytest_env/include"
    #"${CMAKE_SOURCE_DIR}/envs/sygnalytest_env/Library/include"
    #"${VCPKG_INSTALLED_DIR}/x64-windows/include"
#)

# Link directories
#link_directories(
    #"${CMAKE_SOURCE_DIR}/envs/sygnalytest_env/Library/lib"
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib"
#)
target_link_libraries(signals PRIVATE matplot)

# Link libraries
#target_link_libraries(sygnalytest PRIVATE
    #matplot
    #pybind11::module
    #fftw3
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/fftw3.lib"
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/fftw3f.lib"
    #fftw3f
    #fftw3l
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/fftw3l.lib"
    #jpeg
   # "${VCPKG_INSTALLED_DIR}/x64-windows/lib/jpeg.lib"
    #lapack
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/lapack.lib"
    #libpng16
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/libpng16.lib"
    #PNG::PNG
    #liblzma
   # "${VCPKG_INSTALLED_DIR}/x64-windows/lib/lzma.lib"
   # nodesoup
    #"${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-windows/lib/nodesoup.lib"
    #openblas
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/openblas.lib"
    #tiff
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/tiff.lib"
   #turbojpeg
   #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/turbojpeg.lib"
    #zlib
    #"${VCPKG_INSTALLED_DIR}/x64-windows/lib/zlib.lib"




