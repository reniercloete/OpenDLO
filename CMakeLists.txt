cmake_minimum_required (VERSION 3.0.0 FATAL_ERROR)

# Project name
project (OpenDLO)

add_compile_definitions( TRILIBRARY NO_TIMER ANSI_DECLARATORS PARDISO_BARRIER)

set (DLO_COIN_OR_DIR_DESC "Path to the COIN-OR repo")
if (NOT DEFINED DLO_COIN_OR_DIR)
  set (DLO_COIN_OR_DIR "" CACHE PATH ${DLO_COIN_OR_DIR_DESC})
  get_filename_component (DLO_COIN_OR_DIR "${DLO_COIN_OR_DIR}" ABSOLUTE)
else()
  file (TO_CMAKE_PATH "${DLO_COIN_OR_DIR}" DLO_COIN_OR_DIR)
  set (DLO_COIN_OR_DIR "${DLO_COIN_OR_DIR}" CACHE PATH "${DLO_COIN_OR_DIR_DESC}" FORCE)
endif()

set (DLO_MKL_DIR_DESC "Path to Intel MKL")
if (NOT DEFINED DLO_MKL_DIR)
  set (DLO_MKL_DIR "" CACHE PATH ${DLO_MKL_DIR_DESC})
  get_filename_component (DLO_MKL_DIR "${DLO_MKL_DIR}" ABSOLUTE)
else()
  file (TO_CMAKE_PATH "${DLO_MKL_DIR}" DLO_MKL_DIR)
  set (DLO_MKL_DIR "${DLO_MKL_DIR}" CACHE PATH "${DLO_MKL_DIR_DESC}" FORCE)
endif()

set (DLO_MOSEK_DIR_DESC "Path to Mosek")
if (NOT DEFINED DLO_MOSEK_DIR)
  set (DLO_MOSEK_DIR "" CACHE PATH ${DLO_MOSEK_DIR_DESC})
  get_filename_component (DLO_MOSEK_DIR "${DLO_MOSEK_DIR}" ABSOLUTE)
else()
  file (TO_CMAKE_PATH "${DLO_MOSEK_DIR}" DLO_MOSEK_DIR)
  set (DLO_MOSEK_DIR "${DLO_MOSEK_DIR}" CACHE PATH "${DLO_MOSEK_DIR_DESC}" FORCE)
endif()

set (DLO_GLFW_DIR_DESC "Path to GLFW")
if (NOT DEFINED DLO_GLFW_DIR)
  set (DLO_GLFW_DIR "" CACHE PATH ${DLO_GLFW_DIR_DESC})
  get_filename_component (DLO_GLFW_DIR "${DLO_GLFW_DIR}" ABSOLUTE)
else()
  file (TO_CMAKE_PATH "${DLO_GLFW_DIR}" DLO_GLFW_DIR)
  set (DLO_GLFW_DIR "${DLO_GLFW_DIR}" CACHE PATH "${DLO_GLFW_DIR_DESC}" FORCE)
endif()

# Add executable
ADD_LIBRARY (OpenDLOLib STATIC ${SOURCE_FILES}
	src/CoinDLOSolver.cpp
	src/CoinDLOSolver.h
	src/Constants.h
	src/DLOSolver.cpp
	src/DLOSolver.h
	src/Domain.cpp
	src/Domain.h
	src/Edge.cpp
	src/Edge.h
	src/Enums.h
	src/Line2D.cpp
	src/Line2D.h
	src/MosekDLOSolver.cpp
	src/MosekDLOSolver.h
	src/Node.h
	src/Point2D.cpp
	src/Point2D.h
	src/Poly2D.cpp
	src/Poly2D.h
	src/triangle.cpp
	src/triangle.h
	src/Vector2d.h
	)

target_include_directories(OpenDLOLib PRIVATE ${DLO_COIN_OR_DIR}/Clp/src ${DLO_COIN_OR_DIR}/CoinUtils/src)
target_include_directories(OpenDLOLib PRIVATE ${DLO_MKL_DIR}/mkl/latest/include)
target_include_directories(OpenDLOLib PRIVATE ${DLO_MOSEK_DIR})

# Add executable
add_executable (OpenDLOTest
  src/OpenDLOTest.cpp
  ${DLO_GLFW_DIR}/deps/glad_gl.c
)

target_include_directories(OpenDLOTest PRIVATE ${DLO_COIN_OR_DIR}/Clp/src ${DLO_COIN_OR_DIR}/CoinUtils/src)
target_include_directories(OpenDLOTest PRIVATE ${DLO_MKL_DIR}/mkl/latest/include)
target_include_directories(OpenDLOTest PRIVATE ${DLO_MOSEK_DIR})
target_include_directories(OpenDLOTest PRIVATE ${DLO_GLFW_DIR}/include)
target_include_directories(OpenDLOTest PRIVATE ${DLO_GLFW_DIR}/deps)

target_link_libraries(OpenDLOTest debug ${DLO_GLFW_DIR}/build/src/Debug/glfw3.lib )
target_link_libraries(OpenDLOTest debug OpenDLOLib)
target_link_libraries(OpenDLOTest debug ${DLO_COIN_OR_DIR}/Clp/MSVisualStudio/v10/x64-v143-Debug/libClp.lib)
target_link_libraries(OpenDLOTest debug ${DLO_COIN_OR_DIR}/Clp/MSVisualStudio/v10/x64-v143-Debug/libCoinUtils.lib)
target_link_libraries(OpenDLOTest debug ${DLO_MKL_DIR}/mkl/latest/lib/intel64/mkl_intel_lp64.lib)
target_link_libraries(OpenDLOTest debug ${DLO_MKL_DIR}/mkl/latest/lib/intel64/mkl_intel_thread.lib)
target_link_libraries(OpenDLOTest debug ${DLO_MKL_DIR}/mkl/latest/lib/intel64/mkl_core.lib)
target_link_libraries(OpenDLOTest debug ${DLO_MKL_DIR}/compiler/latest/windows/compiler/lib/intel64_win/libiomp5md.lib)

target_link_libraries(OpenDLOTest optimized ${DLO_GLFW_DIR}/build/src/Release/glfw3.lib)
target_link_libraries(OpenDLOTest optimized OpenDLOLib)
target_link_libraries(OpenDLOTest optimized ${DLO_COIN_OR_DIR}/Clp/MSVisualStudio/v10/x64-v143-Release/libClp.lib)
target_link_libraries(OpenDLOTest optimized ${DLO_COIN_OR_DIR}/Clp/MSVisualStudio/v10/x64-v143-Release/libCoinUtils.lib)
target_link_libraries(OpenDLOTest debug ${DLO_MKL_DIR}/mkl/latest/lib/intel64/mkl_intel_lp64.lib)
target_link_libraries(OpenDLOTest debug ${DLO_MKL_DIR}/mkl/latest/lib/intel64/mkl_intel_thread.lib)
target_link_libraries(OpenDLOTest debug ${DLO_MKL_DIR}/mkl/latest/lib/intel64/mkl_core.lib)
target_link_libraries(OpenDLOTest debug ${DLO_MKL_DIR}/compiler/latest/windows/compiler/lib/intel64_win/libiomp5md.lib)



