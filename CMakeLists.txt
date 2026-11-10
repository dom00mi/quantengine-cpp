cmake_minimum_required(VERSION 3.21)
project(QuantEngine VERSION 0.2.1 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

option(QUANTENGINE_BUILD_EXAMPLES "Build sample pricing applications" ON)

add_library(quantengine STATIC
  src/analytics.cpp
  src/calibration.cpp
  src/curves.cpp
  src/market.cpp
  src/model_calibration.cpp
  src/models.cpp
  src/multifactor.cpp
  src/numerics.cpp
  src/optimisation.cpp
  src/pde.cpp
  src/rates.cpp
  src/risk.cpp
  src/volatility_surface.cpp)
add_library(QuantEngine::QuantEngine ALIAS quantengine)
target_compile_features(quantengine PUBLIC cxx_std_20)
target_include_directories(quantengine PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(quantengine PROPERTIES EXPORT_NAME QuantEngine)

if(MSVC)
  target_compile_options(quantengine PRIVATE /W4)
else()
  target_compile_options(quantengine PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(QUANTENGINE_BUILD_EXAMPLES)
  add_executable(pricing_workbench examples/pricing_workbench.cpp)
  target_link_libraries(pricing_workbench PRIVATE QuantEngine::QuantEngine)
  add_executable(expanded_workbench examples/expanded_workbench.cpp)
  target_link_libraries(expanded_workbench PRIVATE QuantEngine::QuantEngine)
endif()

if(BUILD_TESTING)
  add_executable(quantengine_tests tests/test_quantengine.cpp)
  target_link_libraries(quantengine_tests PRIVATE QuantEngine::QuantEngine)
  add_test(NAME quantengine_regression COMMAND quantengine_tests)
  add_executable(quantengine_extended_tests tests/test_extended.cpp)
  target_link_libraries(quantengine_extended_tests PRIVATE QuantEngine::QuantEngine)
  add_test(NAME quantengine_extended COMMAND quantengine_extended_tests)
endif()

install(TARGETS quantengine EXPORT QuantEngineTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
configure_package_config_file(cmake/QuantEngineConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/QuantEngineConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QuantEngine)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/QuantEngineConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/QuantEngineConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/QuantEngineConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QuantEngine)
install(EXPORT QuantEngineTargets FILE QuantEngineTargets.cmake
  NAMESPACE QuantEngine:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QuantEngine)
