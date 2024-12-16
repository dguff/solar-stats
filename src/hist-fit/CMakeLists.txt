######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists
# @created     : Sunday Dec 15, 2024 16:33:49 CET
######################################################################

# Add the shared library
add_library(sstats-hist-fit SHARED ${SSTATS_HISTFIT_DIR}/MSHistFit.cxx)
target_link_libraries(sstats-hist-fit PRIVATE ${ROOT_LIBRARIES} sstats-core)
target_include_directories(sstats-hist-fit
  PUBLIC
  $<BUILD_INTERFACE:${SSTATS_HISTFIT_DIR}>
  $<INSTALL_INTERFACE:"${SSTATS_INC_DEST}/hist-fit">
)
# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__sstats-hist-fit
  LINKDEF ${SSTATS_HISTFIT_DIR}/m-stats-hist-fit_LinkDef.h
  MODULE sstats-hist-fit
)

# Executable source
set(EXECUTABLE_NAME solar_stats)
set(EXECUTABLE_SOURCE ${SSTATS_HISTFIT_DIR}/m-stats-hist-fit.cxx)

# Add the executable
add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SOURCE})
target_link_libraries(${EXECUTABLE_NAME} PUBLIC sstats-hist-fit sstats-core ${ROOT_LIBRARIES} ROOT::Minuit)
set_target_properties(${EXECUTABLE_NAME} PROPERTIES
  BUILD_WITH_INSTALL_RPATH TRUE
  INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
)

# Installation for library headers and executables
install(TARGETS sstats-hist-fit LIBRARY DESTINATION ${SSTATS_LIB_DIR})
install(TARGETS ${EXECUTABLE_NAME} RUNTIME DESTINATION ${SSTATS_BIN_DIR})


