######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists.txt
# @created     : Sunday Dec 15, 2024 15:53:26 CET
######################################################################

set(SSTATS_CORE_HEADERS 
  ${SSTATS_CORE_DIR}/MSConfig.h
  ${SSTATS_CORE_DIR}/MSDataPoint.h
  ${SSTATS_CORE_DIR}/MSDataPointVector.h
  ${SSTATS_CORE_DIR}/MSMath.h
  ${SSTATS_CORE_DIR}/MSMinimizer.h
  ${SSTATS_CORE_DIR}/MSModel.h
  ${SSTATS_CORE_DIR}/MSModelPulls.h
  ${SSTATS_CORE_DIR}/MSModelTHnBMLF.h
  ${SSTATS_CORE_DIR}/MSObject.h
  ${SSTATS_CORE_DIR}/MSParameter.h
  ${SSTATS_CORE_DIR}/MSPDFBuilderTHn.h
  ${SSTATS_CORE_DIR}/MSTHnHandler.h
  )

set(SSTATS_CORE_SOURCES
  ${SSTATS_CORE_DIR}/MSConfig.cxx
  ${SSTATS_CORE_DIR}/MSDataPoint.cxx
  ${SSTATS_CORE_DIR}/MSMath.cxx
  ${SSTATS_CORE_DIR}/MSMinimizer.cxx
  ${SSTATS_CORE_DIR}/MSModel.cxx
  ${SSTATS_CORE_DIR}/MSModelPulls.cxx
  ${SSTATS_CORE_DIR}/MSModelTHnBMLF.cxx
  ${SSTATS_CORE_DIR}/MSParameter.cxx
  ${SSTATS_CORE_DIR}/MSPDFBuilderTHn.cxx
  ${SSTATS_CORE_DIR}/MSTHnHandler.cxx
  )

# SLArReadoutSystemConfig dictionary 
add_library(sstats-core SHARED
  ${SSTATS_CORE_SOURCES} ${SSTATS_CORE_HEADERS}
)
target_link_libraries(sstats-core PRIVATE ${ROOT_LIBRARIES})
target_include_directories(sstats-core
  PUBLIC
  $<BUILD_INTERFACE:${SSTATS_CORE_DIR}>
  $<INSTALL_INTERFACE:"${SSTATS_INC_DEST}/core">
)
ROOT_GENERATE_DICTIONARY(G__sstats-core 
  ${SSTATS_CORE_DIR}/*.h
  MODULE sstats-core
  LINKDEF ${SSTATS_CORE_DIR}/m-stats-core_LinkDef.h)

#----------------------------------------------------------------------
# Include files in header and sources list
set(solar-stats-sources ${solar-stats-sources} ${SSTATS_CORE_SOURCES} PARENT_SCOPE)
set(solar-stats-headers ${solar-stats-headers} ${SSTATS_CORE_HEADERS} PARENT_SCOPE)

#----------------------------------------------------------------------
# Install instructions

set(sstats_core_libs
  sstats-core
  )

set(sstats_core_resources 
  ${CMAKE_CURRENT_BINARY_DIR}/G__sstats-core.cxx
  ${CMAKE_CURRENT_BINARY_DIR}/libsstats-core.rootmap 
  ${CMAKE_CURRENT_BINARY_DIR}/libsstats-core_rdict.pcm
  )

foreach(lib IN LISTS sstats_core_libs)
  set_target_properties(${lib} PROPERTIES 
    BUILD_RPATH ${CMAKE_CURRENT_BINARY_DIR}
    INSTALL_RPATH ${SSTATS_LIB_DIR}
    )

  install(TARGETS ${lib}
    EXPORT                     "${PROJECT_NAME}Targets" 
    LIBRARY DESTINATION        "${SSTATS_LIB_DIR}" 
    RUNTIME DESTINATION        "${SSTATS_LIB_DIR}")
endforeach()

install(FILES ${SSTATS_CORE_HEADERS} DESTINATION "${SSTATS_INC_DEST}/core")
install(FILES ${sstats_core_resources} DESTINATION "${SSTATS_LIB_DIR}")
