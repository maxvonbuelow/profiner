find_path(NVBIT_INCLUDE_DIR nvbit.h
          HINTS
          "${CMAKE_CURRENT_SOURCE_DIR}/../nvbit"
          PATH_SUFFIXES libraw
         )

find_library(NVBIT_LIBRARIES NAMES nvbit
             HINTS
             "${CMAKE_CURRENT_SOURCE_DIR}/../nvbit"
            )

if(NVBIT_INCLUDE_DIR AND NVBIT_LIBRARIES)
	add_library(NVBit INTERFACE IMPORTED)
	set_property(TARGET NVBit PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${NVBIT_INCLUDE_DIR}")
	set_property(TARGET NVBit PROPERTY INTERFACE_LINK_LIBRARIES "${NVBIT_LIBRARIES}")
endif()
