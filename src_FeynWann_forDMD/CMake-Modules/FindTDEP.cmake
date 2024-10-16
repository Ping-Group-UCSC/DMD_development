find_library(TDEP_OLLE_LIBRARY NAMES olle PATHS ${TDEP_PATH} ${TDEP_PATH}/lib ${TDEP_PATH}/lib64 NO_DEFAULT_PATH)
find_library(TDEP_OLLE_LIBRARY NAMES olle)

find_library(TDEP_FLAP_LIBRARY NAMES flap PATHS ${TDEP_PATH} ${TDEP_PATH}/lib ${TDEP_PATH}/lib64 NO_DEFAULT_PATH)
find_library(TDEP_FLAP_LIBRARY NAMES flap)

find_path(TDEP_INCLUDE_DIR type_crystalstructure.mod ${TDEP_PATH} ${TDEP_PATH}/include ${TDEP_PATH}/inc/libolle NO_DEFAULT_PATH)
find_path(TDEP_INCLUDE_DIR type_crystalstructure.mod)

if(TDEP_OLLE_LIBRARY AND TDEP_FLAP_LIBRARY AND TDEP_INCLUDE_DIR)
	set(TDEP_FOUND TRUE)
	set(TDEP_LIBRARIES ${TDEP_OLLE_LIBRARY} ${TDEP_FLAP_LIBRARY})
endif()

if(TDEP_FOUND)
	if(NOT TDEP_FIND_QUIETLY)
		message(STATUS "Found TDEP: ${TDEP_LIBRARY}")
	endif()
else()
	if(TDEP_FIND_REQUIRED)
		if(NOT TDEP_OLLE_LIBRARY)
			message(FATAL_ERROR "Could not find TDEP libolle (Add -D TDEP_PATH=<path> to the cmake commandline for a non-standard installation)")
		endif()
		if(NOT TDEP_FLAP_LIBRARY)
			message(FATAL_ERROR "Could not find TDEP libflap (Add -D TDEP_PATH=<path> to the cmake commandline for a non-standard installation)")
		endif()
		if(NOT TDEP_INCLUDE_DIR)
			message(FATAL_ERROR "Could not find TDEP modules (Add -D TDEP_PATH=<path> to the cmake commandline for a non-standard installation)")
		endif()
	endif()
endif()
