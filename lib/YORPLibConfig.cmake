if (EXISTS /home/bebe0705/.am_fortuna)
	set(IS_FORTUNA ON)
	message("-- This is Fortuna")

else()
	set(IS_FORTUNA OFF)
endif()


if(${IS_FORTUNA})
	set(YORPLIB_INCLUDE_HEADER /home/bebe0705/libs/local/include/YORPLib/)
	set(YORPLIB_LIBRARY /home/bebe0705/libs/local/lib/libYORPLIB.so)

else()
	set(YORPLIB_INCLUDE_HEADER /usr/local/include/YORPLib/)

	if (APPLE)
		set(YORPLIB_LIBRARY /usr/local/lib/libYORPLIB.dylib)
	elseif(UNIX AND NOT APPLE)
		set(YORPLIB_LIBRARY /usr/local/lib/libYORPLIB.so)
	else()
		message(FATAL_ERROR "Unsupported platform")
	endif()
endif()


message("-- Found YORPLib: " ${YORPLIB_LIBRARY})
