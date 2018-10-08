##############################################################################
# this macro logs simultaneously to screen and to file
##############################################################################
macro( pdekit_log line )
   message(STATUS ${line})
   file(APPEND ${PROJECT_LOG_FILE} "${line}\n")
endmacro()

##############################################################################
# this macro logs just to file
##############################################################################
macro( pdekit_log_file line )
   file(APPEND ${PROJECT_LOG_FILE} "${line}\n")
endmacro()
##############################################################################

##############################################################################
# this macro logs to screen if we are verbose
##############################################################################
macro( pdekit_log_verbose line )
    if( PDEKIT_CMAKE_VERBOSE )
      pdekit_log( ${line} )
    else()
      pdekit_log_file( ${line} )
    endif()
endmacro()
##############################################################################


