##############################################################################
# adds a path to search for when searching for include files
##############################################################################
macro( pdekit_add_trial_include_path INCPATH )
 if( EXISTS ${INCPATH})
  list( APPEND TRIAL_INCLUDE_PATHS ${INCPATH} )
 endif( EXISTS ${INCPATH})
endmacro( pdekit_add_trial_include_path )
##############################################################################

##############################################################################
# adds a path to search for when searching for library files
##############################################################################
macro( pdekit_add_trial_library_path LIBPATH )
 if( EXISTS ${LIBPATH})
   list( APPEND TRIAL_LIBRARY_PATHS ${LIBPATH} )
 endif( EXISTS ${LIBPATH})
endmacro( pdekit_add_trial_library_path )
##############################################################################

##############################################################################
# sets a path to search for when searching for INCLUDE files
##############################################################################
macro( pdekit_set_trial_include_path INCPATHS )
 set( TRIAL_INCLUDE_PATHS "" )
 foreach( path ${INCPATHS} )
   if( EXISTS ${path})
     list( APPEND TRIAL_INCLUDE_PATHS ${path} )
   endif( EXISTS ${path})
 endforeach( path ${INCPATHS} )
endmacro( pdekit_set_trial_include_path )
##############################################################################

##############################################################################
# sets a path to search for when searching for library files
##############################################################################
macro( pdekit_set_trial_library_path LIBPATHS )
 set( TRIAL_LIBRARY_PATHS "" )
 foreach( path ${LIBPATHS} )
   if( EXISTS ${path})
     list( APPEND TRIAL_LIBRARY_PATHS ${path} )
   endif( EXISTS ${path})
 endforeach( path ${LIBPATHS} )
endmacro( pdekit_set_trial_library_path )
##############################################################################

