# Checks if the argument value is one of possible values contained in list
# This version is case SENSITIVE
# allowed_values ... list of values that are allowed
# result         ... result of the check (true or false)

function(check_if_arg_value_is_valid_strict allowed_values checked_value result)
  set(${result} FALSE PARENT_SCOPE)
  foreach(value ${allowed_values})
    if ( ${value} STREQUAL ${checked_value} )
      set(${result} TRUE PARENT_SCOPE)
    endif()
  endforeach(value)
endfunction(check_if_arg_value_is_valid_strict)

# Checks if the argument value is one of possible values contained in list
# This version is case INSENSITIVE
# allowed_values ... list of values that are allowed
# result         ... result of the check (true or false)

function(check_if_arg_value_is_valid allowed_values checked_value result)
  set(${result} FALSE PARENT_SCOPE)
  string( TOLOWER "${checked_value}" checked_value_lower )
  foreach(value ${allowed_values})
    string( TOLOWER "${value}" value_lower )
    if ( ${value_lower} STREQUAL ${checked_value_lower} )
      set(${result} TRUE PARENT_SCOPE)
    endif()
  endforeach(value)
endfunction(check_if_arg_value_is_valid)
