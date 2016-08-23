#ifndef PARMOON__UTILITIES
#define PARMOON__UTILITIES

#include <string>

namespace utilities
{

/// @brief return the host name of the computer ParMooN is running on
std::string get_host_name();

/// @brief return the date and time as a string
std::string get_date_and_time();

}

#endif// PARMOON__UTILITIES