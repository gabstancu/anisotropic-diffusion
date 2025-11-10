#ifndef HELPER_HPP
#define HELPER_HPP

#include <cstdlib> // for std::getenv
#include <filesystem>

inline std::string get_home_directory () 
{
    const char* home = std::getenv("HOME");
    return home ? std::string(home) : "";
}


inline std::string get_current_working_directory() 
{
    return std::filesystem::current_path().string();
}


#endif // HELPER_HPP
