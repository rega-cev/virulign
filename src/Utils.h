// This may look like C code, but it's really -*- C++ -*-
#ifndef UTILS_H_ 
#define UTILS_H_ 

#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>

template <typename T>
T lexical_cast(const std::string& s)
{
    std::stringstream ss(s);

    T result;
    if ((ss >> result).fail() || !(ss >> std::ws).eof())
    {
        throw std::bad_cast();
    }

    return result;
}

template <typename T>
std::string to_string(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::string to_upper_copy(const std::string& s);

bool ends_with(const std::string& s, const std::string& p);

#endif // UTILS_H_ 
