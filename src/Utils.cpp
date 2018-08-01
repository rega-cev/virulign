#include "Utils.h"

#include <string>
#include <algorithm>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

std::string to_upper_copy(const std::string& s)
{
  std::string copy = s;
  std::transform(copy.begin(), copy.end(), copy.begin(), ::toupper);
  return copy;
}

bool ends_with(const std::string& s, const std::string& p)
{
  if (p.size() > s.size())
    return false;
  else {
    for (unsigned i = 0; i < p.size(); ++i) {
      if (p[i] != s[s.size() - p.size() + i])
	return false;
    }
    return true;
  }
}

long long current_time_ms()
{
#ifdef _WIN32
  static LARGE_INTEGER s_frequency;
  static BOOL s_use_qpc = QueryPerformanceFrequency(&s_frequency);
  //if there is no high resolution time stamp, use GetTickCount()
  if (s_use_qpc) {
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);
    return (1000LL * now.QuadPart) / s_frequency.QuadPart;
  } else {
    return GetTickCount();
  }  
#else
  struct timeval  tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ; 
#endif
}

std::string format_time(const long long& milliseconds)
{
  int seconds = (int) (milliseconds / 1000) % 60 ;
  int minutes = (int) ((milliseconds / (1000*60)) % 60);
  int hours   = (int) ((milliseconds / (1000*60*60)) % 24);

  std::stringstream ss;
  if (hours != 0)
    ss << hours << "h";
  if (minutes != 0 || hours != 0)
    ss << minutes << "m";
  ss << seconds << "s";

  return ss.str();
}
