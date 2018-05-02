#include "Utils.h"

#include <string>
#include <algorithm>

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
