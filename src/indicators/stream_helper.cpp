
/*
Activity Indicators for Modern C++
https://github.com/p-ranav/indicators

Licensed under the MIT License <http://opensource.org/licenses/MIT>.
SPDX-License-Identifier: MIT
Copyright (c) 2019 Dawid Pilarski <dawid.pilarski@panicsoftware.com>.

Permission is hereby  granted, free of charge, to any  person obtaining a copy
of this software and associated  documentation files (the "Software"), to deal
in the Software  without restriction, including without  limitation the rights
to  use, copy,  modify, merge,  publish, distribute,  sublicense, and/or  sell
copies  of  the Software,  and  to  permit persons  to  whom  the Software  is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE  IS PROVIDED "AS  IS", WITHOUT WARRANTY  OF ANY KIND,  EXPRESS OR
IMPLIED,  INCLUDING BUT  NOT  LIMITED TO  THE  WARRANTIES OF  MERCHANTABILITY,
FITNESS FOR  A PARTICULAR PURPOSE AND  NONINFRINGEMENT. IN NO EVENT  SHALL THE
AUTHORS  OR COPYRIGHT  HOLDERS  BE  LIABLE FOR  ANY  CLAIM,  DAMAGES OR  OTHER
LIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE  OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "display_width.cpp"
#include "setting.cpp"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

#include <cassert>
#include <cmath>

namespace indicators {
namespace details {

inline std::ostream &write_duration(std::ostream &os, std::chrono::nanoseconds ns) {
  using namespace std;
  using namespace std::chrono;
  using days = duration<int, ratio<86400>>;
  char fill = os.fill();
  os.fill('0');
  auto d = duration_cast<days>(ns);
  ns -= d;
  auto h = duration_cast<hours>(ns);
  ns -= h;
  auto m = duration_cast<minutes>(ns);
  ns -= m;
  auto s = duration_cast<seconds>(ns);
  if (d.count() > 0)
    os << setw(2) << d.count() << "d:";
  if (h.count() > 0)
    os << setw(2) << h.count() << "h:";
  os << setw(2) << m.count() << "m:" << setw(2) << s.count() << 's';
  os.fill(fill);
  return os;
}

class BlockProgressScaleWriter {
public:
  BlockProgressScaleWriter(std::ostream &os, size_t bar_width) : os(os), bar_width(bar_width) {}

  std::ostream &write(float progress) {
    std::string fill_text{"█"};
    std::vector<std::string> lead_characters{" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉"};
    auto value = std::min(1.0f, std::max(0.0f, progress / 100.0f));
    auto whole_width = std::floor(value * bar_width);
    auto remainder_width = fmod((value * bar_width), 1.0f);
    auto part_width = std::floor(remainder_width * lead_characters.size());
    std::string lead_text = lead_characters[size_t(part_width)];
    if ((bar_width - whole_width - 1) < 0)
      lead_text = "";
    for (size_t i = 0; i < whole_width; ++i)
      os << fill_text;
    os << lead_text;
    for (size_t i = 0; i < (bar_width - whole_width - 1); ++i)
      os << " ";
    return os;
  }

private:
  std::ostream &os;
  size_t bar_width = 0;
};

class ProgressScaleWriter {
public:
  ProgressScaleWriter(std::ostream &os, size_t bar_width, const std::string &fill,
                      const std::string &lead, const std::string &remainder)
      : os(os), bar_width(bar_width), fill(fill), lead(lead), remainder(remainder) {}

  std::ostream &write(float progress) {
    auto pos = static_cast<size_t>(progress * bar_width / 100.0);
    for (size_t i = 0, current_display_width = 0; i < bar_width;) {
      std::string next;

      if (i < pos) {
        next = fill;
        current_display_width = unicode::display_width(fill);
      } else if (i == pos) {
        next = lead;
        current_display_width = unicode::display_width(lead);
      } else {
        next = remainder;
        current_display_width = unicode::display_width(remainder);
      }

      i += current_display_width;

      if (i > bar_width) {
        // `next` is larger than the allowed bar width
        // fill with empty space instead
        os << std::string((bar_width - (i - current_display_width)), ' ');
        break;
      }

      os << next;
    }
    return os;
  }

private:
  std::ostream &os;
  size_t bar_width = 0;
  std::string fill;
  std::string lead;
  std::string remainder;
};

class IndeterminateProgressScaleWriter {
public:
  IndeterminateProgressScaleWriter(std::ostream &os, size_t bar_width, const std::string &fill,
                                   const std::string &lead)
      : os(os), bar_width(bar_width), fill(fill), lead(lead) {}

  std::ostream &write(size_t progress) {
    for (size_t i = 0; i < bar_width;) {
      std::string next;
      size_t current_display_width = 0;

      if (i < progress) {
        next = fill;
        current_display_width = unicode::display_width(fill);
      } else if (i == progress) {
        next = lead;
        current_display_width = unicode::display_width(lead);
      } else {
        next = fill;
        current_display_width = unicode::display_width(fill);
      }

      i += current_display_width;

      if (i > bar_width) {
        // `next` is larger than the allowed bar width
        // fill with empty space instead
        os << std::string((bar_width - (i - current_display_width)), ' ');
        break;
      }

      os << next;
    }
    return os;
  }

private:
  std::ostream &os;
  size_t bar_width = 0;
  std::string fill;
  std::string lead;
};

} // namespace details
} // namespace indicators
