//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.

#include "HmmUtils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace asmc
{

std::vector<bool> subsetXorVec(const std::vector<bool>& v1, const std::vector<bool>& v2, const unsigned long from,
                               const unsigned long to) noexcept
{
  const auto min_to = std::min<unsigned long>(v1.size(), to);

  assert(v1.size() == v2.size());
  assert(from < min_to);

  std::vector<bool> ret(min_to - from);

  for (unsigned long i = from; i < min_to; ++i) {
    ret[i - from] = v1[i] ^ v2[i];
  }

  return ret;
}

std::vector<bool> subsetAndVec(const std::vector<bool>& v1, const std::vector<bool>& v2, const unsigned long from,
                               const unsigned long to) noexcept
{
  const auto min_to = std::min<unsigned long>(v1.size(), to);

  assert(v1.size() == v2.size());
  assert(from < min_to);

  std::vector<bool> ret(min_to - from);

  for (unsigned i = from; i < min_to; i++) {
    ret[i - from] = v1[i] & v2[i];
  }

  return ret;
}

float roundMorgans(const float value, const int precision, const float min) noexcept
{
  assert(precision >= 0);
  assert(min > 0.f);

  if (value <= min) {
    return min;
  }

  const float correction = 10.f - static_cast<float>(precision);
  const float L10 = std::max<float>(0.f, floorf(log10f(value)) + correction);
  const float factor = powf(10.f, 10.f - L10);

  return roundf(value * factor) / factor;
}

int roundPhysical(const int value, const int precision) noexcept
{
  assert(precision >= 0);
  assert(value > -2); // Since HMM for sequence uses distance-1, it can be -1

  if (value <= 1) {
    return 1;
  }

  int L10 = std::max<int>(0, static_cast<int>(floor(log10(value))) - precision);
  int factor = static_cast<int>(pow(10, L10));

  return static_cast<int>(round(value / static_cast<double>(factor))) * factor;
}

} // namespace asmc
