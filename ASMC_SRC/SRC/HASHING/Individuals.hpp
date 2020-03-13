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

#ifndef ASMC_HASHING_INDIVIDUALS_HPP
#define ASMC_HASHING_INDIVIDUALS_HPP

#include <array>
#include <bitset>
#include <cassert>
#include <string>

template <int WORD_SIZE, int CONST_READ_AHEAD> class Individuals
{
  unsigned mIdNum;
  std::array<std::bitset<WORD_SIZE>, CONST_READ_AHEAD> mHap;

public:
  explicit Individuals(const unsigned idNum) : mIdNum(idNum)
  {
    for (auto w = 0; w < CONST_READ_AHEAD; w++) {
      clear(w);
    }
  }

  void clear(const int w)
  {
    assert(w >= 0);
    mHap[w % CONST_READ_AHEAD].reset();
  }

  void setMarker(const int w, const std::size_t bit)
  {
    assert(w >= 0);
    assert(bit < WORD_SIZE);
    mHap[w % CONST_READ_AHEAD].set(bit);
  }

  unsigned long getWordHash(const int w)
  {
    assert(w >= 0);
    return mHap[w % CONST_READ_AHEAD].to_ulong();
  }

  std::string getWordString(const int w)
  {
    assert(w >= 0);
    return mHap[w % CONST_READ_AHEAD].to_string();
  }

  unsigned int getIdNum() const
  {
    return mIdNum;
  }
};

#endif // ASMC_HASHING_INDIVIDUALS_HPP
