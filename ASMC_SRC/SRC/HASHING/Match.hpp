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

#ifndef ASMC_HASHING_MATCH_HPP
#define ASMC_HASHING_MATCH_HPP

#include <algorithm>
#include <array>
#include <utility>

#include "HASHING/Utils.hpp"
#include "HMM.hpp"

/**
 * Match object that does \\todo not clear to me what this does...
 *
 * @tparam WORD_SIZE the word size
 */
template <int WORD_SIZE> class Match
{
private:
  std::array<int, 2> mInterval = {0, 0};
  unsigned mGaps = 0u;

public:
  explicit Match(const int i = 0) : mInterval{i, i}
  {
  }

  // pair : identifiers for the corresponding Individuals in all_ind
  void print(std::pair<unsigned, unsigned> p, const double PAR_MIN_MATCH, const std::vector<float>& geneticPositions,
             HMM& hmm)
  {
    double mlen = asmc::cmBetween(mInterval[0], mInterval[1], geneticPositions, WORD_SIZE);
    if (mlen >= PAR_MIN_MATCH) {
      const int from = mInterval[0] * WORD_SIZE;
      const int to = mInterval[1] * WORD_SIZE + WORD_SIZE - 1;
      hmm.decodeFromGERMLINE(p.first, p.second, from, to);
    }
  }

  void extend(const int w)
  {
    mInterval[1] = std::max<int>(w, mInterval[1]);
  }

  void addGap()
  {
    mGaps++;
  }

  const array<int, 2>& getInterval() const
  {
    return mInterval;
  }

  unsigned int getGaps() const
  {
    return mGaps;
  }
};

#endif // ASMC_HASHING_MATCH_HPP
