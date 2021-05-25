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

#ifndef FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP
#define FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP

#include <Eigen/Core>
#include <fmt/format.h>

#include <vector>

struct DecodePairsReturnStruct {

private:
  bool inUse = false;

  std::size_t numWritten = 0ul;

  Eigen::Array<unsigned int, Eigen::Dynamic, 4, Eigen::RowMajor> indices;
  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> posteriors;
  Eigen::Array<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MAPs;

public:
  void resize(const std::vector<uint>& individualsA, const std::vector<uint>& individualsB, long int numSites)
  {
    Eigen::Index resizeTo = 0ll;
    for (auto idx = 0ul; idx < individualsA.size(); ++idx) {
      if (individualsA[idx] == individualsB[idx]) {
        resizeTo += 1ll;
      } else {
        resizeTo += 4ll;
      }
    }

    indices.resize(resizeTo, Eigen::NoChange);
    posteriors.resize(resizeTo, numSites);
    inUse = true;
  }

  void incrementNumWritten()
  {
    numWritten += 1;
  }

  [[nodiscard]] const Eigen::Array<unsigned int, Eigen::Dynamic, 4, Eigen::RowMajor>& getIndices() const
  {
    return indices;
  }

  [[nodiscard]] const Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getPosteriors() const
  {
    return posteriors;
  }

  [[nodiscard]] const Eigen::Array<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getMAPs() const
  {
    return MAPs;
  }

  [[nodiscard]] Eigen::Array<unsigned int, Eigen::Dynamic, 4, Eigen::RowMajor>& getModifiableIndices()
  {
    return indices;
  }

  [[nodiscard]] Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getModifiablePosteriors()
  {
    return posteriors;
  }

  [[nodiscard]] Eigen::Array<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getModifiableMAPs()
  {
    return MAPs;
  }

  [[nodiscard]] std::size_t getNumWritten() const
  {
    return numWritten;
  }

  [[nodiscard]] bool isInUse() const
  {
    return inUse;
  }
};

#endif // FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP
