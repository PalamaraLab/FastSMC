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

#ifndef ASMC_FASTSMC_HPP
#define ASMC_FASTSMC_HPP

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "HMM.hpp"

namespace ASMC
{

class FastSMC
{

private:
  int mHashingWordSize = 64;
  int mConstReadAhead = 10;
  bool mHaploid = true;

public:

  FastSMC() = default;

  void run(const DecodingParams& params, const Data& data, HMM& hmm);
};

} // namespace ASMC

#endif // ASMC_FASTSMC_HPP
