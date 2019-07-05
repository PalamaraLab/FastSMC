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


#ifndef ASMC_HPP
#define ASMC_HPP

#include "DecodingQuantities.hpp"
#include "Data.hpp"
#include "DecodingParams.hpp"

DecodingReturnValues run(std::string haps_file_root, std::string decoding_quant_file,
         std::string out_file_root, DecodingModeOverall mode,
         int jobs, int job_index,
         float skip_csfs_distance,
         bool compress, bool use_ancestral,
         bool posterior_sums, bool major_minor_posterior_sums);

#endif
