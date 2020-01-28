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

#ifndef HMMUTILS_HPP
#define HMMUTILS_HPP

#include <limits>
#include <vector>

namespace asmc
{

/**
 * Return the elementwise xor of two vectors, or optionally the elementwise xor of the subset of two vectors between
 * indices `from` and `to`. If `from` and `to` are left as defaults, the returned vector is the full elementwise xor,
 * the same length as each input.
 *
 * It is assumed that both input vectors are the same length and that, if specified, from < to.
 *
 * This is used for decoding as xor --> 1 if heterozygous.
 *
 * @param v1 the first vector
 * @param v2 the second vector
 * @param from the starting index of a subset (default 0)
 * @param to the ending index of a subset (default equivalent to v1.size())
 * @return a vector of length v1.size() (or to - from) with the elementwise xor of v1 and v2 (or a subset thereof).
 */
std::vector<bool> subsetXorVec(const std::vector<bool>& v1, const std::vector<bool>& v2, unsigned long from = 0ul,
                               unsigned long to = std::numeric_limits<unsigned long>::max()) noexcept;

/**
 * Return the elementwise and of two vectors, or optionally the elementwise and of the subset of two vectors between
 * indices `from` and `to`. If `from` and `to` are left as defaults, the returned vector is the full elementwise and,
 * the same length as each input.
 *
 * It is assumed that both input vectors are the same length and that, if specified, from < to.
 *
 * This is used for decoding to distinguish homozygous minor/derived from homozygous major/ancestral.
 *
 * @param v1 the first vector
 * @param v2 the second vector
 * @param from the starting index of a subset (default 0)
 * @param to the ending index of a subset (default equivalent to v1.size())
 * @return a vector of length v1.size() (or to - from) with the elementwise and of v1 and v2 (or a subset thereof).
 */
std::vector<bool> subsetAndVec(const std::vector<bool>& v1, const std::vector<bool>& v2, unsigned long from = 0ul,
                               unsigned long to = std::numeric_limits<unsigned long>::max()) noexcept;

} // namespace asmc

#endif // HMMUTILS_HPP
