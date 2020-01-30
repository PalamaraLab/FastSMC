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

#include <cassert>
#include <iostream>
#include <limits>
#include <numeric>
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

/**
 * Round a float to an integer number of multiples of 10^-10. The number of significant digits in the rounded value is
 * controlled by the precision, with (precision + 1) being the number of significant digits in the rounded value.
 *
 * This is very similar to just rounding to (precision + 1) significant figures.
 *
 * This function never returns a value lower than the given minimum value.
 *
 * For example:
 *  - 0.1 will be unchanged no matter the precision.
 *  - 0.77426 will round to 0.8 (with precision 0), to 0.77 (with precision 1), to 0.774 (with precision 2), etc.
 *
 * @param value the value to round
 * @param precision one less than the number of significant figures to keep
 * @param min the minimum allowed value after rounding
 * @return the rounded value
 */
float roundMorgans(float value, int precision, float min) noexcept;

/**
 * Round an integer to (precision + 1) significant figures.
 *
 * This function never returns a value smaller than 1.
 *
 * For example:
 *  - -1, 0 and 1 will all round to 1
 *  - 1234 will round to 1000 (with precision 0), to 1200 (with precision 1), to 1230 (with precision 2), etc.
 *
 * @param value the value to round
 * @param precision
 * @return the rounded value
 */
int roundPhysical(int value, int precision) noexcept;

/**
 * Print a tab-separated vector, to std::cout by default
 *
 * @tparam T the vector value type
 * @param v the vector to print
 * @param os the ostream to send output to, default is std::cout
 */
template <typename T> void printVector(const std::vector<T>& v, std::ostream& os = std::cout)
{
  for (const auto& x : v) {
    os << '\t' << x;
  }
  os << '\n';
}

/**
 * Helper to print a percentage of time spent in certain functions.
 *
 * @param str description of the function (forward, backward etc)
 * @param fracTime fraction of time spent in that function
 * @param os the ostream to send output to, default is std::cout
 */
void printPctTime(const std::string& str, double fracTime, std::ostream& os = std::cout);

/**
 * Calculate the sum of elements in a vector
 *
 * @tparam Summable the vector value type
 * @param v the vector to sum
 * @return the sum of the elements of the vector
 */
template <typename Summable> Summable getSumOfVector(const std::vector<Summable>& v)
{
  return std::accumulate(v.cbegin(), v.cend(), Summable{0});
}

/**
 * Calculate the element-wise product of a vector with a scalar and return the new vector. The existing vector remains
 * unaltered.
 *
 * @tparam RealType the vector and scalar value type
 * @param vec a vector
 * @param val a scalar
 * @return a vector the element-wise multiplication of vec with val
 */
template <typename RealType>
std::vector<RealType> elementWiseMultVectorScalar(const std::vector<RealType>& vec, RealType val)
{
  std::vector<RealType> ret;
  ret.reserve(vec.size());

  for (const auto& x : vec) {
    ret.push_back(x * val);
  }

  return ret;
}

/**
 * Calculate the element-wise product of a vector with a another vector and return the new vector. The existing vectors
 * remain unaltered.
 *
 * @tparam RealType the value type of both vectors
 * @param v1 the first vector
 * @param v2 the second vector
 * @return a vector the element-wise multiplication of v1 with v2
 */
template <typename RealType>
std::vector<RealType> elementWiseMultVectorVector(const std::vector<RealType>& v1, const std::vector<RealType>& v2)
{
  assert(v1.size() == v2.size());

  std::vector<RealType> ret;
  ret.reserve(v1.size());

  for (auto i = 0ul; i < v1.size(); ++i) {
    ret.push_back(v1[i] * v2[i]);
  }

  return ret;
}

/**
 * Calculate the element-wise product of a matrix with a another matrix and return the new matrix. The existing matrices
 * remain unaltered.
 *
 * @tparam RealType the value type of both matrices
 * @param mat1 the first matrix
 * @param mat2 the second matrix
 * @return a matrix the element-wise multiplication of mat1 with mat2
 */
template <typename RealType>
std::vector<std::vector<RealType>> elementWiseMultMatrixMatrix(const std::vector<std::vector<RealType>>& mat1,
                                                               const std::vector<std::vector<RealType>>& mat2)
{
  assert(mat1.size() == mat2.size());
  assert(mat1.size() > 0);
  assert(mat1.at(0).size() == mat2.at(0).size());

  std::vector<std::vector<RealType>> ret(mat1.size(), std::vector<RealType>(mat1[0].size()));
  for (auto i = 0ul; i < mat1.size(); ++i) {
    for (auto j = 0ul; j < mat1[0].size(); ++j) {
      ret[i][j] = mat1[i][j] * mat2[i][j];
    }
  }
  return ret;
}

/**
 * Normalize each column of a matrix and return the new matrix. The existing matrix remains unaltered.
 *
 * It is assumed that no column of the matrix sums to zero.
 *
 * @tparam RealType the value type of the matrix
 * @param matrix the matrix
 * @return a matrix with normalised rows (that sum to 1)
 */
template <typename RealType>
std::vector<std::vector<RealType>> normalizeMatrixColumns(const std::vector<std::vector<RealType>>& matrix)
{
  assert(matrix.size() > 0);
  std::vector<std::vector<RealType>> ret(matrix.size(), std::vector<RealType>(matrix[0].size()));

  for (auto colIdx = 0ul; colIdx < matrix[0].size(); ++colIdx) {

    auto sum = RealType{0};
    for (const auto& row : matrix) {
      sum += row[colIdx];
    }
    assert(sum != RealType{0});

    for (auto rowIdx = 0ul; rowIdx < matrix.size(); ++rowIdx) {
      ret[rowIdx][colIdx] = matrix[rowIdx][colIdx] / sum;
    }
  }
  return ret;
}

/**
 * Fill a matrix column with a vector, modifying the original matrix.
 *
 * @tparam RealType the value type of the matrix and vector
 * @param matrix the matrix to modify a column of
 * @param vec the vector fill a column of the matrix
 * @param colIdx the index of the column to fill
 */
template <typename RealType>
void fillMatrixColumn(std::vector<std::vector<RealType>>& matrix, const std::vector<RealType>& vec,
                      unsigned long colIdx)
{
  assert(matrix.size() == vec.size());
  assert(matrix.size() > 0);
  assert(colIdx < matrix.at(0).size());

  for (auto rowIdx = 0ul; rowIdx < vec.size(); ++rowIdx) {
    matrix[rowIdx][colIdx] = vec[rowIdx];
  }
}

/**
 * Calculate scaling factors for a contiguous array of curBatchSize * numStates floats. The computed scaling factors
 * are normalised across all states for each item in the batch.
 *
 * This method takes three buffers; one with data, one for storing calculated scaling factors, and one for temporarily
 * storing intermediate sums. It is assumed that all three buffers are appropriately allocated and deallocated outside
 * this function.
 *
 * The sums buffer is set to zeros in this method.
 *
 * @param vec the buffer of length curBatchSize * numStates containing the data
 * @param scalings a buffer of length curBatchSize to write calculated scaling factors into
 * @param sums a buffer of length curBatchSize for storing intermediate sums
 * @param batchSize the number of items in the batch
 * @param numStates the number of states over which to normalise
 */
void calculateScalingBatch(float* vec, float* scalings, float* sums, int batchSize, int numStates);

/**
 * Apply scaling factors to a contiguous array of curBatchSize * numStates floats. The scale factors are applied per
 * item to each state, and each state is separated in the data array by a stride length of batchSize.
 *
 * This method takes two buffers; one with data that is modified in-place, and one with scale factors that is just read
 * from.
 *
 * @param vec the buffer of length curBatchSize * numStates containing the data
 * @param scalings a buffer of length curBatchSize containing the prescribed scale factors
 * @param batchSize the number of items in the batch
 * @param numStates the number of states over which to normalise
 */
void applyScalingBatch(float* vec, float* scalings, int batchSize, int numStates);

} // namespace asmc

#endif // HMMUTILS_HPP
