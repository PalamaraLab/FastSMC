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

  bool m_storeFullPosteriors = false;
  bool m_storeSumOfPosteriors = false;
  bool m_storePerPairPosteriors = false;
  bool m_storePerPairMAPs = false;

  std::size_t numWritten = 0ul;

  /// The pair indices: iInd, iHap, jInd, jHap
  Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajor> m_perPairIndices;

  /// The full set of posteriors: for each pair this is a (states * numSites) matrix
  Eigen::Array<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Dynamic, 1> m_perPairPosteriors;

  /// The sum of all posteriors in perPairPosteriors: a (states * numSites) matrix
  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_sumOfPosteriors;

  /// Posterior means: each row is an array of length numSites
  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_perPairPosteriorMeans;

  Eigen::Array<float, 1, Eigen::Dynamic, Eigen::RowMajor> m_minPosteriorMeans;
  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> m_argminPosteriorMeans;

  Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_perPairMAPs;

  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> m_minMAPs;
  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> m_argminMAPs;

public:
  void initialise(const std::vector<uint>& individualsA, const std::vector<uint>& individualsB, long int numSites,
                  long int numStates, bool fullPosteriors = false, bool sumOfPosteriors = false,
                  bool perPairPosteriors = false, bool perPairMAPs = false)
  {
    m_storeFullPosteriors = fullPosteriors;
    m_storeSumOfPosteriors = sumOfPosteriors;
    m_storePerPairPosteriors = perPairPosteriors;
    m_storePerPairMAPs = perPairMAPs;

    Eigen::Index numPairsToDecode = 0ll;
    for (auto idx = 0ul; idx < individualsA.size(); ++idx) {
      if (individualsA[idx] == individualsB[idx]) {
        numPairsToDecode += 1ll;
      } else {
        numPairsToDecode += 4ll;
      }
    }

    m_perPairIndices.resize(numPairsToDecode, Eigen::NoChange);

    if (m_storeFullPosteriors) {
      m_perPairPosteriors.resize(numPairsToDecode);
      for (Eigen::Index i = 0ll; i < m_perPairPosteriors.size(); ++i) {
        m_perPairPosteriors(i).resize(numStates, numSites);
      }
    }

    if (m_storeSumOfPosteriors) {
      m_sumOfPosteriors.resize(numStates, numSites);
      m_sumOfPosteriors.setZero();
    }

    if (m_storePerPairPosteriors) {
      m_perPairPosteriorMeans.resize(numPairsToDecode, numSites);
      m_minPosteriorMeans.resize(numSites);
      m_argminPosteriorMeans.resize(numSites);
    }

    if (m_storePerPairMAPs) {
      m_perPairMAPs.resize(numPairsToDecode, numSites);
      m_minMAPs.resize(numSites);
      m_argminMAPs.resize(numSites);
    }
  }

  void incrementNumWritten()
  {
    numWritten += 1;
  }

  void finaliseCalculations()
  {
    for (Eigen::Index siteIdx = 0ll; siteIdx < m_perPairPosteriorMeans.cols(); ++siteIdx) {
      Eigen::Index argmin{};
      m_minPosteriorMeans(siteIdx) = m_perPairPosteriorMeans.col(siteIdx).minCoeff(&argmin);
      m_argminPosteriorMeans(siteIdx) = static_cast<int>(argmin);
    }

    for (Eigen::Index siteIdx = 0ll; siteIdx < m_perPairMAPs.cols(); ++siteIdx) {
      Eigen::Index argmin{};
      m_minMAPs(siteIdx) = m_perPairMAPs.col(siteIdx).minCoeff(&argmin);
      m_argminMAPs(siteIdx) = static_cast<int>(argmin);
    }
  }

  [[nodiscard]] const Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajor>& getPerPairIndices() const
  {
    return m_perPairIndices;
  }

  [[nodiscard]] const Eigen::Array<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Dynamic, 1>& getPerPairPosteriors() const
  {
    return m_perPairPosteriors;
  }

  [[nodiscard]] const Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getSumOfPosteriors() const
  {
    return m_sumOfPosteriors;
  }

  [[nodiscard]] const Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getPerPairPosteriorMeans() const
  {
    return m_perPairPosteriorMeans;
  }

  [[nodiscard]] const Eigen::Array<float, 1, Eigen::Dynamic, Eigen::RowMajor>& getMinPosteriorMeans() const
  {
    return m_minPosteriorMeans;
  }

  [[nodiscard]] const Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>& getArgminPosteriorMeans() const
  {
    return m_argminPosteriorMeans;
  }

  [[nodiscard]] const Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getPerPairMAPs() const
  {
    return m_perPairMAPs;
  }

  [[nodiscard]] const Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>& getMinMAPs() const
  {
    return m_minMAPs;
  }

  [[nodiscard]] const Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>& getArgminMAPs() const
  {
    return m_argminMAPs;
  }

  [[nodiscard]] Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajor>& getModifiablePerPairIndices()
  {
    return m_perPairIndices;
  }

  [[nodiscard]] Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getModifiableSumOfPosteriors()
  {
    return m_sumOfPosteriors;
  }

  [[nodiscard]] Eigen::Array<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Dynamic, 1>& getModifiablePerPairPosteriors()
  {
    return m_perPairPosteriors;
  }

  [[nodiscard]] Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getModifiablePerPairPosteriorMeans()
  {
    return m_perPairPosteriorMeans;
  }

  [[nodiscard]] Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getModifiablePerPairMAPs()
  {
    return m_perPairMAPs;
  }

  [[nodiscard]] std::size_t getNumWritten() const
  {
    return numWritten;
  }
};

#endif // FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP
