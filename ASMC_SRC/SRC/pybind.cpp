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

#include <pybind11/pybind11.h>
#include <vector>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "ASMC.hpp"
#include "Individual.hpp"
#include "HMM.hpp"
#include "Data.hpp"
#include "DecodingQuantities.hpp"
#include "DecodingParams.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<bool>)
PYBIND11_MAKE_OPAQUE(std::vector<float>)
//PYBIND11_MAKE_OPAQUE(std::vector <std::vector <float> >)
PYBIND11_MAKE_OPAQUE(std::vector<Individual>)
PYBIND11_MAKE_OPAQUE(std::vector<PairObservations>)
PYBIND11_MAKE_OPAQUE(std::unordered_map<float, std::vector<float>>)
PYBIND11_MAKE_OPAQUE(std::unordered_map<int, std::vector<float>>)
PYBIND11_MAKE_OPAQUE(DecodingQuantities)
PYBIND11_MAKE_OPAQUE(DecodingReturnValues)
PYBIND11_MAKE_OPAQUE(PairObservations)
PYBIND11_MAKE_OPAQUE(Individual)
PYBIND11_MAKE_OPAQUE(Data)
PYBIND11_MAKE_OPAQUE(HMM)

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(pyASMC, m) {
    py::enum_<DecodingModeOverall>(m, "DecodingModeOverall", py::arithmetic())
        .value("sequence", DecodingModeOverall::sequence)
        .value("array", DecodingModeOverall::array);
    py::enum_<DecodingMode>(m, "DecodingMode", py::arithmetic())
        .value("sequenceFolded", DecodingMode::sequenceFolded)
        .value("arrayFolded", DecodingMode::arrayFolded)
        .value("sequence", DecodingMode::sequence)
        .value("array", DecodingMode::array);
    py::bind_vector<std::vector<bool>>(m, "VectorBool");
    py::bind_vector<std::vector<float>>(m, "VectorFloat");
    py::bind_vector<std::vector<Individual>>(m, "VectorIndividual");
    py::bind_vector<std::vector<uint>>(m, "VectorUInt");
    py::bind_vector<std::vector<PairObservations>>(m, "VectorPairObservations");
    py::bind_vector<std::vector<std::vector<float>>>(m, "Matrix");
    py::bind_map<std::unordered_map<float, std::vector<float>>>(m, "UMapFloatToVectorFloat");
    py::bind_map<std::unordered_map<int, std::vector<float>>>(m, "UMapIntToVectorFloat");
    py::class_<DecodingReturnValues>(m, "DecodingReturnValues")
        .def_readwrite("sumOverPairs", &DecodingReturnValues::sumOverPairs)
        .def_readwrite("sumOverPairs00", &DecodingReturnValues::sumOverPairs00)
        .def_readwrite("sumOverPairs01", &DecodingReturnValues::sumOverPairs01)
        .def_readwrite("sumOverPairs11", &DecodingReturnValues::sumOverPairs11)
        .def_readwrite("sites", &DecodingReturnValues::sites)
        .def_readwrite("states", &DecodingReturnValues::states)
        .def_readwrite("siteWasFlippedDuringFolding", &DecodingReturnValues::siteWasFlippedDuringFolding);
    py::class_<Individual>(m, "Individual")
        .def(py::init<string, string, int>(),
                "famId"_a = "", "IId"_a = "", "numOfSites"_a = 0)
        .def("setGenotype", &Individual::setGenotype,
                "hap"_a, "pos"_a, "val"_a)
        .def_readwrite("famId", &Individual::famId)
        .def_readwrite("IId", &Individual::IId)
        .def_readwrite("name", &Individual::name)
        .def_readwrite("genotype1", &Individual::genotype1)
        .def_readwrite("genotype2", &Individual::genotype2)
        ;
    py::class_<PairObservations>(m, "PairObservations")
        .def_readwrite("iName", &PairObservations::iName)
        .def_readwrite("jName", &PairObservations::jName)
        .def_readwrite("obsBits", &PairObservations::obsBits)
        .def_readwrite("homMinorBits", &PairObservations::homMinorBits);
    py::class_<DecodingQuantities>(m, "DecodingQuantities")
        .def(py::init<const char*>())
        .def_readwrite("CSFSSamples", &DecodingQuantities::CSFSSamples)
        .def_readwrite("states", &DecodingQuantities::states)
        .def_readwrite("initialStateProb", &DecodingQuantities::initialStateProb)
        .def_readwrite("expectedTimes", &DecodingQuantities::expectedTimes)
        .def_readwrite("discretization", &DecodingQuantities::discretization)
        .def_readwrite("timeVector", &DecodingQuantities::timeVector)
        .def_readwrite("columnRatios", &DecodingQuantities::columnRatios)
        .def_readwrite("classicEmissionTable", &DecodingQuantities::classicEmissionTable)
        .def_readwrite("compressedEmissionTable", &DecodingQuantities::compressedEmissionTable)
        .def_readwrite("Dvectors", &DecodingQuantities::Dvectors)
        .def_readwrite("Bvectors", &DecodingQuantities::Bvectors)
        .def_readwrite("Uvectors", &DecodingQuantities::Uvectors)
        .def_readwrite("rowRatioVectors", &DecodingQuantities::rowRatioVectors)
        .def_readwrite("homozygousEmissionMap", &DecodingQuantities::homozygousEmissionMap)
        .def_readwrite("CSFSmap", &DecodingQuantities::CSFSmap)
        .def_readwrite("foldedCSFSmap", &DecodingQuantities::foldedCSFSmap)
        .def_readwrite("ascertainedCSFSmap", &DecodingQuantities::ascertainedCSFSmap)
        .def_readwrite("foldedAscertainedCSFSmap", &DecodingQuantities::foldedAscertainedCSFSmap)
        ;
    py::class_<DecodingParams>(m, "DecodingParams")
        .def(py::init<string, string, string,
                int, int, string,
                bool, bool, bool, bool,
                float, bool, bool, bool,
                string, bool, bool>(),
                "hapsFileRoot"_a, "decodingQuantFile"_a,
                "outFileRoot"_a = "",
                "jobs"_a = 1, "jobInd"_a = 1,
                "decodingModeString"_a = "array",
                "decodingSequence"_a = false,
                "usingCSFS"_a = true,
                "compress"_a = false,
                "useAncestral"_a = false,
                "skipCSFSdistance"_a = 0.f,
                "noBatches"_a = false,
                "doPosteriorSums"_a = false,
                "doPerPairPosteriorMean"_a = false,
                "expectedCoalTimesFile"_a = "",
                "withinOnly"_a = false,
                "doMajorMinorPosteriorSums"_a = false)
        .def_readwrite("hapsFileRoot", &DecodingParams::hapsFileRoot)
        .def_readwrite("decodingQuantFile", &DecodingParams::decodingQuantFile)
        .def_readwrite("outFileRoot", &DecodingParams::outFileRoot)
        .def_readwrite("jobs", &DecodingParams::jobs)
        .def_readwrite("jobInd", &DecodingParams::jobInd)
        .def_readwrite("decodingModeString", &DecodingParams::decodingModeString)
        .def_readwrite("decodingMode", &DecodingParams::decodingMode)
        .def_readwrite("decodingSequence", &DecodingParams::decodingSequence)
        .def_readwrite("foldData", &DecodingParams::foldData)
        .def_readwrite("usingCSFS", &DecodingParams::usingCSFS)
        .def_readwrite("compress", &DecodingParams::compress)
        .def_readwrite("useAncestral", &DecodingParams::useAncestral)
        .def_readwrite("skipCSFSdistance", &DecodingParams::skipCSFSdistance)
        .def_readwrite("noBatches", &DecodingParams::noBatches)
        .def_readwrite("doPosteriorSums", &DecodingParams::doPosteriorSums)
        .def_readwrite("doPerPairPosteriorMean", &DecodingParams::doPerPairPosteriorMean)
        .def_readwrite("expectedCoalTimesFile", &DecodingParams::expectedCoalTimesFile)
        .def_readwrite("withinOnly", &DecodingParams::withinOnly)
        .def_readwrite("doMajorMinorPosteriorSums", &DecodingParams::doMajorMinorPosteriorSums)
        ;

    py::class_<Data>(m, "Data")
        .def(py::init<std::string, int, int, bool, bool>(),
             "hapsFileRoot"_a, "sites"_a, "totalSamplesBound"_a,
             "foldToMinorAlleles"_a, "decodingUsesCSFS"_a)
        .def_static("countHapLines", &Data::countHapLines)
        .def_readwrite("FamIDList", &Data::FamIDList)
        .def_readwrite("IIDList", &Data::IIDList)
        .def_readwrite("famAndIndNameList", &Data::famAndIndNameList)
        .def_readwrite("individuals", &Data::individuals)
        .def_readwrite("sampleSize", &Data::sampleSize)
        .def_readwrite("haploidSampleSize", &Data::haploidSampleSize)
        .def_readwrite("sites", &Data::sites)
        .def_readwrite("totalSamplesBound", &Data::totalSamplesBound)
        .def_readwrite("decodingUsesCSFS", &Data::decodingUsesCSFS)
        .def_readwrite("geneticPositions", &Data::geneticPositions)
        .def_readwrite("physicalPositions", &Data::physicalPositions)
        .def_readwrite("siteWasFlippedDuringFolding", &Data::siteWasFlippedDuringFolding)
        .def_readwrite("recRateAtMarker", &Data::recRateAtMarker)
        .def_readwrite("undistinguishedCounts", &Data::undistinguishedCounts)
        ;
    py::class_<HMM>(m, "HMM")
        .def(py::init<Data&, const DecodingQuantities&, DecodingParams&, bool, int>())
        .def("decode", &HMM::decode)
        .def("decodeAll", &HMM::decodeAll)
        .def("decodeSummarize", &HMM::decodeSummarize)
        .def("getDecodingReturnValues", &HMM::getDecodingReturnValues)
        .def("decodePair", &HMM::decodePair)
        .def("decodePairs", &HMM::decodePairs)
        .def("getBatchBuffer", &HMM::getBatchBuffer)
        .def("finishDecoding", &HMM::finishDecoding)
        ;
    m.def("makePairObs", &makePairObs);
    m.def("asmc", &run, "Runs ASMC on HAPS files",
          "hapsFileRoot"_a, "decodingQuantFile"_a,
          "outFileRoot"_a = "", "mode"_a = DecodingModeOverall::array,
          "jobs"_a = 0, "jobInd"_a = 0,
          "skipCSFSdistance"_a = 0,
          "compress"_a = false, "useAncestral"_a = false,
          "doPosteriorSums"_a = true, "doMajorMinorPosteriorSums"_a = false);
}
