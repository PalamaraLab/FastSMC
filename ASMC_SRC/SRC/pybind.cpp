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

PYBIND11_MAKE_OPAQUE(std::vector<bool>)
PYBIND11_MAKE_OPAQUE(std::vector<float>)
PYBIND11_MAKE_OPAQUE(std::vector <std::vector <float> >)

#include "ASMC.hpp"
#include "Individual.hpp"
#include "HMM.hpp"
#include "Data.hpp"
#include "DecodingQuantities.hpp"
#include "DecodingParams.hpp"

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
    py::bind_vector<std::vector<std::vector<float>>>(m, "Matrix");
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
                py::arg("_famId") = "", py::arg("_IId") = "",
                py::arg("numOfSites") = 0)
        .def("setGenotype", &Individual::setGenotype);
    py::class_<PairObservations>(m, "PairObservations")
        .def_readwrite("iName", &PairObservations::iName)
        .def_readwrite("jName", &PairObservations::jName)
        .def_readwrite("obsBits", &PairObservations::obsBits)
        .def_readwrite("homMinorBits", &PairObservations::homMinorBits);
    py::class_<DecodingQuantities>(m, "DecodingQuantities")
        .def(py::init<const char*>())
        .def_readwrite("CSFSSamples", &DecodingQuantities::CSFSSamples)
        ;
    py::class_<DecodingParams>(m, "DecodingParams")
        .def(py::init<string, string, string,
                int, int, string, DecodingMode,
                bool, bool, bool, bool, bool,
                float, bool, bool, bool,
                string, bool, bool>(),
                "hapsFileRoot"_a, "decodingQuantFile"_a,
                "outFileRoot"_a = "",
                "jobs"_a = 1, py::arg("jobInd") = 1,
                "decodingModeString"_a = "array",
                "decodingMode"_a = DecodingMode::arrayFolded,
                "decodingSequence"_a = false,
                "foldData"_a = true,
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
        .def_static("countHapLines", &Data::countHapLines);
    py::class_<HMM>(m, "HMM")
        .def(py::init<Data&, const DecodingQuantities&, DecodingParams&, bool, int>())
        .def("decode", &HMM::decode);
    m.def("asmc", &run, "Runs ASMC on HAPS files",
          "hapsFileRoot"_a, "decodingQuantFile"_a,
          "outFileRoot"_a = "", "mode"_a = DecodingModeOverall::array,
          "jobs"_a = 0, "jobInd"_a = 0,
          "skipCSFSdistance"_a = 0,
          "compress"_a = false, "useAncestral"_a = false,
          "doPosteriorSums"_a = true, "doMajorMinorPosteriorSums"_a = false);
}
