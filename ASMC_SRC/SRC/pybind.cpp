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
#include "ASMC.hpp"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(asmc, m) {
    py::enum_<DecodingModeOverall>(m, "DecodingModeOverall", py::arithmetic())
        .value("sequence", DecodingModeOverall::sequence)
        .value("array", DecodingModeOverall::array);
    py::class_<DecodingReturnValues>(m, "DecodingReturnValues")
        .def_readwrite("sumOverPairs", &DecodingReturnValues::sumOverPairs)
        .def_readwrite("sumOverPairs00", &DecodingReturnValues::sumOverPairs00)
        .def_readwrite("sumOverPairs01", &DecodingReturnValues::sumOverPairs01)
        .def_readwrite("sumOverPairs11", &DecodingReturnValues::sumOverPairs11)
        .def_readwrite("sites", &DecodingReturnValues::sites)
        .def_readwrite("states", &DecodingReturnValues::states)
        .def_readwrite("siteWasFlippedDuringFolding", &DecodingReturnValues::siteWasFlippedDuringFolding);
    m.def("asmc", &run, "Runs ASMC on HAPS files",
          "haps_file_root"_a, "decoding_quant_file"_a,
          "out_file_root"_a, "mode"_a,
          "jobs"_a = 0, "job_index"_a = 0,
          "skip_csfs_distance"_a = 0,
          "compress"_a = false, "use_ancestral"_a = false);
}
