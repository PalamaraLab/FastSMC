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

PYBIND11_MODULE(asmc, m) {
    py::enum_<DecodingMode>(m, "DecodingMode", py::arithmetic())
        .value("sequenceFolded", DecodingMode::sequenceFolded)
        .value("arrayFolded", DecodingMode::arrayFolded)
        .value("sequence", DecodingMode::sequence)
        .value("array", DecodingMode::array);
    py::class_<DecodingReturnValues>(m, "DecodingReturnValues");
    py::class_<DecodingQuantities>(m, "DecodingQuantities");
    py::class_<Data>(m, "Data");
    m.def("asmc", &run, "Runs ASMC on HAPS files");
}
