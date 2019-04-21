/*
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "pybind11/stl.h"

#include "ndarray/pybind11.h"

#include "lsst/geom/AffineTransform.h"
#include "lsst/utils/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace geom {

void wrapAffineTransform(utils::python::WrapperCollection & wrappers) {
    wrappers.wrapType(
        py::class_<AffineTransform, std::shared_ptr<AffineTransform>>(wrappers.module, "AffineTransform"),
        [](auto & mod, auto & cls) mutable {

            // Parameters enum is really only used as integer constants.
            cls.attr("XX") = py::cast(int(AffineTransform::Parameters::XX));
            cls.attr("YX") = py::cast(int(AffineTransform::Parameters::YX));
            cls.attr("XY") = py::cast(int(AffineTransform::Parameters::XY));
            cls.attr("YY") = py::cast(int(AffineTransform::Parameters::YY));
            cls.attr("X") = py::cast(int(AffineTransform::Parameters::X));
            cls.attr("Y") = py::cast(int(AffineTransform::Parameters::Y));

            /* Constructors */
            cls.def(py::init<>());
            cls.def(py::init<Eigen::Matrix3d const &>(), "matrix"_a);
            cls.def(py::init<Eigen::Matrix2d const &>(), "linear"_a);
            cls.def(py::init<Eigen::Vector2d const &>(), "translation"_a);
            cls.def(py::init<Eigen::Matrix2d const &, Eigen::Vector2d const &>(),
                    "linear"_a, "translation"_a);
            cls.def(py::init<LinearTransform const &>(), "linear"_a);
            cls.def(py::init<Extent2D const &>(), "translation"_a);
            cls.def(py::init<LinearTransform const &, Extent2D const &>(), "linear"_a, "translation"_a);

            /* Operators and special methods */
            cls.def("__mul__", &AffineTransform::operator*, py::is_operator());
            cls.def("__call__",
                    py::overload_cast<Point2D const &>(&AffineTransform::operator(), py::const_));
            cls.def("__call__",
                    py::overload_cast<Extent2D const &>(&AffineTransform::operator(), py::const_));
            cls.def("__call__",
                    [](py::object self, py::object x, py::object y) mutable {
                        return py::make_tuple(self.attr("applyX")(x, y),
                                              self.attr("applyY")(x, y));
                    },
                    "x"_a, "y"_a);
            cls.def("__setitem__", [](AffineTransform &self, int i, double value) {
                if (i < 0 || i > 5) {
                    PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d", i);
                    throw py::error_already_set();
                }
                self[i] = value;
            });
            cls.def("__getitem__", [](AffineTransform const &self, int row, int col) {
                if (row < 0 || row > 2 || col < 0 || col > 2) {
                    PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d, %d", row, col);
                    throw py::error_already_set();
                }
                return (self.getMatrix())(row, col);
            });
            cls.def("__getitem__", [](AffineTransform const &self, int i) {
                if (i < 0 || i > 5) {
                    PyErr_Format(PyExc_IndexError, "Invalid index for AffineTransform: %d", i);
                    throw py::error_already_set();
                }
                return self[i];
            });
            cls.def("__str__", [](AffineTransform const &self) {
                return py::str(py::cast(self.getMatrix())); }
            );
            cls.def("__repr__", [](AffineTransform const &self) {
                return py::str("AffineTransform(\n{}\n)").format(py::cast(self.getMatrix()));
            });
            cls.def("__reduce__", [cls](AffineTransform const &self) {
                return py::make_tuple(cls, py::make_tuple(py::cast(self.getMatrix())));
            });

            /* Members */
            cls.def("inverted", &AffineTransform::inverted);
            cls.def("isIdentity", &AffineTransform::isIdentity);
            cls.def("getTranslation", (Extent2D & (AffineTransform::*)()) & AffineTransform::getTranslation);
            cls.def("getLinear", (LinearTransform & (AffineTransform::*)()) & AffineTransform::getLinear);
            cls.def("getMatrix", &AffineTransform::getMatrix);
            cls.def("getParameterVector", &AffineTransform::getParameterVector);
            cls.def("setParameterVector", &AffineTransform::setParameterVector);
            cls.def("applyX", py::vectorize(&AffineTransform::applyX), "x"_a, "y"_a);
            cls.def("applyY", py::vectorize(&AffineTransform::applyY), "x"_a, "y"_a);
            cls.def_static("makeScaling", py::overload_cast<double>(&AffineTransform::makeScaling));
            cls.def_static("makeScaling", py::overload_cast<double, double>(&AffineTransform::makeScaling));
            cls.def_static("makeRotation", &AffineTransform::makeRotation, "angle"_a);
            cls.def_static("makeTranslation", &AffineTransform::makeTranslation, "translation"_a);

            /* Non-members */
            mod.def("makeAffineTransformFromTriple", makeAffineTransformFromTriple);
        }
    );
}

}  // namespace geom
}  // namespace lsst
