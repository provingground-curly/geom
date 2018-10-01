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
#include "pybind11/stl.h"

#include "lsst/geom/Box.h"
#include "lsst/utils/python.h"

namespace py = pybind11;
using namespace py::literals;

namespace lsst {
namespace geom {

void wrapBox(utils::python::WrapperCollection & wrappers) {
    wrappers.wrapType(
        py::class_<Box2I, std::shared_ptr<Box2I>>(wrappers.module, "Box2I"),
        [](auto & mod, auto & cls) mutable {

            cls.attr("Point") = mod.attr("Point2I");
            cls.attr("Extent") = mod.attr("Extent2I");

            py::enum_<Box2I::EdgeHandlingEnum>(cls, "EdgeHandlingEnum")
                    .value("EXPAND", Box2I::EdgeHandlingEnum::EXPAND)
                    .value("SHRINK", Box2I::EdgeHandlingEnum::SHRINK)
                    .export_values();

            cls.def(py::init<>());
            cls.def(py::init<Point2I const &, Point2I const &, bool>(), "minimum"_a, "maximum"_a,
                    "invert"_a = true);
            cls.def(py::init<Point2I const &, Extent2I const &, bool>(), "corner"_a, "dimensions"_a,
                    "invert"_a = true);
            cls.def(py::init<Box2D const &, Box2I::EdgeHandlingEnum>(), "other"_a,
                    "edgeHandling"_a = Box2I::EXPAND);
            cls.def(py::init<Box2I const &>(), "other"_a);

            cls.def("__eq__", [](Box2I const &self, Box2I const &other) { return self == other; },
                    py::is_operator());
            cls.def("__ne__", [](Box2I const &self, Box2I const &other) { return self != other; },
                    py::is_operator());

            cls.def_static("makeCenteredBox", &Box2I::makeCenteredBox, "center"_a, "size"_a);
            cls.def("swap", &Box2I::swap);
            cls.def("getMin", &Box2I::getMin);
            cls.def("getMinX", &Box2I::getMinX);
            cls.def("getMinY", &Box2I::getMinY);
            cls.def("getMax", &Box2I::getMax);
            cls.def("getMaxX", &Box2I::getMaxX);
            cls.def("getMaxY", &Box2I::getMaxY);
            cls.def("getBegin", &Box2I::getBegin);
            cls.def("getBeginX", &Box2I::getBeginX);
            cls.def("getBeginY", &Box2I::getBeginY);
            cls.def("getEnd", &Box2I::getEnd);
            cls.def("getEndX", &Box2I::getEndX);
            cls.def("getEndY", &Box2I::getEndY);
            cls.def("getDimensions", &Box2I::getDimensions);
            cls.def("getWidth", &Box2I::getWidth);
            cls.def("getHeight", &Box2I::getHeight);
            cls.def("getArea", &Box2I::getArea);
            cls.def("isEmpty", &Box2I::isEmpty);
            cls.def("contains", (bool (Box2I::*)(Point2I const &) const) & Box2I::contains);
            cls.def("contains", (bool (Box2I::*)(Box2I const &) const) & Box2I::contains);
            cls.def("__contains__", (bool (Box2I::*)(Point2I const &) const) & Box2I::contains);
            cls.def("__contains__", (bool (Box2I::*)(Box2I const &) const) & Box2I::contains);
            cls.def("overlaps", &Box2I::overlaps);
            cls.def("grow", (void (Box2I::*)(int)) & Box2I::grow);
            cls.def("grow", (void (Box2I::*)(Extent2I const &)) & Box2I::grow);
            cls.def("shift", &Box2I::shift);
            cls.def("flipLR", &Box2I::flipLR);
            cls.def("flipTB", &Box2I::flipTB);
            cls.def("include", (void (Box2I::*)(Point2I const &)) & Box2I::include);
            cls.def("include", (void (Box2I::*)(Box2I const &)) & Box2I::include);
            cls.def("clip", &Box2I::clip);
            cls.def("getCorners", &Box2I::getCorners);
            cls.def("toString", &Box2I::toString);
            cls.def("__repr__", [](Box2I const &self) {
                return py::str("Box2I(minimum={}, dimensions={})")
                    .format(py::repr(py::cast(self.getMin())), py::repr(py::cast(self.getDimensions())));
            });
            cls.def("__str__", [](Box2I const &self) {
                return py::str("(minimum={}, maximum={})")
                    .format(py::str(py::cast(self.getMin())), py::str(py::cast(self.getMax())));
            });
            cls.def("__reduce__", [cls](Box2I const &self) {
                return py::make_tuple(cls, make_tuple(py::cast(self.getMin()), py::cast(self.getMax())));
            });
            cls.def("getSlices", [](Box2I const &self) {
                return py::make_tuple(py::slice(self.getBeginY(), self.getEndY(), 1),
                                      py::slice(self.getBeginX(), self.getEndX(), 1));
            });

            mod.attr("BoxI") = cls;
        }
    );

    wrappers.wrapType(
        py::class_<Box2D, std::shared_ptr<Box2D>>(wrappers.module, "Box2D"),
        [](auto & mod, auto & cls) mutable {

            cls.attr("Point") = mod.attr("Point2D");
            cls.attr("Extent") = mod.attr("Extent2D");

            cls.attr("EPSILON") = py::float_(Box2D::EPSILON);
            cls.attr("INVALID") = py::float_(Box2D::INVALID);

            cls.def(py::init<>());
            cls.def(py::init<Point2D const &, Point2D const &, bool>(), "minimum"_a, "maximum"_a,
                    "invert"_a = true);
            cls.def(py::init<Point2D const &, Extent2D const &, bool>(), "corner"_a, "dimensions"_a,
                    "invert"_a = true);
            cls.def(py::init<Box2I const &>());
            cls.def(py::init<Box2D const &>());

            cls.def("__eq__", [](Box2D const &self, Box2D const &other) { return self == other; },
                    py::is_operator());
            cls.def("__ne__", [](Box2D const &self, Box2D const &other) { return self != other; },
                    py::is_operator());

            cls.def_static("makeCenteredBox", &Box2D::makeCenteredBox, "center"_a, "size"_a);
            cls.def("swap", &Box2D::swap);
            cls.def("getMin", &Box2D::getMin);
            cls.def("getMinX", &Box2D::getMinX);
            cls.def("getMinY", &Box2D::getMinY);
            cls.def("getMax", &Box2D::getMax);
            cls.def("getMaxX", &Box2D::getMaxX);
            cls.def("getMaxY", &Box2D::getMaxY);
            cls.def("getDimensions", &Box2D::getDimensions);
            cls.def("getWidth", &Box2D::getWidth);
            cls.def("getHeight", &Box2D::getHeight);
            cls.def("getArea", &Box2D::getArea);
            cls.def("getCenter", &Box2D::getCenter);
            cls.def("getCenterX", &Box2D::getCenterX);
            cls.def("getCenterY", &Box2D::getCenterY);
            cls.def("isEmpty", &Box2D::isEmpty);
            cls.def("contains", (bool (Box2D::*)(Point2D const &) const) & Box2D::contains);
            cls.def("contains", (bool (Box2D::*)(Box2D const &) const) & Box2D::contains);
            cls.def("__contains__", (bool (Box2D::*)(Point2D const &) const) & Box2D::contains);
            cls.def("__contains__", (bool (Box2D::*)(Box2D const &) const) & Box2D::contains);
            cls.def("overlaps", &Box2D::overlaps);
            cls.def("grow", (void (Box2D::*)(double)) & Box2D::grow);
            cls.def("grow", (void (Box2D::*)(Extent2D const &)) & Box2D::grow);
            cls.def("overlaps", (void (Box2D::*)(Extent2D const &)) & Box2D::overlaps);
            cls.def("shift", &Box2D::shift);
            cls.def("flipLR", &Box2D::flipLR);
            cls.def("flipTB", &Box2D::flipTB);
            cls.def("include", (void (Box2D::*)(Point2D const &)) & Box2D::include);
            cls.def("include", (void (Box2D::*)(Box2D const &)) & Box2D::include);
            cls.def("clip", &Box2D::clip);
            cls.def("getCorners", &Box2D::getCorners);
            cls.def("toString", &Box2D::toString);
            cls.def("__repr__", [](Box2D const &self) {
                return py::str("Box2D(minimum={}, dimensions={})")
                    .format(py::repr(py::cast(self.getMin())), py::repr(py::cast(self.getDimensions())));
            });
            cls.def("__str__", [](Box2D const &self) {
                return py::str("(minimum={}, maximum={})")
                    .format(py::str(py::cast(self.getMin())), py::str(py::cast(self.getMax())));
            });
            cls.def("__reduce__", [cls](Box2D const &self) {
                return py::make_tuple(cls, make_tuple(py::cast(self.getMin()), py::cast(self.getMax())));
            });

            mod.attr("BoxD") = cls;
        }
    );
}

}  // namespace geom
}  // namespace lsst
