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

#include "lsst/utils/python.h"

#include "lsst/geom/Angle.h"

namespace py = pybind11;

namespace lsst {
namespace geom {

using PyAngle = py::class_<Angle>;
using PyAngleUnit = py::class_<AngleUnit>;

namespace {

template <typename OtherT>
void declareAngleComparisonOperators(PyAngle& cls) {
    cls.def("__eq__", [](Angle const& self, OtherT const& other) { return self == other; },
            py::is_operator());
    cls.def("__ne__", [](Angle const& self, OtherT const& other) { return self != other; },
            py::is_operator());
    cls.def("__le__", [](Angle const& self, OtherT const& other) { return self <= other; },
            py::is_operator());
    cls.def("__ge__", [](Angle const& self, OtherT const& other) { return self >= other; },
            py::is_operator());
    cls.def("__lt__", [](Angle const& self, OtherT const& other) { return self < other; }, py::is_operator());
    cls.def("__gt__", [](Angle const& self, OtherT const& other) { return self > other; }, py::is_operator());
}

} // anonymous

void wrapAngle(utils::python::WrapperCollection & wrappers) {
    wrappers.wrapType(
        PyAngleUnit(wrappers.module, "AngleUnit"),
        [](auto & mod, auto & cls) mutable {
            cls.def("__eq__", [](AngleUnit const& self, AngleUnit const& other) { return self == other; },
                    py::is_operator());
            cls.def("__ne__", [](AngleUnit const& self, AngleUnit const& other) { return !(self == other); },
                    py::is_operator());
            cls.def("_mul", [](AngleUnit const& self, double other) { return other * self; },
                    py::is_operator());
            cls.def("_rmul", [](AngleUnit const& self, double other) { return other * self; },
                    py::is_operator());
            mod.attr("radians") = py::cast(radians);
            mod.attr("degrees") = py::cast(degrees);
            mod.attr("hours") = py::cast(hours);
            mod.attr("arcminutes") = py::cast(arcminutes);
            mod.attr("arcseconds") = py::cast(arcseconds);
            mod.attr("milliarcseconds") = py::cast(milliarcseconds);
        }
    );

    wrappers.wrapType(
        PyAngle(wrappers.module, "Angle"),
        [](auto & mod, auto & cls) mutable {
            cls.def(py::init<double, AngleUnit>(), py::arg("val"), py::arg("units") = radians);
            cls.def(py::init<>());
            declareAngleComparisonOperators<Angle>(cls);
            declareAngleComparisonOperators<double>(cls);
            declareAngleComparisonOperators<int>(cls);
            cls.def("__mul__", [](Angle const& self, double other) { return self * other; },
                    py::is_operator());
            cls.def("__mul__", [](Angle const& self, int other) { return self * other; },
                    py::is_operator());
            cls.def("__rmul__", [](Angle const& self, double other) { return self * other; },
                    py::is_operator());
            cls.def("__rmul__", [](Angle const& self, int other) { return self * other; },
                    py::is_operator());
            cls.def("__imul__", [](Angle& self, double other) { return self *= other; });
            cls.def("__imul__", [](Angle& self, int other) { return self *= other; });
            cls.def("__add__", [](Angle const& self, Angle const& other) { return self + other; },
                    py::is_operator());
            cls.def("__sub__", [](Angle const& self, Angle const& other) { return self - other; },
                    py::is_operator());
            cls.def("__neg__", [](Angle const& self) { return -self; }, py::is_operator());
            cls.def("__iadd__", [](Angle& self, Angle const& other) { return self += other; });
            cls.def("__isub__", [](Angle& self, Angle const& other) { return self -= other; });
            cls.def("__truediv__", [](Angle const& self, double other) { return self / other; },
                    py::is_operator());
            // Without an explicit wrapper, Python lets Angle / Angle -> Angle
            cls.def("__truediv__", [](Angle const& self, Angle const& other) {
                throw py::type_error("unsupported operand type(s) for /: 'Angle' and 'Angle'");
            });
            cls.def("__float__", &Angle::operator double);
            cls.def("__abs__", [](Angle const& self) { return std::abs(self.asRadians()) * radians; });

            cls.def("__reduce__", [cls](Angle const& self) {
                return py::make_tuple(cls, py::make_tuple(py::cast(self.asRadians())));
            });
            utils::python::addOutputOp(cls, "__str__");
            utils::python::addOutputOp(cls, "__repr__");
            cls.def("asAngularUnits", &Angle::asAngularUnits);
            cls.def("asRadians", &Angle::asRadians);
            cls.def("asDegrees", &Angle::asDegrees);
            cls.def("asHours", &Angle::asHours);
            cls.def("asArcminutes", &Angle::asArcminutes);
            cls.def("asArcseconds", &Angle::asArcseconds);
            cls.def("asMilliarcseconds", &Angle::asMilliarcseconds);
            cls.def("wrap", &Angle::wrap);
            cls.def("wrapCtr", &Angle::wrapCtr);
            cls.def("wrapNear", &Angle::wrapNear);
            cls.def("separation", &Angle::separation);
            mod.def("isAngle", isAngle<Angle>);
            mod.def("isAngle", isAngle<double>);
        }
    );

    wrappers.wrapFunctions(
        [](auto & mod) mutable {
            mod.attr("PI") = py::float_(PI);
            mod.attr("TWOPI") = py::float_(TWOPI);
            mod.attr("HALFPI") = py::float_(HALFPI);
            mod.attr("ONE_OVER_PI") = py::float_(ONE_OVER_PI);
            mod.attr("SQRTPI") = py::float_(SQRTPI);
            mod.attr("INVSQRTPI") = py::float_(INVSQRTPI);
            mod.attr("ROOT2") = py::float_(ROOT2);
            mod.def("degToRad", degToRad);
            mod.def("radToDeg", radToDeg);
            mod.def("radToArcsec", radToArcsec);
            mod.def("radToMas", radToMas);
            mod.def("arcsecToRad", arcsecToRad);
            mod.def("masToRad", masToRad);
        }
    );
}

}  // namespace geom
}  // namespace lsst
