#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

__all__ = ["CoordinateExpr", "Extent", "ExtentI", "ExtentD",
           "Point", "PointI", "PointD"]

from lsst.utils import TemplateMeta

from . import _geom


def _coordinateStr(self):
    return "({})".format(", ".join("%0.5g" % v for v in self))


def _coordinateRepr(self):
    return "{}({})".format(type(self).__name__,
                           ", ".join("%0.10g" % v for v in self))


def _coordinateReduce(self):
    return (type(self), tuple(self))


class CoordinateExpr(metaclass=TemplateMeta):
    """Abstract base class and factory for CoordinateExpr objects.
    """
    TEMPLATE_PARAMS = ("dimensions", )

    __str__ = _coordinateStr
    __repr__ = _coordinateRepr
    __reduce__ = _coordinateReduce


CoordinateExpr.register(2, _geom.CoordinateExpr2)
CoordinateExpr.register(3, _geom.CoordinateExpr3)


class Extent(metaclass=TemplateMeta):
    """Abstract base class and factory for Extent objects.
    """
    TEMPLATE_PARAMS = ("dtype", "dimensions")
    TEMPLATE_DEFAULTS = (None, 2)

    __str__ = _coordinateStr
    __repr__ = _coordinateRepr
    __reduce__ = _coordinateReduce


Extent.register((int, 2), _geom.Extent2I)
Extent.register((float, 2), _geom.Extent2D)
Extent.register((int, 3), _geom.Extent3I)
Extent.register((float, 3), _geom.Extent3D)
ExtentI = _geom.Extent2I
ExtentD = _geom.Extent2D


class Point(metaclass=TemplateMeta):
    """Abstract base class and factory for Point objects.
    """
    TEMPLATE_PARAMS = ("dtype", "dimensions")
    TEMPLATE_DEFAULTS = (None, 2)

    __str__ = _coordinateStr
    __repr__ = _coordinateRepr
    __reduce__ = _coordinateReduce


Point.register((int, 2), _geom.Point2I)
Point.register((float, 2), _geom.Point2D)
Point.register((int, 3), _geom.Point3I)
Point.register((float, 3), _geom.Point3D)
PointI = _geom.Point2I
PointD = _geom.Point2D
