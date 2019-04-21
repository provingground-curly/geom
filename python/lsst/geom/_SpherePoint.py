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

__all__ = []

from lsst.utils import continueClass

from ._geom import SpherePoint, _toUnitX, _toUnitY, _toUnitZ


def _pickExactlyOne(a, b, message):
    if a is None and b is not None:
        return b
    if a is not None and b is None:
        return a
    raise ValueError(message)


@continueClass  # noqa F811
class SpherePoint:

    def __iter__(self):
        for i in (0, 1):
            yield self[i]

    def __repr__(self):
        argList = ["%r*geom.degrees" % (pos.asDegrees(),) for pos in self]
        return "SpherePoint(%s)" % (", ".join(argList))

    @staticmethod
    def toUnitXYZ(self, *, longitude=None, latitude=None, ra=None, dec=None, units):
        """Compute the unit 3-vectors (as separate arrays) corresponding to
        arrays of longitude and latitude.

        Parameters
        ----------
        longitude : `float` or `numpy.ndarray`
            Longitude coordinate of input points.
        latitude : `float` or `numpy.ndarray`
            Latitude coordinate of input points.
        ra : `float` or `numpy.ndarray`
            Synonym for `longitude`.
        dec : `float` or `numpy.ndarray`
            Synonym for `latitude`.
        units : `AngleUnit`
            Angle unit for inputs.

        Returns
        -------
        x : `float` or numpy.ndarray`
            X coordinates of unit 3-vectors.
        y : `float` or numpy.ndarray`
            Y coordinates of unit 3-vectors.
        z : `float` or numpy.ndarray`
            Z coordinates of unit 3-vectors.
        """
        factor = (1.0*units).asRadians()
        lon = factor*_pickExactlyOne(longitude, ra, "Exactly one of ra and longitude must be provided.")
        lat = factor*_pickExactlyOne(latitude, dec, "Exactly one of dec and latitude must be provided.")
        return _toUnitX(lon, lat), _toUnitY(lon, lat), _toUnitZ(lon, lat)
