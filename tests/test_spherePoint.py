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

# -*- python -*-
"""
Unit tests for SpherePoint

Run with:
   python testSpherePoint.py
or
   python
   >>> import testSpherePoint
   >>> testSpherePoint.run()
"""

import copy
import math
import re
import unittest

import numpy as np
from numpy.testing import assert_allclose

import lsst.utils.tests
import lsst.sphgeom
import lsst.geom as geom
import lsst.pex.exceptions as pexEx

from lsst.geom import degrees, radians, SpherePoint
from numpy import nan, inf


class SpherePointTestSuite(lsst.utils.tests.TestCase):

    def setUp(self):
        self._dataset = SpherePointTestSuite.positions()
        self._poleLatitudes = [
            geom.HALFPI*geom.radians,
            6.0*geom.hours,
            90.0*geom.degrees,
            5400.0*geom.arcminutes,
            324000.0*geom.arcseconds,
        ]

    @property
    def pointSet(self):
        for lon, lat in self._dataset:
            for point in (
                SpherePoint(lon, lat),
                SpherePoint(lon.asDegrees(), lat.asDegrees(), degrees),
                SpherePoint(lon.asRadians(), lat.asRadians(), radians),
            ):
                yield point

    @staticmethod
    def positions():
        """Provide valid coordinates for nominal-case testing.

        Returns
        -------
        positions : `iterable`
            An iterable of pairs of Angles, each representing the
            longitude and latitude (in that order) of a test point.
        """
        nValidPoints = 100
        rng = np.random.RandomState(42)
        ra = rng.uniform(0.0, 360.0, nValidPoints)
        dec = rng.uniform(-90.0, 90.0, nValidPoints)

        points = list(zip(ra*degrees, dec*degrees))
        # Ensure corner cases are tested.
        points += [
            (0.0*degrees, 0.0*degrees),
            (geom.PI*radians, -6.0*degrees),
            (42.0*degrees, -90.0*degrees),
            (172.0*degrees, geom.HALFPI*radians),
            (360.0*degrees, 45.0*degrees),
            (-278.0*degrees, -42.0*degrees),
            (765.0*degrees, 0.25*geom.PI*radians),
            (180.0*degrees, nan*radians),
            (inf*degrees, 45.0*degrees),
            (nan*degrees, -8.3*degrees),
        ]
        return points

    def testLonLatConstructorErrors(self):
        """Test if the longitude, latitude constructors handle invalid input
        """
        # Latitude should be checked for out-of-range.
        for lat in self._poleLatitudes:
            with self.assertRaises(pexEx.InvalidParameterError):
                SpherePoint(0.0*degrees, self.nextUp(lat))
            with self.assertRaises(pexEx.InvalidParameterError):
                SpherePoint(0.0, self.nextUp(lat).asDegrees(), degrees)
            with self.assertRaises(pexEx.InvalidParameterError):
                SpherePoint(0.0*degrees, self.nextDown(-lat))
            with self.assertRaises(pexEx.InvalidParameterError):
                SpherePoint(0.0, self.nextDown(-lat).asDegrees(), degrees)

        # Longitude should not be checked for out of range.
        SpherePoint(360.0*degrees, 45.0*degrees)
        SpherePoint(360.0, 45.0, degrees)
        SpherePoint(-42.0*degrees, 45.0*degrees)
        SpherePoint(-42.0, 45.0, degrees)
        SpherePoint(391.0*degrees, 45.0*degrees)
        SpherePoint(391.0, 45.0, degrees)

        # Infinite latitude is not allowed.
        with self.assertRaises(pexEx.InvalidParameterError):
            SpherePoint(-42.0*degrees, inf*degrees)
        with self.assertRaises(pexEx.InvalidParameterError):
            SpherePoint(-42.0, inf, degrees)
        with self.assertRaises(pexEx.InvalidParameterError):
            SpherePoint(-42.0*degrees, -inf*degrees)
        with self.assertRaises(pexEx.InvalidParameterError):
            SpherePoint(-42.0, -inf, degrees)

    def testVector3dConstructor(self):
        # test poles
        for z in (-11.3, -1.1, 0.1, 2.5):  # arbitrary non-zero values
            sp = SpherePoint(lsst.sphgeom.Vector3d(0.0, 0.0, z))
            self.assertTrue(sp.atPole())
            self.assertEqual(sp.getLongitude().asRadians(), 0.0)
            if z < 0:
                self.assertAnglesAlmostEqual(sp.getLatitude(), -90 * degrees)
            else:
                self.assertAnglesAlmostEqual(sp.getLatitude(), 90 * degrees)

        spx = SpherePoint(lsst.sphgeom.Vector3d(11.1, 0.0, 0.0))
        self.assertAnglesAlmostEqual(spx.getLongitude(), 0.0 * degrees)
        self.assertAnglesAlmostEqual(spx.getLatitude(), 0.0 * degrees)

        spy = SpherePoint(lsst.sphgeom.Vector3d(0.0, 234234.5, 0.0))
        self.assertAnglesAlmostEqual(spy.getLongitude(), 90.0 * degrees)
        self.assertAnglesAlmostEqual(spy.getLatitude(), 0.0 * degrees)

        spxy = SpherePoint(lsst.sphgeom.Vector3d(7.5, -7.5, 0.0))
        self.assertAnglesAlmostEqual(spxy.getLongitude(), -45.0 * degrees)
        self.assertAnglesAlmostEqual(spxy.getLatitude(), 0.0 * degrees)

        spxz = SpherePoint(lsst.sphgeom.Vector3d(100.0, 0.0, -100.0))
        self.assertAnglesAlmostEqual(spxz.getLongitude(), 0.0 * degrees)
        self.assertAnglesAlmostEqual(spxz.getLatitude(), -45.0 * degrees)

        # Only one singularity: a vector of all zeros
        with self.assertRaises(pexEx.InvalidParameterError):
            SpherePoint(lsst.sphgeom.Vector3d(0.0, 0.0, 0.0))

    def testDefaultConstructor(self):
        sp = SpherePoint()
        self.assertTrue(math.isnan(sp.getLongitude()))
        self.assertTrue(math.isnan(sp.getLatitude()))
        self.assertFalse(sp.isFinite())

    def testCopyConstructor(self):
        sp = SpherePoint(-42.0*degrees, 45.0*degrees)
        spcopy = SpherePoint(sp)
        self.assertEqual(sp, spcopy)

    def testInitNArgFail(self):
        """Test incorrect calls to the SpherePoint constructor
        """
        with self.assertRaises(TypeError):
            SpherePoint("Rotund", "Bovine")
        with self.assertRaises(TypeError):
            SpherePoint(42)
        with self.assertRaises(TypeError):
            SpherePoint("ICRS", 34.0, -56.0)
        with self.assertRaises(TypeError):
            SpherePoint(34.0, -56.0)  # missing units

    def testGetLongitudeValue(self):
        """Test if getLongitude() and getRa() return the expected value.
        """
        for lon, lat in self._dataset:
            for point in (
                SpherePoint(lon, lat),
                SpherePoint(lon.asDegrees(), lat.asDegrees(), degrees),
                SpherePoint(lon.asRadians(), lat.asRadians(), radians),
            ):
                self.assertIsInstance(point.getLongitude(), geom.Angle)
                # Behavior for non-finite points is undefined; depends on internal
                # data representation
                if point.isFinite():
                    self.assertGreaterEqual(point.getLongitude().asDegrees(), 0.0)
                    self.assertLess(point.getLongitude().asDegrees(), 360.0)

                    # Longitude not guaranteed to match input at pole
                    if not point.atPole():
                        # assertAnglesAlmostEqual handles angle wrapping internally
                        self.assertAnglesAlmostEqual(lon, point.getLongitude())
                        self.assertAnglesAlmostEqual(lon, point.getRa())

        # Vector construction should return valid longitude even in edge cases.
        point = SpherePoint(lsst.sphgeom.Vector3d(0.0, 0.0, -1.0))
        self.assertGreaterEqual(point.getLongitude().asDegrees(), 0.0)
        self.assertLess(point.getLongitude().asDegrees(), 360.0)

    def testGetPosition(self):
        for sp in self.pointSet:
            for units in (degrees, geom.hours, radians):
                point = sp.getPosition(units)
                expectedPoint = [val.asAngularUnits(units) for val in sp]
                assert_allclose(point, expectedPoint, atol=1e-15)

    def testTicket1394(self):
        """Regression test for Ticket 1761.

        Checks that negative longitudes within epsilon of lon=0 lead
        are correctly bounded and rounded.
        """
        # The problem was that the coordinate is less than epsilon
        # close to RA == 0 and bounds checking was getting a
        # negative RA.
        point = SpherePoint(lsst.sphgeom.Vector3d(
            0.6070619982, -1.264309928e-16, 0.7946544723))

        self.assertEqual(point[0].asDegrees(), 0.0)

    def testGetLatitudeValue(self):
        """Test if getLatitude() and getDec() return the expected value.
        """
        for lon, lat in self._dataset:
            for point in (
                SpherePoint(lon, lat),
                SpherePoint(lon.asDegrees(), lat.asDegrees(), degrees),
                SpherePoint(lon.asRadians(), lat.asRadians(), radians),
            ):
                self.assertIsInstance(point.getLatitude(), geom.Angle)
                # Behavior for non-finite points is undefined; depends on internal
                # data representation
                if point.isFinite():
                    self.assertGreaterEqual(point.getLatitude().asDegrees(), -90.0)
                    self.assertLessEqual(point.getLatitude().asDegrees(), 90.0)
                    self.assertAnglesAlmostEqual(lat, point.getLatitude())
                    self.assertAnglesAlmostEqual(lat, point.getDec())

    def testGetVectorValue(self):
        """Test if getVector() returns the expected value.

        The test includes conformance to vector-angle conventions.
        """
        for lon, lat, vector in [
            (0.0*degrees, 0.0*degrees, lsst.sphgeom.Vector3d(1.0, 0.0, 0.0)),
            (90.0*degrees, 0.0*degrees, lsst.sphgeom.Vector3d(0.0, 1.0, 0.0)),
            (0.0*degrees, 90.0*degrees, lsst.sphgeom.Vector3d(0.0, 0.0, 1.0)),
        ]:
            for point in (
                SpherePoint(lon, lat),
                SpherePoint(lon.asDegrees(), lat.asDegrees(), degrees),
                SpherePoint(lon.asRadians(), lat.asRadians(), radians),
            ):
                newVector = point.getVector()
                self.assertIsInstance(newVector, lsst.sphgeom.UnitVector3d)
                for oldElement, newElement in zip(vector, newVector):
                    self.assertAlmostEqual(oldElement, newElement)

                # Convert back to spherical.
                newLon, newLat = SpherePoint(newVector)
                self.assertAlmostEqual(newLon.asDegrees(), lon.asDegrees())
                self.assertAlmostEqual(newLat.asDegrees(), lat.asDegrees())

        # Try some un-normalized ones, too.
        pointList = [
            ((0.0, 0.0), lsst.sphgeom.Vector3d(1.3, 0.0, 0.0)),
            ((90.0, 0.0), lsst.sphgeom.Vector3d(0.0, 1.2, 0.0)),
            ((0.0, 90.0), lsst.sphgeom.Vector3d(0.0, 0.0, 2.3)),
            ((0.0, 0.0), lsst.sphgeom.Vector3d(0.5, 0.0, 0.0)),
            ((90.0, 0.0), lsst.sphgeom.Vector3d(0.0, 0.7, 0.0)),
            ((0.0, 90.0), lsst.sphgeom.Vector3d(0.0, 0.0, 0.9)),
        ]

        for lonLat, vector in pointList:
            # Only convert from vector to spherical.
            point = SpherePoint(vector)
            newLon, newLat = point
            self.assertAlmostEqual(lonLat[0], newLon.asDegrees())
            self.assertAlmostEqual(lonLat[1], newLat.asDegrees())
            vector = lsst.sphgeom.Vector3d(point.getVector())
            self.assertAlmostEqual(1.0, vector.getSquaredNorm())

        # Ill-defined points should be all NaN after normalization
        cleanValues = [0.5, -0.3, 0.2]
        badValues = [nan, inf, -inf]
        for i in range(3):
            for badValue in badValues:
                values = cleanValues[:]
                values[i] = badValue
                nonFiniteVector = lsst.sphgeom.Vector3d(*values)
                for element in SpherePoint(nonFiniteVector).getVector():
                    self.assertTrue(math.isnan(element))

    def testToUnitXZY(self):
        """Test that the numpy-vectorized transformation from (lat, lon) to
        (x, y, z) matches SpherePoint.getVector().
        """
        for units in (degrees, radians):
            scale = float(180.0*degrees)/float(1.0*units)
            lon = scale*np.random.rand(5, 3)
            lat = scale*(np.random.rand(5, 3) - 0.5)
        x, y, z = SpherePoint.toUnitXYZ(longitude=lon, latitude=lat, units=units)
        for i in range(lon.shape[0]):
            for j in range(lon.shape[1]):
                s = SpherePoint(lon[i, j], lat[i, j], units)
                u1 = s.getVector()
                u2 = lsst.sphgeom.UnitVector3d(x=x[i, j], y=y[i, j], z=z[i, j])
                self.assertFloatsAlmostEqual(np.array(u1, dtype=float), np.array(u2, dtype=float))

    def testTicket1761(self):
        """Regression test for Ticket 1761.

        Checks for math errors caused by unnormalized vectors.
        """
        refPoint = SpherePoint(lsst.sphgeom.Vector3d(0, 1, 0))

        point1 = SpherePoint(lsst.sphgeom.Vector3d(0.1, 0.1, 0.1))
        point2 = SpherePoint(lsst.sphgeom.Vector3d(0.6, 0.6, 0.6))
        sep1 = refPoint.separation(point1)
        sep2 = refPoint.separation(point2)
        sepTrue = 54.735610317245339*degrees

        self.assertAnglesAlmostEqual(sepTrue, sep1)
        self.assertAnglesAlmostEqual(sepTrue, sep2)

    def testAtPoleValue(self):
        """Test if atPole() returns the expected value.
        """
        poleList = \
            [SpherePoint(42.0*degrees, lat) for lat in self._poleLatitudes] + \
            [SpherePoint(42.0, lat.asDegrees(), degrees) for lat in self._poleLatitudes] + \
            [SpherePoint(42.0*degrees, -lat) for lat in self._poleLatitudes] + \
            [SpherePoint(42.0, -lat.asDegrees(), degrees) for lat in self._poleLatitudes] + \
            [
                SpherePoint(lsst.sphgeom.Vector3d(0.0, 0.0, 1.0)),
                SpherePoint(lsst.sphgeom.Vector3d(0.0, 0.0, -1.0)),
            ]
        nonPoleList = \
            [SpherePoint(42.0*degrees, self.nextDown(lat)) for lat in self._poleLatitudes] + \
            [SpherePoint(42.0, self.nextDown(lat).asDegrees(), degrees) for lat in self._poleLatitudes] + \
            [SpherePoint(42.0*degrees, self.nextUp(-lat)) for lat in self._poleLatitudes] + \
            [SpherePoint(42.0, self.nextUp(-lat).asDegrees(), degrees)
             for lat in self._poleLatitudes] + \
            [
                SpherePoint(lsst.sphgeom.Vector3d(9.9e-7, 0.0, 1.0)),
                SpherePoint(lsst.sphgeom.Vector3d(9.9e-7, 0.0, -1.0)),
                SpherePoint(0.0*degrees, nan*degrees),
            ]

        for pole in poleList:
            self.assertIsInstance(pole.atPole(), bool)
            self.assertTrue(pole.atPole())

        for nonPole in nonPoleList:
            self.assertIsInstance(nonPole.atPole(), bool)
            self.assertFalse(nonPole.atPole())

    def testIsFiniteValue(self):
        """Test if isFinite() returns the expected value.
        """
        finiteList = [
            SpherePoint(0.0*degrees, -90.0*degrees),
            SpherePoint(0.0, -90.0, degrees),
            SpherePoint(lsst.sphgeom.Vector3d(0.1, 0.2, 0.3)),
        ]
        nonFiniteList = [
            SpherePoint(0.0*degrees, nan*degrees),
            SpherePoint(0.0, nan, degrees),
            SpherePoint(nan*degrees, 0.0*degrees),
            SpherePoint(nan, 0.0, degrees),
            SpherePoint(inf*degrees, 0.0*degrees),
            SpherePoint(inf, 0.0, degrees),
            SpherePoint(-inf*degrees, 0.0*degrees),
            SpherePoint(-inf, 0.0, degrees),
            SpherePoint(lsst.sphgeom.Vector3d(nan, 0.2, 0.3)),
            SpherePoint(lsst.sphgeom.Vector3d(0.1, inf, 0.3)),
            SpherePoint(lsst.sphgeom.Vector3d(0.1, 0.2, -inf)),
        ]

        for finite in finiteList:
            self.assertIsInstance(finite.isFinite(), bool)
            self.assertTrue(finite.isFinite())

        for nonFinite in nonFiniteList:
            self.assertIsInstance(nonFinite.isFinite(), bool)
            self.assertFalse(nonFinite.isFinite())

    def testGetItemError(self):
        """Test if indexing correctly handles invalid input.
        """
        point = SpherePoint(lsst.sphgeom.Vector3d(1.0, 1.0, 1.0))

        with self.assertRaises(IndexError):
            point[2]
        with self.assertRaises(IndexError):
            point[-3]

    def testGetItemValue(self):
        """Test if indexing returns the expected value.
        """
        for point in self.pointSet:
            self.assertIsInstance(point[-2], geom.Angle)
            self.assertIsInstance(point[-1], geom.Angle)
            self.assertIsInstance(point[0], geom.Angle)
            self.assertIsInstance(point[1], geom.Angle)

            if not math.isnan(point.getLongitude().asRadians()):
                self.assertEqual(point.getLongitude(), point[-2])
                self.assertEqual(point.getLongitude(), point[0])
            else:
                self.assertTrue(math.isnan(point[-2].asRadians()))
                self.assertTrue(math.isnan(point[0].asRadians()))
            if not math.isnan(point.getLatitude().asRadians()):
                self.assertEqual(point.getLatitude(), point[-1])
                self.assertEqual(point.getLatitude(), point[1])
            else:
                self.assertTrue(math.isnan(point[-1].asRadians()))
                self.assertTrue(math.isnan(point[1].asRadians()))

    def testEquality(self):
        """Test if tests for equality treat SpherePoints as values.
        """
        # (In)equality is determined by value, not identity.
        # See DM-2347, DM-2465. These asserts are testing the
        # functionality of `==` and `!=` and should not be changed.
        for lon1, lat1 in self._dataset:
            point1 = SpherePoint(lon1, lat1)
            self.assertIsInstance(point1 == point1, bool)
            self.assertIsInstance(point1 != point1, bool)
            if point1.isFinite():
                self.assertTrue(point1 == point1)
                self.assertFalse(point1 != point1)

                pointCopy = copy.deepcopy(point1)
                self.assertIsNot(pointCopy, point1)
                self.assertEqual(pointCopy, point1)
                self.assertEqual(point1, pointCopy)
                self.assertFalse(pointCopy != point1)
                self.assertFalse(point1 != pointCopy)
            else:
                self.assertFalse(point1 == point1)
                self.assertTrue(point1 != point1)

            for lon2, lat2 in self._dataset:
                point2 = SpherePoint(lon2, lat2)
                if lon1 == lon2 and lat1 == lat2 and point1.isFinite() and point2.isFinite():
                    # note: the isFinite checks are needed because if longitude is infinite
                    # then the resulting SpherePoint has nan as its longitude, due to wrapping
                    self.assertFalse(point2 != point1)
                    self.assertFalse(point1 != point2)
                    self.assertTrue(point2 == point1)
                    self.assertTrue(point1 == point2)
                else:
                    self.assertTrue(point2 != point1)
                    self.assertTrue(point1 != point2)
                    self.assertFalse(point2 == point1)
                    self.assertFalse(point1 == point2)

        # Test for transitivity (may be assumed by algorithms).
        for delta in [10.0**(0.1*x) for x in range(-150, -49, 5)]:
            self.checkTransitive(delta*radians)

    def checkTransitive(self, delta):
        """Test if equality is transitive even for close points.

        This test prevents misuse of approximate floating-point
        equality -- if `__eq__` is implemented using AFP, then this
        test will fail for some value of `delta`. Testing multiple
        values is recommended.

        Parameters
        ----------
        delta : `number`
            The separation, in degrees, at which point equality may
            become intransitive.
        """
        for lon, lat in self._dataset:
            point1 = SpherePoint(lon - delta, lat)
            point2 = SpherePoint(lon, lat)
            point3 = SpherePoint(lon + delta, lat)

            self.assertTrue(point1 != point2 or point2 !=
                            point3 or point1 == point3)
            self.assertTrue(point3 != point1 or point1 !=
                            point2 or point3 == point2)
            self.assertTrue(point2 == point3 or point3 !=
                            point1 or point2 == point1)

    def testBearingToValueOnEquator(self):
        """Test if bearingTo() returns the expected value from a point on the equator
        """
        lon0 = 90.0
        lat0 = 0.0   # These tests only work from the equator.
        arcLen = 10.0

        trials = [
            # Along celestial equator
            dict(lon=lon0, lat=lat0, bearing=0.0,
                 lonEnd=lon0+arcLen, latEnd=lat0),
            # Along a meridian
            dict(lon=lon0, lat=lat0, bearing=90.0,
                 lonEnd=lon0, latEnd=lat0+arcLen),
            # 180 degree arc (should go to antipodal point)
            dict(lon=lon0, lat=lat0, bearing=45.0,
                 lonEnd=lon0+180.0, latEnd=-lat0),
            #
            dict(lon=lon0, lat=lat0, bearing=45.0,
                 lonEnd=lon0+90.0, latEnd=lat0 + 45.0),
            dict(lon=lon0, lat=lat0, bearing=225.0,
                 lonEnd=lon0-90.0, latEnd=lat0 - 45.0),
            dict(lon=lon0, lat=np.nextafter(-90.0, inf),
                 bearing=90.0, lonEnd=lon0, latEnd=0.0),
            dict(lon=lon0, lat=np.nextafter(-90.0, inf),
                 bearing=0.0, lonEnd=lon0 + 90.0, latEnd=0.0),
            # Argument at a pole should work
            dict(lon=lon0, lat=lat0, bearing=270.0, lonEnd=lon0, latEnd=-90.0),
            # Support for non-finite values
            dict(lon=lon0, lat=nan, bearing=nan, lonEnd=lon0, latEnd=45.0),
            dict(lon=lon0, lat=lat0, bearing=nan, lonEnd=nan, latEnd=90.0),
            dict(lon=inf, lat=lat0, bearing=nan, lonEnd=lon0, latEnd=42.0),
            dict(lon=lon0, lat=lat0, bearing=nan, lonEnd=-inf, latEnd=42.0),
        ]

        for trial in trials:
            origin = SpherePoint(trial['lon']*degrees, trial['lat']*degrees)
            end = SpherePoint(trial['lonEnd']*degrees, trial['latEnd']*degrees)
            bearing = origin.bearingTo(end)

            self.assertIsInstance(bearing, geom.Angle)
            if origin.isFinite() and end.isFinite():
                self.assertGreaterEqual(bearing.asDegrees(), 0.0)
                self.assertLess(bearing.asDegrees(), 360.0)
            if origin.separation(end).asDegrees() != 180.0:
                if not math.isnan(trial['bearing']):
                    self.assertAlmostEqual(
                        trial['bearing'], bearing.asDegrees(), 12)
                else:
                    self.assertTrue(math.isnan(bearing.asRadians()))

    def testBearingToValueSameLongitude(self):
        """Test that bearingTo() returns +/- 90 for two points on the same longitude
        """
        for longDeg in (0, 55, 270):
            for lat0Deg in (-90, -5, 0, 44, 90):
                sp0 = SpherePoint(longDeg, lat0Deg, degrees)
                for lat1Deg in (-90, -41, 1, 41, 90):
                    if lat0Deg == lat1Deg:
                        continue
                    sp1 = SpherePoint(longDeg, lat1Deg, degrees)
                    if sp0.atPole() and sp1.atPole():
                        # the points are at opposite poles; any bearing may be returned
                        continue
                    bearing = sp0.bearingTo(sp1)
                    if lat1Deg > lat0Deg:
                        self.assertAnglesAlmostEqual(bearing, 90 * degrees)
                    else:
                        self.assertAnglesAlmostEqual(bearing, -90 * degrees)

    def testBearingToFromPole(self):
        """Test if bearingTo() returns the expected value from a point at a pole
        """
        for long0Deg in (0, 55, 270):
            for atSouthPole in (False, True):
                lat0Deg = -90 if atSouthPole else 90
                sp0 = SpherePoint(long0Deg, lat0Deg, degrees)
                for long1Deg in (0, 55, 270):
                    for lat1Deg in (-89, 0, 89):
                        sp1 = SpherePoint(long1Deg, lat1Deg, degrees)
                        desiredBearing = ((long1Deg - long0Deg) - 90) * degrees
                        if atSouthPole:
                            desiredBearing *= -1
                        measuredBearing = sp0.bearingTo(sp1)
                        self.assertAnglesAlmostEqual(desiredBearing, measuredBearing)

    def testBearingToValueSingular(self):
        """White-box test: bearingTo() may be unstable if points are near opposite poles.

        This test is motivated by an error analysis of the `bearingTo`
        implementation. It may become irrelevant if the implementation
        changes.
        """
        southPole = SpherePoint(0.0*degrees, self.nextUp(-90.0*degrees))
        northPoleSame = SpherePoint(0.0*degrees, self.nextDown(90.0*degrees))
        # Don't let it be on exactly the opposite side.
        northPoleOpposite = SpherePoint(
            180.0*degrees, self.nextDown(northPoleSame.getLatitude()))

        self.assertAnglesAlmostEqual(southPole.bearingTo(northPoleSame),
                                     geom.HALFPI*geom.radians)
        self.assertAnglesAlmostEqual(southPole.bearingTo(northPoleOpposite),
                                     (geom.PI + geom.HALFPI)*geom.radians)

    def testSeparationValueGeneric(self):
        """Test if separation() returns the correct value.
        """
        # This should cover arcs over the meridian, across the pole, etc.
        # Do not use sphgeom as an oracle, in case SpherePoint uses it
        # internally.
        for lon1, lat1 in self._dataset:
            point1 = SpherePoint(lon1, lat1)
            x1, y1, z1 = SpherePointTestSuite.toVector(lon1, lat1)
            for lon2, lat2 in self._dataset:
                point2 = SpherePoint(lon2, lat2)
                if lon1 != lon2 or lat1 != lat2:
                    # Numerically unstable at small angles, but that's ok.
                    x2, y2, z2 = SpherePointTestSuite.toVector(lon2, lat2)
                    expected = math.acos(x1*x2 + y1*y2 + z1*z2)
                else:
                    expected = 0.0

                sep = point1.separation(point2)
                self.assertIsInstance(sep, geom.Angle)
                if point1.isFinite() and point2.isFinite():
                    self.assertGreaterEqual(sep.asDegrees(), 0.0)
                    self.assertLessEqual(sep.asDegrees(), 180.0)
                    self.assertAlmostEqual(expected, sep.asRadians())
                    self.assertAnglesAlmostEqual(
                        sep, point2.separation(point1))
                else:
                    self.assertTrue(math.isnan(sep.asRadians()))
                    self.assertTrue(math.isnan(
                        point2.separation(point1).asRadians()))

    def testSeparationValueAbsolute(self):
        """Test if separation() returns specific values.
        """
        # Test from "Meeus, p. 110" (test originally written for coord::Coord;
        # don't know exact reference)
        spica = SpherePoint(201.2983, -11.1614, degrees)
        arcturus = SpherePoint(213.9154, 19.1825, degrees)

        # Verify to precision of quoted distance and positions.
        self.assertAlmostEqual(
            32.7930, spica.separation(arcturus).asDegrees(), 4)

        # Verify small angles: along a constant ra, add an arcsec to spica dec.
        epsilon = 1.0*geom.arcseconds
        spicaPlus = SpherePoint(spica.getLongitude(),
                                spica.getLatitude() + epsilon)

        self.assertAnglesAlmostEqual(epsilon, spicaPlus.separation(spica))

    def testSeparationPoles(self):
        """White-box test: all representations of a pole should have the same distance to another point.
        """
        southPole1 = SpherePoint(-30.0, -90.0, degrees)
        southPole2 = SpherePoint(183.0, -90.0, degrees)
        regularPoint = SpherePoint(42.0, 45.0, degrees)
        expectedSep = (45.0 + 90.0)*degrees

        self.assertAnglesAlmostEqual(
            expectedSep, southPole1.separation(regularPoint))
        self.assertAnglesAlmostEqual(
            expectedSep, regularPoint.separation(southPole1))
        self.assertAnglesAlmostEqual(
            expectedSep, southPole2.separation(regularPoint))
        self.assertAnglesAlmostEqual(
            expectedSep, regularPoint.separation(southPole2))

    @staticmethod
    def toVector(longitude, latitude):
        """Converts a set of spherical coordinates to a 3-vector.

        The conversion shall not be performed by any library, to ensure
        that the test case does not duplicate the code being tested.

        Parameters
        ----------
        longitude : `Angle`
            The longitude (right ascension, azimuth, etc.) of the
            position.
        latitude : `Angle`
            The latitude (declination, elevation, etc.) of the
            position.

        Returns
        -------
        x, y, z : `number`
            Components of the unit vector representation of
            `(longitude, latitude)`
        """
        alpha = longitude.asRadians()
        delta = latitude.asRadians()
        if math.isnan(alpha) or math.isinf(alpha) or math.isnan(delta) or math.isinf(delta):
            return (nan, nan, nan)

        x = math.cos(alpha)*math.cos(delta)
        y = math.sin(alpha)*math.cos(delta)
        z = math.sin(delta)
        return (x, y, z)

    def testRotatedValue(self):
        """Test if rotated() returns the expected value.
        """
        # Try rotating about the equatorial pole (ie. along a parallel).
        longitude = 90.0
        latitudes = [0.0, 30.0, 60.0]
        arcLen = 10.0
        pole = SpherePoint(0.0*degrees, 90.0*degrees)
        for latitude in latitudes:
            point = SpherePoint(longitude*degrees, latitude*degrees)
            newPoint = point.rotated(pole, arcLen*degrees)

            self.assertIsInstance(newPoint, SpherePoint)
            self.assertAlmostEqual(
                longitude + arcLen, newPoint.getLongitude().asDegrees())
            self.assertAlmostEqual(
                latitude, newPoint.getLatitude().asDegrees())

        # Try with pole = vernal equinox and rotate up the 90 degree meridian.
        pole = SpherePoint(0.0*degrees, 0.0*degrees)
        for latitude in latitudes:
            point = SpherePoint(longitude*degrees, latitude*degrees)
            newPoint = point.rotated(pole, arcLen*degrees)

            self.assertAlmostEqual(
                longitude, newPoint.getLongitude().asDegrees())
            self.assertAlmostEqual(
                latitude + arcLen, newPoint.getLatitude().asDegrees())

        # Test accuracy close to coordinate pole
        point = SpherePoint(90.0*degrees, np.nextafter(90.0, -inf)*degrees)
        newPoint = point.rotated(pole, 90.0*degrees)
        self.assertAlmostEqual(270.0, newPoint.getLongitude().asDegrees())
        self.assertAlmostEqual(90.0 - np.nextafter(90.0, -inf),
                               newPoint.getLatitude().asDegrees())

        # Generic pole; can't predict position, but test for rotation
        # invariant.
        pole = SpherePoint(283.5*degrees, -23.6*degrees)
        for lon, lat in self._dataset:
            point = SpherePoint(lon, lat)
            dist = point.separation(pole)
            newPoint = point.rotated(pole, -32.4*geom.radians)

            self.assertNotAlmostEqual(point.getLongitude().asDegrees(),
                                      newPoint.getLongitude().asDegrees())
            self.assertNotAlmostEqual(point.getLatitude().asDegrees(),
                                      newPoint.getLatitude().asDegrees())
            self.assertAnglesAlmostEqual(dist, newPoint.separation(pole))

        # Non-finite values give undefined rotations
        for latitude in latitudes:
            point = SpherePoint(longitude*degrees, latitude*degrees)
            nanPoint = point.rotated(pole, nan*degrees)
            infPoint = point.rotated(pole, inf*degrees)

            self.assertTrue(math.isnan(nanPoint.getLongitude().asRadians()))
            self.assertTrue(math.isnan(nanPoint.getLatitude().asRadians()))
            self.assertTrue(math.isnan(infPoint.getLongitude().asRadians()))
            self.assertTrue(math.isnan(infPoint.getLatitude().asRadians()))

        # Non-finite points rotate into non-finite points
        for point in [
            SpherePoint(-inf*degrees, 1.0*radians),
            SpherePoint(32.0*degrees, nan*radians),
        ]:
            newPoint = point.rotated(pole, arcLen*degrees)
            self.assertTrue(math.isnan(nanPoint.getLongitude().asRadians()))
            self.assertTrue(math.isnan(nanPoint.getLatitude().asRadians()))
            self.assertTrue(math.isnan(infPoint.getLongitude().asRadians()))
            self.assertTrue(math.isnan(infPoint.getLatitude().asRadians()))

        # Rotation around non-finite poles undefined
        for latitude in latitudes:
            point = SpherePoint(longitude*degrees, latitude*degrees)
            for pole in [
                SpherePoint(-inf*degrees, 1.0*radians),
                SpherePoint(32.0*degrees, nan*radians),
            ]:
                newPoint = point.rotated(pole, arcLen*degrees)
                self.assertTrue(math.isnan(
                    nanPoint.getLongitude().asRadians()))
                self.assertTrue(math.isnan(nanPoint.getLatitude().asRadians()))
                self.assertTrue(math.isnan(
                    infPoint.getLongitude().asRadians()))
                self.assertTrue(math.isnan(infPoint.getLatitude().asRadians()))

    def testRotatedAlias(self):
        """White-box test: all representations of a pole should rotate into the same point.
        """
        longitudes = [0.0, 90.0, 242.0]
        latitude = 90.0
        arcLen = 10.0
        pole = SpherePoint(90.0*degrees, 0.0*degrees)
        for longitude in longitudes:
            point = SpherePoint(longitude*degrees, latitude*degrees)
            newPoint = point.rotated(pole, arcLen*degrees)

            self.assertAlmostEqual(0.0, newPoint.getLongitude().asDegrees())
            self.assertAlmostEqual(80.0, newPoint.getLatitude().asDegrees())

    def testOffsetValue(self):
        """Test if offset() returns the expected value.
        """
        # This should cover arcs over the meridian, across the pole, etc.
        for lon1, lat1 in self._dataset:
            point1 = SpherePoint(lon1, lat1)
            for lon2, lat2 in self._dataset:
                if lon1 == lon2 and lat1 == lat2:
                    continue
                point2 = SpherePoint(lon2, lat2)
                bearing = point1.bearingTo(point2)
                distance = point1.separation(point2)

                # offsetting point1 by bearing and distance should produce the same result as point2
                newPoint = point1.offset(bearing, distance)
                self.assertIsInstance(newPoint, SpherePoint)
                self.assertSpherePointsAlmostEqual(point2, newPoint)
                if newPoint.atPole():
                    self.assertAnglesAlmostEqual(newPoint.getLongitude(), 0*degrees)

                # measuring the separation and bearing from point1 to the new point
                # should produce the requested separation and bearing
                measuredDistance = point1.separation(newPoint)
                self.assertAnglesAlmostEqual(measuredDistance, distance)
                if abs(measuredDistance.asDegrees() - 180) > 1e-5:
                    # The two points are not opposite each other on the sphere,
                    # so the bearing has a well defined value
                    measuredBearing = point1.bearingTo(newPoint)
                    self.assertAnglesAlmostEqual(measuredBearing, bearing)

                # offset by a negative amount in the opposite direction should produce the same result
                newPoint2 = point1.offset(bearing + 180 * degrees, -distance)
                self.assertIsInstance(newPoint2, SpherePoint)
                # check angular separation (longitude is checked below)
                self.assertSpherePointsAlmostEqual(newPoint, newPoint2)

                if point1.isFinite() and point2.isFinite():
                    if not point2.atPole():
                        self.assertAnglesAlmostEqual(
                            point2.getLongitude(), newPoint.getLongitude())
                        self.assertAnglesAlmostEqual(
                            point2.getLongitude(), newPoint2.getLongitude())
                    self.assertAnglesAlmostEqual(
                        point2.getLatitude(), newPoint.getLatitude())
                    self.assertAnglesAlmostEqual(
                        point2.getLatitude(), newPoint2.getLatitude())
                else:
                    self.assertTrue(math.isnan(
                        newPoint.getLongitude().asRadians()))
                    self.assertTrue(math.isnan(
                        newPoint2.getLongitude().asRadians()))
                    self.assertTrue(math.isnan(
                        newPoint.getLatitude().asRadians()))
                    self.assertTrue(math.isnan(
                        newPoint2.getLatitude().asRadians()))

        # Test precision near the poles
        lon = 123.0*degrees
        almostPole = SpherePoint(lon, self.nextDown(90.0*degrees))
        goSouth = almostPole.offset(-90.0*degrees, 90.0*degrees)
        self.assertAnglesAlmostEqual(lon, goSouth.getLongitude())
        self.assertAnglesAlmostEqual(0.0*degrees, goSouth.getLatitude())
        goEast = almostPole.offset(0.0*degrees, 90.0*degrees)
        self.assertAnglesAlmostEqual(lon + 90.0*degrees, goEast.getLongitude())
        self.assertAnglesAlmostEqual(0.0*degrees, goEast.getLatitude())

    def testOffsetTangentPlane(self):
        """Test offsets on a tangent plane (good for small angles)"""

        c0 = SpherePoint(0.0, 0.0, geom.degrees)

        for dRaDeg in (0.0123, 0.0, -0.0321):
            dRa = dRaDeg*geom.degrees
            for dDecDeg in (0.0543, 0.0, -0.0987):
                dDec = dDecDeg*geom.degrees
                c1 = SpherePoint(dRa, dDec)

                offset = c0.getTangentPlaneOffset(c1)

                # This more-or-less works for small angles because c0 is 0,0
                expectedOffset = [
                    math.tan(dRa.asRadians())*geom.radians,
                    math.tan(dDec.asRadians())*geom.radians,
                ]

                for i in range(2):
                    self.assertAnglesAlmostEqual(offset[i], expectedOffset[i])

    def testIterResult(self):
        """Test if iteration returns the expected values.
        """
        for point in self.pointSet:
            if not point.isFinite():
                continue

            # Test mechanics directly
            it = iter(point)
            self.assertEqual(point.getLongitude(), next(it))
            self.assertEqual(point.getLatitude(), next(it))
            with self.assertRaises(StopIteration):
                next(it)

            # Intended use case
            lon, lat = point
            self.assertEqual(point.getLongitude(), lon)
            self.assertEqual(point.getLatitude(), lat)

    def testStrValue(self):
        """Test if __str__ produces output consistent with its spec.

        This is necessarily a loose test, as the behavior of __str__
        is (deliberately) incompletely specified.
        """
        for point in self.pointSet:
            numbers = re.findall(r'(?:\+|-)?(?:[\d.]+|nan|inf)', str(point))
            self.assertEqual(2, len(numbers),
                             "String '%s' should have exactly two coordinates." % (point,))

            # Low precision to allow for only a few digits in string.
            if not math.isnan(point.getLongitude().asRadians()):
                self.assertAlmostEqual(
                    point.getLongitude().asDegrees(), float(numbers[0]), delta=1e-6)
            else:
                self.assertRegex(numbers[0], r'-?nan')
            if not math.isnan(point.getLatitude().asRadians()):
                self.assertAlmostEqual(
                    point.getLatitude().asDegrees(), float(numbers[1]), delta=1e-6)
                # Latitude must be signed
                self.assertTrue(numbers[1].startswith("+") or
                                numbers[1].startswith("-"))
            else:
                # Some C++ compilers will output NaN with a sign, others won't
                self.assertRegex(numbers[1], r'(?:\+|-)?nan')

    def testReprValue(self):
        """Test if __repr__ is a machine-readable representation.
        """
        for point in self.pointSet:
            pointRepr = repr(point)
            self.assertIn("degrees", pointRepr)
            self.assertEqual(2, len(pointRepr.split(",")))

            spcopy = eval(pointRepr)
            self.assertAnglesAlmostEqual(
                point.getLongitude(), spcopy.getLongitude())
            self.assertAnglesAlmostEqual(
                point.getLatitude(), spcopy.getLatitude())

    def testAverageSpherePoint(self):
        """Test the averageSpherePoint function"""

        def checkCircle(center, start, numPts, maxSep=1.0e-9*geom.arcseconds):
            """Generate points in a circle; test that average is in the center
            """
            coords = []
            deltaAngle = 360*degrees / numPts
            for ii in range(numPts):
                new = start.rotated(center, ii*deltaAngle)
                coords.append(new)
            result = geom.averageSpherePoint(coords)
            self.assertSpherePointsAlmostEqual(center, result, maxSep=maxSep)

        for numPts in (2, 3, 120):
            for center, start in (
                    # RA=0=360 border
                    (SpherePoint(0, 0, geom.degrees), SpherePoint(5, 0, geom.degrees)),
                    # North pole
                    (SpherePoint(0, 90, geom.degrees), SpherePoint(0, 85, geom.degrees)),
                    # South pole
                    (SpherePoint(0, -90, geom.degrees), SpherePoint(0, -85, geom.degrees)),
            ):
                checkCircle(center=center, start=start, numPts=numPts)

    def nextUp(self, angle):
        """Returns the smallest angle that is larger than the argument.
        """
        return np.nextafter(angle.asRadians(), inf)*radians

    def nextDown(self, angle):
        """Returns the largest angle that is smaller than the argument.
        """
        return np.nextafter(angle.asRadians(), -inf)*radians


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
