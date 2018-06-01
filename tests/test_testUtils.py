#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import math
import unittest

# importing lsst.geom adds assert methods to lsst.utils.tests.TestCase,
# and these are most or all of what this module tests
import lsst.geom
import lsst.utils.tests


class TestTestUtils(lsst.utils.tests.TestCase):
    """Test test methods added to lsst.utils.tests.TestCase
    """
    def testAssertAnglesAlmostEqual(self):
        """Test assertAnglesAlmostEqual"""
        for angDeg in (0, 45, -75):
            ang0 = angDeg*lsst.geom.degrees
            self.assertAnglesAlmostEqual(
                ang0,
                ang0 + 0.01*lsst.geom.arcseconds,
                maxDiff=0.010001*lsst.geom.arcseconds,
            )
            with self.assertRaises(AssertionError):
                self.assertAnglesAlmostEqual(
                    ang0,
                    ang0 + 0.01*lsst.geom.arcseconds,
                    maxDiff=0.009999*lsst.geom.arcseconds,
                )

            self.assertAnglesAlmostEqual(
                ang0,
                ang0 - 0.01*lsst.geom.arcseconds,
                maxDiff=0.010001*lsst.geom.arcseconds,
            )
            with self.assertRaises(AssertionError):
                self.assertAnglesAlmostEqual(
                    ang0,
                    ang0 - 0.01*lsst.geom.arcseconds,
                    maxDiff=0.009999*lsst.geom.arcseconds,
                )

            self.assertAnglesAlmostEqual(
                ang0 - 720*lsst.geom.degrees,
                ang0 + 0.01*lsst.geom.arcseconds,
                maxDiff=0.010001*lsst.geom.arcseconds,
            )
            with self.assertRaises(AssertionError):
                self.assertAnglesAlmostEqual(
                    ang0 - 720*lsst.geom.degrees,
                    ang0 + 0.01*lsst.geom.arcseconds,
                    ignoreWrap=False,
                    maxDiff=0.010001*lsst.geom.arcseconds,
                )
            with self.assertRaises(AssertionError):
                self.assertAnglesAlmostEqual(
                    ang0 - 720*lsst.geom.degrees,
                    ang0 + 0.01*lsst.geom.arcseconds,
                    maxDiff=0.009999*lsst.geom.arcseconds,
                )

            self.assertAnglesAlmostEqual(
                ang0,
                ang0 + 360*lsst.geom.degrees + 0.01*lsst.geom.arcseconds,
                maxDiff=0.010001*lsst.geom.arcseconds,
            )
            with self.assertRaises(AssertionError):
                self.assertAnglesAlmostEqual(
                    ang0,
                    ang0 + 360*lsst.geom.degrees + 0.01*lsst.geom.arcseconds,
                    ignoreWrap=False,
                    maxDiff=0.010001*lsst.geom.arcseconds,
                )
            with self.assertRaises(AssertionError):
                self.assertAnglesAlmostEqual(
                    ang0,
                    ang0 + 360*lsst.geom.degrees + 0.01*lsst.geom.arcseconds,
                    maxDiff=0.009999*lsst.geom.arcseconds,
                )

    def testAssertBoxesAlmostEqual(self):
        """Test assertBoxesAlmostEqual"""
        for min0 in ((0, 0), (-1000.5, 5000.1)):
            min0 = lsst.geom.Point2D(*min0)
            for extent0 in ((2.01, 3.01), (5432, 2342)):
                extent0 = lsst.geom.Extent2D(*extent0)
                box0 = lsst.geom.Box2D(min0, extent0)
                self.assertBoxesAlmostEqual(box0, box0, maxDiff=1e-7)
                for deltaExtent in ((0.001, -0.001), (2, -3)):
                    deltaExtent = lsst.geom.Extent2D(*deltaExtent)
                    box1 = lsst.geom.Box2D(
                        box0.getMin() + deltaExtent, box0.getMax())
                    radDiff = math.hypot(*deltaExtent)
                    self.assertBoxesAlmostEqual(
                        box0, box1, maxDiff=radDiff*1.00001)
                    with self.assertRaises(AssertionError):
                        self.assertBoxesAlmostEqual(
                            box0, box1, maxDiff=radDiff*0.99999)

                    box2 = lsst.geom.Box2D(
                        box0.getMin() - deltaExtent, box0.getMax())
                    self.assertBoxesAlmostEqual(
                        box0, box2, maxDiff=radDiff*1.00001)
                    with self.assertRaises(AssertionError):
                        self.assertBoxesAlmostEqual(
                            box0, box2, maxDiff=radDiff*0.999999)

                    box3 = lsst.geom.Box2D(
                        box0.getMin(), box0.getMax() + deltaExtent)
                    self.assertBoxesAlmostEqual(
                        box0, box3, maxDiff=radDiff*1.00001)
                    with self.assertRaises(AssertionError):
                        self.assertBoxesAlmostEqual(
                            box0, box3, maxDiff=radDiff*0.999999)

    def testAssertSpherePointsAlmostEqual(self):
        """Test assertSpherePointsAlmostEqual"""
        for raDecDeg in ((45, 45), (-70, 89), (130, -89.5)):
            raDecDeg = [val*lsst.geom.degrees for val in raDecDeg]
            sp0 = lsst.geom.SpherePoint(*raDecDeg)
            self.assertSpherePointsAlmostEqual(
                sp0, sp0, maxSep=1e-7*lsst.geom.arcseconds)
            # make sure specifying msg is acceptable
            self.assertSpherePointsAlmostEqual(
                sp0, sp0, maxSep=1e-7*lsst.geom.arcseconds, msg="any")

            for offAng in (0, 45, 90):
                offAng = offAng*lsst.geom.degrees
                for offDist in (0.001, 0.1):
                    offDist = offDist*lsst.geom.arcseconds
                    sp1 = sp0.offset(bearing=offAng, amount=offDist)
                    self.assertSpherePointsAlmostEqual(
                        sp0, sp1, maxSep=offDist*1.00001)
                    with self.assertRaises(AssertionError):
                        self.assertSpherePointsAlmostEqual(
                            sp0, sp1, maxSep=offDist*0.99999)

                    # make sure msg is appended
                    try:
                        self.assertSpherePointsAlmostEqual(
                            sp0, sp1, maxSep=offDist*0.99999, msg="boo")
                        self.fail("Sphere point lists should be unequal")
                    except AssertionError as e:
                        errMsg = e.args[0]
                    self.assertTrue(errMsg.endswith("boo"))

            # test wraparound in RA
            sp2 = lsst.geom.SpherePoint(
                raDecDeg[0] + 360*lsst.geom.degrees, raDecDeg[1])
            self.assertSpherePointsAlmostEqual(
                sp0, sp2, maxSep=1e-7*lsst.geom.arcseconds)

    def testAssertSpherePointListsAlmostEqual(self):
        """Test assertSpherePointListsAlmostEqual
        """
        splist0 = [lsst.geom.SpherePoint(val[0]*lsst.geom.degrees, val[1]*lsst.geom.degrees)
                   for val in ((45, 45), (-70, 89), (130, -89.5))]
        self.assertSpherePointListsAlmostEqual(splist0, splist0)

        offDist = 1.1 * lsst.geom.arcseconds
        splist1 = [sp0.offset(bearing=bearDeg*lsst.geom.degrees, amount=offDist)
                   for sp0, bearDeg in zip(splist0, (-10, 78, 123))]
        self.assertSpherePointListsAlmostEqual(
            splist0, splist1, maxSep=offDist*1.00001)
        with self.assertRaises(AssertionError):
            self.assertSpherePointListsAlmostEqual(
                splist0, splist1, maxSep=offDist*0.99999)

        # make sure msg is appended
        try:
            self.assertSpherePointListsAlmostEqual(
                splist0, splist1, maxSep=offDist*0.99999, msg="boo")
            self.fail("Sphere point lists should be unequal")
        except AssertionError as e:
            errMsg = e.args[0]
        self.assertTrue(errMsg.endswith("boo"))

    def testAssertPairsAlmostEqual(self):
        """Test assertPairsAlmostEqual"""
        for pair0 in ((-5, 4), (-5, 0.001), (0, 0), (49, 0.1)):
            self.assertPairsAlmostEqual(pair0, pair0, maxDiff=1e-7)
            self.assertPairsAlmostEqual(lsst.geom.Point2D(*pair0),
                                        lsst.geom.Extent2D(*pair0), maxDiff=1e-7)
            for diff in ((0.001, 0), (-0.01, 0.03)):
                pair1 = [pair0[i] + diff[i] for i in range(2)]
                radialDiff = math.hypot(*diff)
                self.assertPairsAlmostEqual(
                    pair0, pair1, maxDiff=radialDiff+1e-7)
                with self.assertRaises(AssertionError):
                    self.assertPairsAlmostEqual(
                        pair0, pair1, maxDiff=radialDiff-1e-7)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
