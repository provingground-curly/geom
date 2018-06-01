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

"""Utilities that should be imported into the lsst.geom namespace when lsst.geom is used

In the case of the assert functions, importing them makes them available in lsst.utils.tests.TestCase
"""
__all__ = []

import math

import numpy as np

import lsst.utils.tests
from .angle import arcseconds


def extraMsg(msg):
    """Format extra error message, if any
    """
    if msg:
        return ": " + msg
    return ""


@lsst.utils.tests.inTestCase
def assertAnglesAlmostEqual(testCase, ang0, ang1, maxDiff=0.001*arcseconds,
                            ignoreWrap=True, msg="Angles differ"):
    r"""Assert that two `~lsst.afw.geom.Angle`\ s are almost equal, ignoring wrap differences by default

    Parameters
    ----------
    testCase : `unittest.TestCase`
        test case the test is part of; an object supporting one method: fail(self, msgStr)
    ang0 : `lsst.afw.geom.Angle`
        angle 0
    ang1 : `an lsst.afw.geom.Angle`
        angle 1
    maxDiff : `an lsst.afw.geom.Angle`
        maximum difference between the two angles
    ignoreWrap : `bool`
        ignore wrap when comparing the angles?
        - if True then wrap is ignored, e.g. 0 and 360 degrees are considered equal
        - if False then wrap matters, e.g. 0 and 360 degrees are considered different
    msg : `str`
        exception message prefix; details of the error are appended after ": "

    Raises
    ------
    AssertionError
        Raised if the difference is greater than ``maxDiff``
    """
    measDiff = ang1 - ang0
    if ignoreWrap:
        measDiff = measDiff.wrapCtr()
    if abs(measDiff) > maxDiff:
        testCase.fail("%s: measured difference %s arcsec > max allowed %s arcsec" %
                      (msg, measDiff.asArcseconds(), maxDiff.asArcseconds()))


@lsst.utils.tests.inTestCase
def assertPairsAlmostEqual(testCase, pair0, pair1, maxDiff=1e-7, msg="Pairs differ"):
    """Assert that two Cartesian points are almost equal.

    Each point can be any indexable pair of two floats, including
    Point2D or Extent2D, a list or a tuple.

    Parameters
    ----------
    testCase : `unittest.TestCase`
        test case the test is part of; an object supporting one method: fail(self, msgStr)
    pair0 : pair of `float`
        pair 0
    pair1 : pair of `floats`
        pair 1
    maxDiff : `float`
        maximum radial separation between the two points
    msg : `str`
        exception message prefix; details of the error are appended after ": "

    Raises
    ------
    AssertionError
        Raised if the radial difference is greater than ``maxDiff``

    Notes
    -----
    .. warning::

       Does not compare types, just compares values.
    """
    if len(pair0) != 2:
        raise RuntimeError("len(pair0)=%s != 2" % (len(pair0),))
    if len(pair1) != 2:
        raise RuntimeError("len(pair1)=%s != 2" % (len(pair1),))

    pairDiff = [float(pair1[i] - pair0[i]) for i in range(2)]
    measDiff = math.hypot(*pairDiff)
    if measDiff > maxDiff:
        testCase.fail("%s: measured radial distance = %s > maxDiff = %s, pair0=(%r, %r), pair1=(%r, %r)" %
                      (msg, measDiff, maxDiff, pair0[0], pair0[1], pair1[0], pair1[1]))


@lsst.utils.tests.inTestCase
def assertPairListsAlmostEqual(testCase, list0, list1, maxDiff=1e-7, msg=None):
    """Assert that two lists of Cartesian points are almost equal

    Each point can be any indexable pair of two floats, including
    Point2D or Extent2D, a list or a tuple.

    Parameters
    ----------
    testCase : `unittest.TestCase`
        test case the test is part of; an object supporting one method: fail(self, msgStr)
    list0 : `list` of pairs of `float`
        list of pairs 0
    list1 : `list` of pairs of `float`
        list of pairs 1
    maxDiff : `float`
        maximum radial separation between the two points
    msg : `str`
        additional information for the error message; appended after ": "

    Raises
    ------
    AssertionError
        Raised if the radial difference is greater than ``maxDiff``

    Notes
    -----
    .. warning::

       Does not compare types, just values.
    """
    testCase.assertEqual(len(list0), len(list1))
    lenList1 = np.array([len(val) for val in list0])
    lenList2 = np.array([len(val) for val in list1])
    testCase.assertTrue(np.all(lenList1 == 2))
    testCase.assertTrue(np.all(lenList2 == 2))

    diffArr = np.array([(val0[0] - val1[0], val0[1] - val1[1])
                        for val0, val1 in zip(list0, list1)], dtype=float)
    sepArr = np.hypot(diffArr[:, 0], diffArr[:, 1])
    badArr = sepArr > maxDiff
    if np.any(badArr):
        maxInd = np.argmax(sepArr)
        testCase.fail("PairLists differ in %s places; max separation is at %s: %s > %s%s" %
                      (np.sum(badArr), maxInd, sepArr[maxInd], maxDiff, extraMsg(msg)))


@lsst.utils.tests.inTestCase
def assertSpherePointsAlmostEqual(testCase, sp0, sp1, maxSep=0.001*arcseconds, msg=""):
    r"""Assert that two `~lsst.afw.geom.SpherePoint`\ s are almost equal

    Parameters
    ----------
    testCase : `unittest.TestCase`
        test case the test is part of; an object supporting one method: fail(self, msgStr)
    sp0 : `lsst.afw.geom.SpherePoint`
        SpherePoint 0
    sp1 : `lsst.afw.geom.SpherePoint`
        SpherePoint 1
    maxSep : `lsst.afw.geom.Angle`
        maximum separation
    msg : `str`
        extra information to be printed with any error message
    """
    if sp0.separation(sp1) > maxSep:
        testCase.fail("Angular separation between %s and %s = %s\" > maxSep = %s\"%s" %
                      (sp0, sp1, sp0.separation(sp1).asArcseconds(), maxSep.asArcseconds(), extraMsg(msg)))


@lsst.utils.tests.inTestCase
def assertSpherePointListsAlmostEqual(testCase, splist0, splist1, maxSep=0.001*arcseconds, msg=None):
    r"""Assert that two lists of `~lsst.afw.geom.SpherePoint`\ s are almost equal

    Parameters
    ----------
    testCase : `unittest.TestCase`
        test case the test is part of; an object supporting one method: fail(self, msgStr)
    splist0 : `list` of `lsst.afw.geom.SpherePoint`
        list of SpherePoints 0
    splist1 : `list` of `lsst.afw.geom.SpherePoint`
        list of SpherePoints 1
    maxSep : `lsst.afw.geom.Angle`
        maximum separation
    msg : `str`
        exception message prefix; details of the error are appended after ": "
    """
    testCase.assertEqual(len(splist0), len(splist1), msg=msg)
    sepArr = np.array([sp0.separation(sp1)
                       for sp0, sp1 in zip(splist0, splist1)])
    badArr = sepArr > maxSep
    if np.any(badArr):
        maxInd = np.argmax(sepArr)
        testCase.fail("SpherePointLists differ in %s places; max separation is at %s: %s\" > %s\"%s" %
                      (np.sum(badArr), maxInd, sepArr[maxInd].asArcseconds(),
                       maxSep.asArcseconds(), extraMsg(msg)))


@lsst.utils.tests.inTestCase
def assertBoxesAlmostEqual(testCase, box0, box1, maxDiff=1e-7, msg="Boxes differ"):
    """Assert that two boxes (`~lsst.afw.geom.Box2D` or `~lsst.afw.geom.Box2I`) are almost equal

    Parameters
    ----------
    testCase : `unittest.TestCase`
        test case the test is part of; an object supporting one method: fail(self, msgStr)
    box0 : `lsst.afw.geom.Box2D` or `lsst.afw.geom.Box2I`
        box 0
    box1 : `lsst.afw.geom.Box2D` or `lsst.afw.geom.Box2I`
        box 1
    maxDiff : `float`
        maximum radial separation between the min points and max points
    msg : `str`
        exception message prefix; details of the error are appended after ": "

    Raises
    ------
    AssertionError
        Raised if the radial difference of the min points or max points is greater than maxDiff

    Notes
    -----
    .. warning::

       Does not compare types, just compares values.
    """
    assertPairsAlmostEqual(testCase, box0.getMin(),
                           box1.getMin(), maxDiff=maxDiff, msg=msg + ": min")
    assertPairsAlmostEqual(testCase, box0.getMax(),
                           box1.getMax(), maxDiff=maxDiff, msg=msg + ": max")
