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

import unittest
import itertools

import numpy as np

import lsst.utils.tests
import lsst.pex.exceptions
from lsst.geom import IntervalI, IntervalD


class IntervalTestData:

    def __init__(self, IntervalClass, points, n=None):
        self.nonsingular = []
        self.singular = []
        self.empty = [IntervalClass()]
        for i, lhs in enumerate(points):
            for j, rhs in enumerate(points):
                if i < j:
                    self.nonsingular.append(IntervalClass(min=lhs, max=rhs))
                elif i > j:
                    self.empty.append(IntervalClass(min=lhs, max=rhs))
                else:
                    self.singular.append(IntervalClass(min=lhs, max=rhs))
        if n is not None:
            self.nonsingular = self.subset(self.nonsingular, n)
            self.singular = self.subset(self.singular, n)
            self.empty = self.subset(self.empty, n)

    @staticmethod
    def subset(seq, n):
        if len(seq) > n:
            return [seq[i] for i in np.random.choice(len(seq), n)]
        return seq

    @property
    def finite(self):
        return itertools.chain(self.nonsingular, self.singular)

    @property
    def all(self):
        return itertools.chain(self.nonsingular, self.singular, self.empty)


class IntervalTestMixin:

    def assertAllTrue(self, iterable):
        seq = list(iterable)
        self.assertEqual(seq, [True]*len(seq))

    def assertAllFalse(self, iterable):
        seq = list(iterable)
        self.assertEqual(seq, [False]*len(seq))

    def testAccessors(self):
        self.assertEqual(self.ab.min, self.a)
        self.assertEqual(self.ab.max, self.b)

    def testEmpty(self):
        self.assertTrue(self.IntervalClass().isEmpty())
        self.assertAllFalse(s.isEmpty() for s in self.intervals.finite)
        self.assertAllTrue(s.isEmpty() for s in self.intervals.empty)

    def testConstructor(self):
        for i in self.intervals.finite:
            with self.subTest(i=i):
                self.assertEqual(i, self.IntervalClass(min=i.min, max=i.max))
                self.assertEqual(i, self.IntervalClass(min=i.min, size=i.size))
                self.assertEqual(i, self.IntervalClass(max=i.max, size=i.size))

    def testFromHull(self):
        for n1, p1 in enumerate(self.points):
            for n2, p2 in enumerate(self.points):
                with self.subTest(n1=n1, p1=p1, n2=n2, p2=p2):
                    seq = list(self.points[n1:n2+1])
                    i = self.IntervalClass(min=p1, max=p2)
                    self.assertEqual(i, self.IntervalClass.fromHull(seq))
                    np.random.shuffle(seq)
                    self.assertEqual(i, self.IntervalClass.fromHull(seq))
                    seq.reverse()
                    self.assertEqual(i, self.IntervalClass.fromHull(seq))

    def testContains(self):
        for lhs in self.intervals.finite:
            for rhs in self.intervals.finite:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertEqual(lhs.contains(rhs), lhs.min <= rhs.min and lhs.max >= rhs.max)
            for rhs in self.intervals.empty:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertTrue(lhs.contains(rhs))
            for rhs in self.points:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertEqual(lhs.contains(rhs), lhs.min <= rhs and lhs.max >= rhs)
            array = np.array(self.points)
            np.testing.assert_array_equal(lhs.contains(array),
                                          np.logical_and(lhs.min <= array, lhs.max >= array))
        for lhs in self.intervals.empty:
            for rhs in self.intervals.finite:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertFalse(lhs.contains(rhs))
            for rhs in self.intervals.empty:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertTrue(lhs.contains(rhs))
            for rhs in self.points:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertFalse(lhs.contains(rhs))

    def testOverlaps(self):
        for lhs in self.intervals.finite:
            for rhs in self.intervals.finite:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertEqual(lhs.overlaps(rhs),
                                     (lhs.min <= rhs.max and lhs.max >= rhs.min) or
                                     (rhs.min <= lhs.max and rhs.max >= lhs.min) or
                                     lhs.contains(rhs) or rhs.contains(lhs))
            for rhs in self.intervals.empty:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertFalse(lhs.overlaps(rhs))
        for lhs in self.intervals.empty:
            for rhs in self.intervals.all:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertFalse(lhs.overlaps(rhs))

    def testEquality(self):
        for lhs in self.intervals.all:
            for rhs in self.intervals.all:
                with self.subTest(lhs=lhs, rhs=rhs):
                    shouldBeEqual = lhs.contains(rhs) and rhs.contains(lhs)
                    self.assertIs(lhs == rhs, shouldBeEqual)
                    self.assertIs(lhs != rhs, not shouldBeEqual)

    def testClippedTo(self):
        for lhs in self.intervals.all:
            for rhs in self.intervals.all:
                with self.subTest(lhs=lhs, rhs=rhs):
                    clipped = lhs.clippedTo(rhs)
                    self.assertTrue(lhs.contains(clipped))
                    self.assertTrue(rhs.contains(clipped))
                    self.assertIs(
                        clipped.isEmpty(),
                        lhs.isEmpty() or rhs.isEmpty() or not lhs.overlaps(rhs)
                    )
                    self.assertIs(clipped == rhs, lhs.contains(rhs))
                    self.assertIs(clipped == lhs, rhs.contains(lhs))

    def testShiftedBy(self):
        for lhs in self.intervals.finite:
            for rhs in self.points:
                with self.subTest(lhs=lhs, rhs=rhs):
                    shifted = lhs.shiftedBy(rhs)
                    self.assertEqual(lhs.size, shifted.size)
                    self.assertEqual(lhs.min + rhs, shifted.min)
                    self.assertEqual(lhs.max + rhs, shifted.max)
        for lhs in self.intervals.empty:
            for rhs in self.points:
                with self.subTest(lhs=lhs, rhs=rhs):
                    shifted = lhs.shiftedBy(rhs)
                    self.assertTrue(shifted.isEmpty())

    def testExpandedTo(self):
        for lhs in self.intervals.all:
            for rhs in self.intervals.all:
                with self.subTest(lhs=lhs, rhs=rhs):
                    expanded = lhs.expandedTo(rhs)
                    self.assertTrue(expanded.contains(lhs))
                    self.assertTrue(expanded.contains(rhs))
                    self.assertIs(
                        expanded.isEmpty(),
                        lhs.isEmpty() and rhs.isEmpty()
                    )
                    self.assertIs(expanded == rhs, rhs.contains(lhs))
                    self.assertIs(expanded == lhs, lhs.contains(rhs))
            for rhs in self.points:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertEqual(lhs.expandedTo(rhs),
                                     lhs.expandedTo(self.IntervalClass(min=rhs, max=rhs)))

    def testDilatedBy(self):
        for lhs in self.intervals.finite:
            for rhs in self.points:
                with self.subTest(lhs=lhs, rhs=rhs):
                    dilated = lhs.dilatedBy(rhs)
                    if not dilated.isEmpty():
                        self.assertEqual(lhs.min - rhs, dilated.min)
                        self.assertEqual(lhs.max + rhs, dilated.max)
        for lhs in self.intervals.empty:
            for rhs in self.points:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertTrue(lhs.dilatedBy(rhs).isEmpty())

    def testErodedBy(self):
        for lhs in self.intervals.all:
            for rhs in self.points:
                with self.subTest(lhs=lhs, rhs=rhs):
                    self.assertEqual(lhs.erodedBy(rhs), lhs.dilatedBy(-rhs))


class IntervalDTestCase(unittest.TestCase, IntervalTestMixin):
    IntervalClass = IntervalD

    def setUp(self):
        self.points = [-1.5, 5.0, 6.75, 8.625]
        self.a = self.points[0]
        self.b = self.points[1]
        self.intervals = IntervalTestData(self.IntervalClass, self.points, n=3)
        self.ab = self.IntervalClass(min=self.a, max=self.b)

    def testCenter(self):
        self.assertEqual(self.ab.center, 0.5*(self.a + self.b))
        for i in self.intervals.finite:
            self.assertEqual(i, self.IntervalClass(center=i.center, size=i.size))


class IntervalITestCase(unittest.TestCase, IntervalTestMixin):
    IntervalClass = IntervalI

    def setUp(self):
        self.points = [-2, 4, 7, 11]
        self.a = self.points[0]
        self.b = self.points[1]
        self.intervals = IntervalTestData(self.IntervalClass, self.points, n=3)
        self.ab = self.IntervalClass(min=self.a, max=self.b)

    def testExtensions(self):
        s = list(range(10))
        i = IntervalI(min=3, max=8)
        self.assertEqual(s[i.slice], list(i))
        self.assertEqual(len(i), i.size)
        np.testing.assert_array_equal(np.array(list(i), dtype=np.int32), i.arange())
        np.testing.assert_array_equal(np.array(list(i), dtype=np.int64), i.arange(dtype=np.int64))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
