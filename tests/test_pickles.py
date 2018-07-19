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

"""
Tests for pickles of some afw types
"""

import unittest
import pickle

import lsst.utils.tests
import lsst.geom


class PickleBase:
    """A test case for pickles"""

    def setUp(self):
        raise NotImplementedError(
            "Need to inherit and create the 'data' element.")

    def tearDown(self):
        del self.data

    def assertPickled(self, new):
        """Assert that the pickled data is the same as the original

        Subclasses should override this method if the particular data
        doesn't support the == operator.
        """
        self.assertEqual(new, self.data)

    def testPickle(self):
        """Test round-trip pickle"""
        pickled = pickle.dumps(self.data)
        newData = pickle.loads(pickled)
        self.assertPickled(newData)


class AngleTestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        self.data = 1.0*lsst.geom.degrees


class Point2DTestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        x, y = 1.0, 1.0
        self.data = lsst.geom.Point2D(x, y)


class Point2ITestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        x, y = 1, 1
        self.data = lsst.geom.Point2I(x, y)


class Point3DTestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        x, y, z = 1.0, 1.0, 1.0
        self.data = lsst.geom.Point3D(x, y, z)


class Point3ITestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        x, y, z = 1, 1, 1
        self.data = lsst.geom.Point3I(x, y, z)


class Extent2DTestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        x, y = 1.0, 1.0
        self.data = lsst.geom.Extent2D(x, y)


class Extent3DTestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        x, y, z = 1, 1, 1
        self.data = lsst.geom.Extent3D(x, y, z)


class Extent2ITestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        x, y = 1, 1
        self.data = lsst.geom.Extent2I(x, y)


class Extent3ITestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        x, y, z = 1, 1, 1
        self.data = lsst.geom.Extent3I(x, y, z)


class Box2DTestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        p, e = lsst.geom.Point2D(1.0, 1.0), lsst.geom.Extent2D(0.5, 0.5)
        self.data = lsst.geom.Box2D(p, e, invert=False)


class Box2ITestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        p, e = lsst.geom.Point2I(1, 2), lsst.geom.Extent2I(1, 1)
        self.data = lsst.geom.Box2I(p, e, invert=False)


class AffineTransformTestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        scale = 2.2
        linear = lsst.geom.LinearTransform().makeScaling(scale)
        dx, dy = 1.1, 3.3
        trans = lsst.geom.Extent2D(dx, dy)
        self.data = lsst.geom.AffineTransform(linear, trans)

    def assertPickled(self, new):
        self.assertListEqual(new.getMatrix().flatten().tolist(),
                             self.data.getMatrix().flatten().tolist())


class LinearTransformTestCase(PickleBase, unittest.TestCase):

    def setUp(self):
        scale = 2.0
        self.data = lsst.geom.LinearTransform().makeScaling(scale)

    def assertPickled(self, new):
        self.assertListEqual(new.getMatrix().flatten().tolist(),
                             self.data.getMatrix().flatten().tolist())


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
