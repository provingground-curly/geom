.. _geom:

####
geom
####

.. _geom-intro:

Introduction
============

Low-level geometry primitives, including Point, Extent, Box,
Angle, SpherePoint, LinearTransform and AffineTransform

.. geom-getting-started:

Getting Started
===============

Sample code showing basic usage of some of the geom classes::

    import math
    import lsst.geom

    # Point, Extent and Box have int and float versions
    # Are are some examples of the integer versions:
    start = lsst.geom.Point2I(0, -1)
    # The end point is not included when constructing an integer box, so...
    dim = lsst.geom.Extent2I(10, 10)
    end = start + dim - lsst.geom.Extent2I(1, 1)
    intBox = lsst.geom.Box2I(start, end)
    assert intBox.getMin() == start
    assert intBox.getMax() == end
    assert intBox.getDimensions() == lsst.geom.Extent2I(10, 10)
    # Boxes can also be constructed from a point and extent
    intBox2 = lsst.geom.Box2I(start, dim)
    assert intBox == intBox2

    # Float boxes represent positions:
    # - the position of the center of pixel 0, 0 is 0.0, 0.0
    # - the position of the lower left corner of pixel 0, 0 is -0.5, -0.5
    floatBox = lsst.geom.Box2D(intBox)
    assert floatBox.getMin() == lsst.geom.Point2D(-0.5, -1.5)
    assert floatBox.getDimensions() == lsst.geom.Extent2D(10, 10)
    assert floatBox.getCenter() == lsst.geom.Point2D(4.5, 3.5)

    # Angles can be constructed by multiplying a float by an AngleUnit,
    # such as `degrees`, `radians` or `arcseconds`.
    # Angles can be added to or subtracted from each other
    # and multiplied or divided by scalars
    rightAngle = 90*lsst.geom.degrees
    assert math.isclose(rightAngle.asRadians(), math.pi/2)
    assert math.isclose((rightAngle - rightAngle).asRadians(), 0)
    assert math.isclose((rightAngle*2).asDegrees(), 180)

    # SpherePoint represents a point on a unit sphere
    # specified by a longitude and latitude (both angles)
    # Unless otherwise noted in code a SpherePoint is ICRS RA, Dec
    # (though the geom package has no support for coordinate systems).
    spherePoint1 = lsst.geom.SpherePoint(10*lsst.geom.degrees, 47*lsst.geom.degrees)
    spherePoint2 = lsst.geom.SpherePoint(10*lsst.geom.degrees, 32*lsst.geom.degrees)
    assert math.isclose(spherePoint1.separation(spherePoint2).asDegrees(), 15)
