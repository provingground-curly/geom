// -*- LSST-C++ -*-
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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Scaling2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Scaling2d_h_INCLUDED

#include "lsst/geom/polynomials/Scaling1d.h"
#include "lsst/geom/Point.h"
#include "lsst/geom/Box.h"

namespace lsst { namespace geom { namespace polynomials {

/**
 *  A 2-d separable affine transform that can be used to map one interval to another.
 *
 *  The transform is represented in each dimension as an additive shift followed by a
 *  multiplicative scaling.  Unlike a full affine transform, Scaling2d cannot
 *  include rotations.
 *
 *  @note This class (and to a lesser extent its 1-d counterpart, Scaling1d)
 *        has a lot in common with geom::AffineTransform, and ideally
 *        they should share code and be highly interoperable.  Doing that well
 *        would require a larger-scale rethink of geom, however, and at
 *        present we don't actually have a use case for that interoperability,
 *        so it's something we should keep in mind for the future, not a high
 *        priority for the present.
 *
 *  @see makeUnitRangeScaling2d
 */
class Scaling2d {
public:

    /// Construct from the given 1-d scalings.
    Scaling2d(Scaling1d const & x, Scaling1d const & y) noexcept : _x(x), _y(y) {}

    /// Default copy constructor.
    Scaling2d(Scaling2d const &) noexcept = default;

    /// Default move constructor.
    Scaling2d(Scaling2d &&) noexcept = default;

    /// Default copy assignment.
    Scaling2d & operator=(Scaling2d const &) noexcept = default;

    /// Default move assignment.
    Scaling2d & operator=(Scaling2d &&) noexcept = default;

    /// Return the 1-d scaling in the X direction.
    Scaling1d const & getX() const noexcept { return _x; }

    /// Return the 1-d scaling in the Y direction.
    Scaling1d const & getY() const noexcept { return _y; }

    /// Apply the transform in the forward direction.
    geom::Point2D applyForward(geom::Point2D const & p) const noexcept {
        return geom::Point2D(getX().applyForward(p.getX()), getY().applyForward(p.getY()));
    }

    /// Apply the inverse of the forward transform.
    geom::Point2D applyInverse(geom::Point2D const & p) const noexcept {
        return geom::Point2D(getX().applyInverse(p.getX()), getY().applyInverse(p.getY()));
    }

    /**
     *  Invert the transform.
     *
     *  If `r = t.inverted()`, then `r.applyForward(p)` is equivalent to
     *  `t.applyInverse(p)` and `r.applyInverse(q)` is equivalent to
     *  `t.applyForward(q)`.
     */
    Scaling2d inverted() const noexcept {
        return Scaling2d(getX().inverted(), getY().inverted());
    }

    /**
     *  Compose two transforms.
     *
     *  If `r = a.then(b)`, then `r.applyForward(p)` is equivalent to
     *  `b.applyForward(a.applyForward(p))` and `r.applyInverse(q)` is
     *  equivalent to `a.applyInverse(b.applyInverse(q))`.
     */
    Scaling2d then(Scaling2d const & second) const noexcept {
        return Scaling2d(getX().then(second.getX()), getY().then(second.getY()));
    }

private:
    Scaling1d _x;
    Scaling1d _y;
};

/**
 *  Return a Scaling1d that maps the given box to [-1, 1]x[-1, 1].
 */
inline Scaling2d makeUnitRangeScaling2d(geom::Box2D const & box) {
    return Scaling2d(
        makeUnitRangeScaling1d(box.getMinX(), box.getMaxX()),
        makeUnitRangeScaling1d(box.getMinY(), box.getMaxY())
    );
}

}}} // namespace lsst::geom::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Scaling2d_h_INCLUDED
