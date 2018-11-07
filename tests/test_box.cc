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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BoxCpp

#include "boost/test/unit_test.hpp"

#include "lsst/utils/tests.h"

#include "lsst/geom/Box.h"

/*
 * Unit tests for C++-only functionality in Box2I and Box2D.
 *
 * See test_box.py for remaining unit tests.
 */
namespace lsst {
namespace geom {

BOOST_AUTO_TEST_CASE(Box2IHash) {
    using Int = Box2I::Element;
    using Point = Box2I::Point;

    utils::assertValidHash<Box2I>();
    utils::assertHashesEqual(Box2I(Point(0, 0), Point(24, 200)), Box2I(Point(0, 0), Extent2I(25, 201)));
    utils::assertHashesEqual(Box2I(), Box2I(Point(24, 200), Point(0, 0), false));
}

BOOST_AUTO_TEST_CASE(Box2DHash) {
    using Real = Box2D::Element;
    using Point = Box2D::Point;
    static const Real nan = std::numeric_limits<Real>::quiet_NaN();

    utils::assertValidHash<Box2D>();
    utils::assertHashesEqual(Box2D(Point(0, 0), Point(24, 20.5)), Box2D(Point(0, 0), Extent2D(24, 20.5)));
    utils::assertHashesEqual(Box2D(), Box2D(Point(24, 20.5), Point(0, 0), false));
    utils::assertHashesEqual(Box2D(), Box2D(Point(nan), Point(42.0, 52.0)));
}

}  // namespace geom
}  // namespace lsst
