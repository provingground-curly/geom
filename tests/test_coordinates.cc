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
#define BOOST_TEST_MODULE ellipses
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop
#include "boost/format.hpp"

#include "lsst/geom/Point.h"
#include "lsst/geom/Extent.h"

// Most Extent and Point operators are tested in coordinates.py in Python, but
// division of negative integers has different behavior in the two languages,
// so we test the C++ version here.
BOOST_AUTO_TEST_CASE(Operators) {
    lsst::geom::Extent2I e1(12, -23);
    BOOST_CHECK_EQUAL(e1 / 4, lsst::geom::Extent2I(e1.getX() / 4, e1.getY() / 4));
    lsst::geom::Extent2I e2(e1);
    e2 /= 3;
    BOOST_CHECK_EQUAL(e2, lsst::geom::Extent2I(e1.getX() / 3, e1.getY() / 3));
}
