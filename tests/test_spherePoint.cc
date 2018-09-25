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
#define BOOST_TEST_MODULE SpherePointCpp

#include "boost/test/unit_test.hpp"

#include "lsst/utils/tests.h"

#include "lsst/sphgeom/Vector3d.h"
#include "lsst/geom/SpherePoint.h"
#include "lsst/pex/exceptions/Exception.h"

/*
 * Unit tests for C++-only functionality in SpherePoint.
 *
 * See test_spherePoint.py for remaining unit tests.
 */
namespace lsst {
namespace geom {

/*
 * Tests whether the result of SpherePoint::SpherePoint(SpherePoint const&)
 * is an identical but independent copy.
 */
BOOST_AUTO_TEST_CASE(SpherePointCopyResult, *boost::unit_test::tolerance(1e-14)) {
    SpherePoint original(sphgeom::Vector3d(0.34, -1.2, 0.97));
    SpherePoint copy(original);

    // Want exact equality, not floating-point equality, for results of copy-construction
    BOOST_TEST(original == copy);

    // Don't compare Angles in case there's aliasing of some sort
    double const copyLon = copy.getLongitude().asDegrees();
    double const copyLat = copy.getLatitude().asDegrees();

    original = SpherePoint(-42 * degrees, 45 * degrees);
    BOOST_TEST(original != copy);
    BOOST_TEST(copy.getLongitude().asDegrees() == copyLon);
    BOOST_TEST(copy.getLatitude().asDegrees() == copyLat);
}

/*
 * Tests whether the result of SpherePoint::SpherePoint(SpherePoint&&)
 * is an identical copy.
 */
BOOST_AUTO_TEST_CASE(SpherePointMoveResult, *boost::unit_test::tolerance(1e-14)) {
    SpherePoint original(sphgeom::Vector3d(0.34, -1.2, 0.97));
    // Don't compare Angles in case there's aliasing of some sort
    double const oldLon = original.getLongitude().asDegrees();
    double const oldLat = original.getLatitude().asDegrees();

    SpherePoint copy(std::move(original));

    BOOST_TEST(copy.getLongitude().asDegrees() == oldLon);
    BOOST_TEST(copy.getLatitude().asDegrees() == oldLat);
}

/*
 * Tests whether SpherePoint::operator=(SpherePoint const&) makes an identical
 * but independent copy.
 */
BOOST_AUTO_TEST_CASE(assignCopyResult, *boost::unit_test::tolerance(1e-14)) {
    SpherePoint original(sphgeom::Vector3d(0.34, -1.2, 0.97));
    // Don't compare Angles in case there's aliasing of some sort
    double const oldLon = original.getLongitude().asDegrees();
    double const oldLat = original.getLatitude().asDegrees();

    SpherePoint copy(45.0 * degrees, -23.5 * degrees);
    // Want exact equality, not floating-point equality, for results of assignment
    BOOST_REQUIRE(original != copy);
    copy = original;
    BOOST_TEST(original == copy);

    original = SpherePoint(-42 * degrees, 45 * degrees);
    BOOST_TEST(original != copy);
    BOOST_TEST(copy.getLongitude().asDegrees() == oldLon);
    BOOST_TEST(copy.getLatitude().asDegrees() == oldLat);
}

/*
 * Tests whether SpherePoint::operator=(SpherePoint const&) makes an identical
 * copy.
 */
BOOST_AUTO_TEST_CASE(assignMoveResult, *boost::unit_test::tolerance(1e-14)) {
    SpherePoint original(sphgeom::Vector3d(0.34, -1.2, 0.97));
    // Don't compare Angles in case there's aliasing of some sort
    double const oldLon = original.getLongitude().asDegrees();
    double const oldLat = original.getLatitude().asDegrees();

    SpherePoint copy(45.0 * degrees, -23.5 * degrees);
    BOOST_REQUIRE(original != copy);
    copy = std::move(original);

    BOOST_TEST(copy.getLongitude().asDegrees() == oldLon);
    BOOST_TEST(copy.getLatitude().asDegrees() == oldLat);
}

/*
 * Tests whether SpherePoint::operator[](size_t) handles invalid indices
 * correctly.
 */
BOOST_AUTO_TEST_CASE(getItemError) {
    SpherePoint point(sphgeom::Vector3d(1.0, 1.0, 1.0));

    BOOST_CHECK_THROW(point[2], pex::exceptions::OutOfRangeError);
    BOOST_CHECK_THROW(point[-1], pex::exceptions::OutOfRangeError);
}

// TODO: add a test for propagation of ostream errors

BOOST_AUTO_TEST_CASE(Hash) {
    utils::assertValidHash<SpherePoint>();

    utils::assertHashesEqual(SpherePoint(0 * degrees, -24 * degrees),
                             SpherePoint(0 * degrees, -24 * degrees));
    utils::assertHashesEqual(SpherePoint(HALFPI * radians, 0.0 * radians),
                             SpherePoint(sphgeom::Vector3d(0, 1, 0)));
}

}  // namespace geom
}  // namespace lsst
