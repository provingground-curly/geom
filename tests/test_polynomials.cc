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
#define BOOST_TEST_MODULE polynomials

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop

#include <vector>
#include <typeinfo>
#include <limits>
#include "Eigen/Core"

#include "lsst/geom/polynomials.h"

using namespace lsst::geom::polynomials;
using lsst::geom::Point2D;
using lsst::geom::Box2D;

namespace {

constexpr double DEFAULT_RTOL = 2*std::numeric_limits<double>::epsilon();

// Wrapper around getType that returns a comparable object with operator<< support,
// so we can use it with BOOST_CHECK_EQUAL
template <typename T>
std::string getType(T const & v) {
    return std::string(typeid(T).name());
}

inline bool compare(double a, double b, double rtol) {
    if (!(std::abs(a - b) <= std::sqrt(a*b)*rtol)) {
        std::cerr << (boost::format("a=%0.16g, b=%0.16g, diff=%0.16g > %0.16g")
                      % a % b % std::abs(a - b) % (std::sqrt(a*b)*rtol)) << std::endl;
        return false;
    }
    return true;
}

inline bool compare(Point2D const & a, Point2D const & b, double rtol) {
    return compare(a.getX(), b.getX(), rtol) && compare(a.getY(), b.getY(), rtol);
}

inline bool compare(std::vector<double> const & a, std::vector<double> const & b, double rtol) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i != a.size(); ++i) {
        if (!compare(a[i], b[i], rtol)) {
            return false;
        }
    }
    return true;
}

#define CUSTOM_CHECK_CLOSE(a, b, rtol) \
    do {                               \
        auto a1 = a;                   \
        auto b1 = b;                   \
        BOOST_CHECK( compare(a1, b1, rtol) );   \
    } while (false)

void testBinomialMatrix(int n) {
    auto factorial = [](int m) {
        double result = 1.0;
        for (int k = 1; k <= m; ++k) {
            result *= k;
        }
        return result;
    };
    BinomialMatrix binomial(n);
    for (int k = 0; k <= n; ++k) {
        CUSTOM_CHECK_CLOSE(binomial(n, k), factorial(n)/(factorial(k)*factorial(n - k)), DEFAULT_RTOL);
    }
}

template <typename Basis, typename Point>
void testBasis(Basis const & basis, Point const & point, std::vector<double> const & coefficients) {

    BOOST_CHECK_EQUAL(basis.size(), coefficients.size());

    // Test that we can call sumWith on:
    //  1) STL random-access containers
    double z1 = basis.sumWith(point, coefficients);
    //  2) STL random-access iterators
    double z2 = basis.sumWith(point, coefficients.begin());
    //  3) Eigen objects
    Eigen::VectorXd coefficients3(basis.size());
    std::copy(coefficients.begin(), coefficients.end(), &coefficients3[0]);
    double z3 = basis.sumWith(point, coefficients3);
    //  4) Eigen expressions
    Eigen::VectorXd coefficients4(basis.size() + 1);
    std::copy(coefficients.begin(), coefficients.end(), &coefficients4[0]);
    coefficients4[basis.size()] = std::numeric_limits<double>::quiet_NaN(); // should be ignored
    double z4 = basis.sumWith(point, coefficients4.head(basis.size())*2.0) / 2.0;
    // All of those evaluations should give the same results.
    BOOST_CHECK_EQUAL(z1, z2);
    BOOST_CHECK_EQUAL(z1, z3);
    BOOST_CHECK_EQUAL(z1, z4);

    // Test that we can call fill on the same:
    //  1) STL random-access containers
    std::vector<double> basis1(basis.size(), 0.0);
    basis.fill(point, basis1);
    //  2) STL random-access iterators
    std::vector<double> basis2(basis.size() + 1, 0.0);
    basis2[basis.size()] = std::numeric_limits<double>::quiet_NaN(); // should be ignored
    basis.fill(point, basis2.begin());
    //  3) Eigen objects
    Eigen::VectorXd basis3(basis.size());
    basis.fill(point, basis3);
    //  4) Eigen *view* expressions (can't assign to rvalue expressions, of course)
    Eigen::VectorXd basis4(basis.size() + 1);
    basis4[basis.size()] = std::numeric_limits<double>::quiet_NaN(); // should be ignored
    basis.fill(point, basis4.head(basis.size()));
    // All of those basis vectors should have the same values
    BOOST_CHECK(std::equal(basis1.begin(), basis1.end(), basis2.begin()));
    BOOST_CHECK(std::equal(basis1.begin(), basis1.end(), &basis3[0]));
    BOOST_CHECK(std::equal(basis1.begin(), basis1.end(), &basis4[0]));

    // sum(basis * coefficients) should be equal to sumWith(), subject to round-off error
    CUSTOM_CHECK_CLOSE(basis3.dot(coefficients3), z1, 5*DEFAULT_RTOL);

    // Test using Function object to do the evaluation.
    typename Basis::Function func(basis, coefficients.begin(), coefficients.end());
    double z5 = func(point);
    BOOST_CHECK_EQUAL(z1, z5);

    // Test that we can access Function coefficients various ways.
    auto asEigen = func.getCoefficients();
    std::size_t n = 0;
    for (auto & coeff : func) {
        BOOST_CHECK_EQUAL(asEigen[n], coeff);
        BOOST_CHECK_EQUAL(&asEigen[n], &coeff);
        BOOST_CHECK_EQUAL(func[n], coeff);
        BOOST_CHECK_EQUAL(&func[n], &coeff);
        ++n;
    }
    BOOST_CHECK_EQUAL(n, func.size());

}

template <typename Basis, typename Point, typename Scaling>
void testScaledBasis(Basis const & basis, Point const & point, std::vector<double> const & coefficients,
                     Scaling const & scaling) {
    auto scaledBasis = basis.scaled(scaling);
    auto scaledPoint = scaling.applyForward(point);

    // Run the regular basis tests on the scaled basis.
    testBasis(scaledBasis, scaledPoint, coefficients);

    double z1 = basis.sumWith(scaledPoint, coefficients);
    double z2 = scaledBasis.sumWith(point, coefficients);
    CUSTOM_CHECK_CLOSE(z1, z2, DEFAULT_RTOL);

    std::vector<double> basis1(basis.size(), 0.0);
    std::vector<double> basis2(scaledBasis.size(), 0.0);
    basis.fill(scaledPoint, basis1);
    scaledBasis.fill(point, basis2);

    CUSTOM_CHECK_CLOSE(basis1, basis2, 2*DEFAULT_RTOL);

    // Test using Function object to do the scaling and evaluation.
    typename Basis::Function func(basis, coefficients.begin(), coefficients.end());
    auto tFunc = func.scaled(scaling);
    double z3 = tFunc(point);
    CUSTOM_CHECK_CLOSE(z1, z3, DEFAULT_RTOL);
}

template <PackingOrder packing>
void testPackedIndex() {
    using Range = PackedIndexRange<packing>;
    using Iterator = PackedIndexIterator<packing>;
    int const order = 6;
    Range const range(Iterator(), Iterator::makeEnd(order));
    std::size_t n = 0;
    std::size_t count = 0;
    for (auto const & index : range) {
        BOOST_CHECK_EQUAL(index.flat, Range::computeIndex(index.nx, index.ny));
        BOOST_CHECK(index.nx + index.ny >= n); // order is strictly increasing throughout iteration.
        n = index.nx + index.ny;
        BOOST_CHECK(n <= order);
        if (index.nx == 0u && packing == PackingOrder::YX) {
            BOOST_CHECK_EQUAL(index.flat, Range::computeOffset(n));
        } else if (index.ny == 0u && packing == PackingOrder::XY) {
            BOOST_CHECK_EQUAL(index.flat, Range::computeOffset(n));
        }
        ++count;
    }
    BOOST_CHECK_EQUAL(count, range.size());
    BOOST_CHECK_EQUAL(count, Range::computeSize(order));
}

} // anonymous

BOOST_AUTO_TEST_CASE(PackedIndex) {
    testPackedIndex<PackingOrder::XY>();
    testPackedIndex<PackingOrder::YX>();
}

BOOST_AUTO_TEST_CASE(scalings1d) {
    double scale = 2.0;
    double shift = -0.5;
    Scaling1d affine(scale, shift);
    BOOST_CHECK_EQUAL(affine.getScale(), scale);
    BOOST_CHECK_EQUAL(affine.getShift(), shift);
    auto inverse = affine.inverted();
    auto identity = affine.then(inverse);
    for (double x = -0.5; x < 2; x += 0.3) {
        CUSTOM_CHECK_CLOSE(identity.applyForward(x), x, DEFAULT_RTOL);
        CUSTOM_CHECK_CLOSE(identity.applyInverse(x), x, DEFAULT_RTOL);
        double y = affine.applyForward(x);
        BOOST_CHECK_EQUAL(y, (x + affine.getShift())*affine.getScale());
        CUSTOM_CHECK_CLOSE(affine.applyInverse(y), x, DEFAULT_RTOL);
        CUSTOM_CHECK_CLOSE(inverse.applyForward(y), x, DEFAULT_RTOL);
        CUSTOM_CHECK_CLOSE(inverse.applyInverse(x), y, DEFAULT_RTOL);
    }

    double min=-0.5, max=2.0;
    auto toUnitRange = makeUnitRangeScaling1d(min, max);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyForward(min), -1.0, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyForward(max), 1.0, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(-1.0), min, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(1.0), max, DEFAULT_RTOL);
}


BOOST_AUTO_TEST_CASE(scalings2d) {
    Scaling1d xs(2.0, 0.5);
    Scaling1d ys(-0.5, -1.0);
    Scaling2d affine(xs, ys);
    BOOST_CHECK_EQUAL(affine.getX().getScale(), xs.getScale());
    BOOST_CHECK_EQUAL(affine.getX().getShift(), xs.getShift());
    BOOST_CHECK_EQUAL(affine.getY().getScale(), ys.getScale());
    BOOST_CHECK_EQUAL(affine.getY().getShift(), ys.getShift());
    auto inverse = affine.inverted();
    auto identity = affine.then(inverse);
    for (double y = -2; y < 2; y += 0.3) {
        for (double x = -2; x < 2; x += 0.3) {
            Point2D p(x, y);
            CUSTOM_CHECK_CLOSE(identity.applyForward(p), p, DEFAULT_RTOL);
            CUSTOM_CHECK_CLOSE(identity.applyInverse(p), p, DEFAULT_RTOL);
            auto z = affine.applyForward(p);
            BOOST_CHECK_EQUAL(z, Point2D(affine.getX().applyForward(x),
                                         affine.getY().applyForward(y)));
            CUSTOM_CHECK_CLOSE(affine.applyInverse(z), p, 5*DEFAULT_RTOL);
            CUSTOM_CHECK_CLOSE(inverse.applyForward(z), p, 5*DEFAULT_RTOL);
            CUSTOM_CHECK_CLOSE(inverse.applyInverse(p), z, 5*DEFAULT_RTOL);
        }
    }

    Box2D box(Point2D(-0.5, -2.0), Point2D(2.0, 1.0));
    auto toUnitRange = makeUnitRangeScaling2d(box);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyForward(box.getMin()).getX(), -1.0, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyForward(box.getMin()).getY(), -1.0, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyForward(box.getMax()).getX(), 1.0, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyForward(box.getMax()).getY(), 1.0, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(Point2D(-1.0, -1.0)).getX(), box.getMinX(), DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(Point2D(-1.0, -1.0)).getY(), box.getMinY(), DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(Point2D(-1.0, 1.0)).getX(), box.getMinX(), DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(Point2D(-1.0, 1.0)).getY(), box.getMaxY(), DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(Point2D(1.0, 1.0)).getX(), box.getMaxX(), DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(Point2D(1.0, 1.0)).getY(), box.getMaxY(), DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(Point2D(1.0, -1.0)).getX(), box.getMaxX(), DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(Point2D(1.0, -1.0)).getY(), box.getMinY(), DEFAULT_RTOL);
}

BOOST_AUTO_TEST_CASE(binomials) {
    testBinomialMatrix(3);
    testBinomialMatrix(5);
}

BOOST_AUTO_TEST_CASE(basis1d) {
    std::vector<double> coefficients = { 4.2, 1.6, -3.0, 0.2, -1.1, 0.8 };
    Scaling1d scaling(2.0, -1.0);
    double point = 1.5;
    double min = -0.5, max=2.0;
    int order = 5;

    // Tests below are self-contained: they don't require any particular
    // input point or coefficients vector or any agreement between them.
    // Instead, they just make sure all basis operations on those inputs
    // (or other inputs derived from them) are self-consistent.

    // regular bases
    testBasis(PolynomialBasis1d(order), point, coefficients);
    testBasis(Chebyshev1Basis1d(order), point, coefficients);

    // scaled once
    testScaledBasis(PolynomialBasis1d(order), point, coefficients, scaling);
    testScaledBasis(Chebyshev1Basis1d(order), point, coefficients, scaling);

    // scaled twice (just to test ScaledBasis1d::scaled() and what it returns)
    testScaledBasis(ScaledPolynomialBasis1d(order, min, max), point, coefficients, scaling);
    testScaledBasis(ScaledChebyshev1Basis1d(order, min, max), point, coefficients, scaling);
}

BOOST_AUTO_TEST_CASE(basis2d) {
    std::vector<double> coefficients = { 4.2, 1.6, -3.0, 0.2, -1.1, 0.8 };
    Scaling2d scaling(
        Scaling1d(2.0, -1.0),
        Scaling1d(0.8, 0.6)
    );
    Point2D point(1.5, -0.3);
    Box2D box(Point2D(-4.0, -3.5), Point2D(2.2, 1.8));
    int order = 2;

    // Tests below are self-contained: they don't require any particular
    // input point or coefficients vector or any agreement between them.
    // Instead, they just make sure all basis operations on those inputs
    // (or other inputs derived from them) are self-consistent.

    // regular bases
    testBasis(PolynomialBasis2dXY(order), point, coefficients);
    testBasis(Chebyshev1Basis2dXY(order), point, coefficients);
    testBasis(PolynomialBasis2dYX(order), point, coefficients);
    testBasis(Chebyshev1Basis2dYX(order), point, coefficients);

    // scaled once
    testScaledBasis(PolynomialBasis2dXY(order), point, coefficients, scaling);
    testScaledBasis(Chebyshev1Basis2dXY(order), point, coefficients, scaling);
    testScaledBasis(PolynomialBasis2dYX(order), point, coefficients, scaling);
    testScaledBasis(Chebyshev1Basis2dYX(order), point, coefficients, scaling);

    // scaled twice (just to test ScaledBasis2d::scaled() and what it returns)
    testScaledBasis(ScaledPolynomialBasis2dXY(order, box), point, coefficients, scaling);
    testScaledBasis(ScaledChebyshev1Basis2dXY(order, box), point, coefficients, scaling);
    testScaledBasis(ScaledPolynomialBasis2dYX(order, box), point, coefficients, scaling);
    testScaledBasis(ScaledChebyshev1Basis2dYX(order, box), point, coefficients, scaling);
}

BOOST_AUTO_TEST_CASE(safeSum) {
    SafeSum<double> s;
    s += 1.0;
    s += 1E100;
    s += 1.0;
    s -= 1E100;
    CUSTOM_CHECK_CLOSE(static_cast<double>(s), 2.0, DEFAULT_RTOL);
}

BOOST_AUTO_TEST_CASE(simplified1d) {
    std::vector<double> const coefficients({ 4.2, 1.6, -3.0, 0.2, -1.1, 0.8 });
    double min = -0.5, max=2.0;
    auto sfunc = ScaledPolynomialFunction1d(
        ScaledPolynomialBasis1d(5, min, max),
        coefficients.begin(), coefficients.end()
    );
    auto func = simplified(sfunc);
    for (double x = min; x < max; x += 0.3) {
        CUSTOM_CHECK_CLOSE(sfunc(x), func(x), 5*DEFAULT_RTOL);
    }
}

BOOST_AUTO_TEST_CASE(simplified2d) {
    std::vector<double> const coefficients({
        4.2,
        1.6, -3.0,
        0.2, -1.1, 0.8,
        1.2, 0.7, 1.9, -0.6,
        5.0, 7.2, -9.1, -8.5, 0.0
    });
    Box2D box(Point2D(-4.0, -3.5), Point2D(2.2, 1.8));
    auto sfunc = ScaledPolynomialFunction2dXY(
        ScaledPolynomialBasis2dXY(4, box),
        coefficients.begin(), coefficients.end()
    );
    auto func = simplified(sfunc);
    for (double x = box.getMinX(); x < box.getMaxX(); x += 0.3) {
        for (double y = box.getMinY(); y < box.getMaxY(); y += 0.3) {
            Point2D point(x, y);
            // First check: safe summation in simplified() itself (always
            // present), but not when evaluating the polynomials.
            // That means we still lose precision in evaluating the polynomials
            // in a way that depends on the relative values in the sum.
            // Tolerance is selected to work for these test values, so think
            // of this as a regression test, not something fundamental.
            CUSTOM_CHECK_CLOSE(sfunc(point), func(point), 100*DEFAULT_RTOL);
            // Second check: safe summation in simplifed() itself *and*
            // evaluating the polynomials.  We still lose a fair amount of
            // precision to operations other than the sum.
            CUSTOM_CHECK_CLOSE(sfunc(point, SumMode::SAFE), func(point, SumMode::SAFE), 50*DEFAULT_RTOL);
        }
    }
}
