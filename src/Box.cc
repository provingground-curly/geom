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

#include <cmath>
#include <limits>

#include "lsst/utils/hashCombine.h"
#include "lsst/geom/Box.h"

namespace lsst {
namespace geom {

Box2I::Box2I(Point2I const& minimum, Point2I const& maximum, bool invert)
        : _minimum(minimum), _dimensions(maximum - minimum) {
    for (int n = 0; n < 2; ++n) {
        if (_dimensions[n] < 0) {
            if (invert) {
                _minimum[n] += _dimensions[n];
                _dimensions[n] = -_dimensions[n];
            } else {
                *this = Box2I();
                return;
            }
        }
    }
    _dimensions += Extent2I(1);
}

Box2I::Box2I(Point2I const& corner, Extent2I const& dimensions, bool invert)
        : _minimum(corner), _dimensions(dimensions) {
    for (int n = 0; n < 2; ++n) {
        if (_dimensions[n] == 0) {
            *this = Box2I();
            return;
        } else if (_dimensions[n] < 0) {
            if (invert) {
                _minimum[n] += (_dimensions[n] + 1);
                _dimensions[n] = -_dimensions[n];
            } else {
                *this = Box2I();
                return;
            }
        }
    }
    if (!isEmpty() && any(getMin().gt(getMax()))) {
        throw LSST_EXCEPT(pex::exceptions::OverflowError,
                          "Box dimensions too large; integer overflow detected.");
    }
}

Box2I::Box2I(Box2D const& other, EdgeHandlingEnum edgeHandling) : _minimum(), _dimensions() {
    if (other.isEmpty()) {
        *this = Box2I();
        return;
    }
    if (!std::isfinite(other.getMinX()) || !std::isfinite(other.getMinY()) ||
        !std::isfinite(other.getMaxX()) || !std::isfinite(other.getMaxY())) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Cannot convert non-finite Box2D to Box2I");
    }
    Point2D fpMin(other.getMin() + Extent2D(0.5));
    Point2D fpMax(other.getMax() - Extent2D(0.5));
    switch (edgeHandling) {
        case EXPAND:
            for (int n = 0; n < 2; ++n) {
                _minimum[n] = static_cast<int>(std::floor(fpMin[n]));
                _dimensions[n] = static_cast<int>(std::ceil(fpMax[n])) + 1 - _minimum[n];
            }
            break;
        case SHRINK:
            for (int n = 0; n < 2; ++n) {
                _minimum[n] = static_cast<int>(std::ceil(fpMin[n]));
                _dimensions[n] = static_cast<int>(std::floor(fpMax[n])) + 1 - _minimum[n];
            }
            break;
    }
}

Box2I Box2I::makeCenteredBox(Point2D const& center, Box2I::Extent const& size) {
    if (!std::isfinite(center[0]) || !std::isfinite(center[1])) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Cannot make Box2I with non-finite center");
    }

    lsst::geom::Point2D corner(center);
    corner.shift(-0.5 * lsst::geom::Extent2D(size));
    // compensate for Box2I's coordinate conventions (where max = min + size - 1)
    corner.shift(lsst::geom::Extent2D(0.5, 0.5));
    return lsst::geom::Box2I(lsst::geom::Point2I(corner), size, false);
}

ndarray::View<boost::fusion::vector2<ndarray::index::Range, ndarray::index::Range> > Box2I::getSlices()
        const {
    return ndarray::view(getBeginY(), getEndY())(getBeginX(), getEndX());
}

bool Box2I::contains(Point2I const& point) const noexcept {
    return getX().contains(point.getX()) && getY().contains(point.getY());
}

bool Box2I::contains(Box2I const& other) const noexcept {
    return getX().contains(other.getX()) && getY().contains(other.getY());
}

bool Box2I::overlaps(Box2I const& other) const noexcept {
    return !isDisjointFrom(other);
}

bool Box2I::isDisjointFrom(Box2I const& other) const noexcept {
    return getX().isDisjointFrom(other.getX()) || getY().isDisjointFrom(other.getY());
}

void Box2I::grow(Extent2I const& buffer) {
    dilateBy(buffer);
}

void Box2I::shift(Extent2I const& offset) {
    shiftBy(offset);
}

void Box2I::flipLR(int xextent) {
    reflectAboutX(xextent - 1);
}

void Box2I::flipTB(int yextent) {
    reflectAboutY(yextent - 1);
}

void Box2I::include(Point2I const& point) {
    if (isEmpty()) {
        _minimum = point;
        _dimensions = Extent2I(1);
        return;
    }
    Point2I maximum(getMax());
    for (int n = 0; n < 2; ++n) {
        if (point[n] < _minimum[n]) {
            _minimum[n] = point[n];
        } else if (point[n] > maximum[n]) {
            maximum[n] = point[n];
        }
    }
    _dimensions = Extent2I(1) + maximum - _minimum;
}

void Box2I::include(Box2I const& other) {
    if (other.isEmpty()) return;
    if (this->isEmpty()) {
        *this = other;
        return;
    }
    Point2I maximum(getMax());
    Point2I const& otherMin = other.getMin();
    Point2I const otherMax = other.getMax();
    for (int n = 0; n < 2; ++n) {
        if (otherMin[n] < _minimum[n]) {
            _minimum[n] = otherMin[n];
        }
        if (otherMax[n] > maximum[n]) {
            maximum[n] = otherMax[n];
        }
    }
    _dimensions = Extent2I(1) + maximum - _minimum;
}

void Box2I::clip(Box2I const& other) noexcept {
    if (isEmpty()) return;
    if (other.isEmpty()) {
        *this = Box2I();
        return;
    }
    Point2I maximum(getMax());
    Point2I const& otherMin = other.getMin();
    Point2I const otherMax = other.getMax();
    for (int n = 0; n < 2; ++n) {
        if (otherMin[n] > _minimum[n]) {
            _minimum[n] = otherMin[n];
        }
        if (otherMax[n] < maximum[n]) {
            maximum[n] = otherMax[n];
        }
    }
    if (any(maximum.lt(_minimum))) {
        *this = Box2I();
        return;
    }
    _dimensions = Extent2I(1) + maximum - _minimum;
}

Box2I Box2I::dilatedBy(Extent const& buffer) const {
    return Box2I(getX().dilatedBy(buffer.getX()), getY().dilatedBy(buffer.getY()));
}

Box2I Box2I::shiftedBy(Extent const& offset) const {
    return Box2I(getX().shiftedBy(offset.getX()), getY().shiftedBy(offset.getY()));
}

Box2I Box2I::reflectedAboutX(Element x) const {
    return Box2I(getX().reflectedAbout(x), getY());
}

Box2I Box2I::reflectedAboutY(Element y) const {
    return Box2I(getX(), getY().reflectedAbout(y));
}

Box2I Box2I::expandedTo(Point const & other) const {
    return Box2I(getX().expandedTo(other.getX()), getY().expandedTo(other.getY()));
}

Box2I Box2I::expandedTo(Box2I const & other) const {
    return Box2I(getX().expandedTo(other.getX()), getY().expandedTo(other.getY()));
}

Box2I & Box2I::clipTo(Box2I const& other) noexcept {
    // Existing clip() is noexcept, so there's no advantage to delegating to
    // the interval implementation as long as Box2I's data members (min, dims)
    // instead of (x, y) intervals.
    clip(other);
    return *this;
}

bool Box2I::operator==(Box2I const& other) const noexcept {
    return other._minimum == this->_minimum && other._dimensions == this->_dimensions;
}

bool Box2I::operator!=(Box2I const& other) const noexcept {
    return other._minimum != this->_minimum || other._dimensions != this->_dimensions;
}

std::size_t Box2I::hash_value() const noexcept {
    // Completely arbitrary seed
    return utils::hashCombine(17, _minimum, _dimensions);
}

std::vector<Point2I> Box2I::getCorners() const {
    std::vector<Point2I> retVec;
    retVec.push_back(getMin());
    retVec.push_back(Point2I(getMaxX(), getMinY()));
    retVec.push_back(getMax());
    retVec.push_back(Point2I(getMinX(), getMaxY()));
    return retVec;
}

double const Box2D::EPSILON = std::numeric_limits<double>::epsilon() * 2;

double const Box2D::INVALID = std::numeric_limits<double>::quiet_NaN();

Box2D::Box2D() noexcept : _minimum(INVALID), _maximum(INVALID) {}

Box2D::Box2D(Point2D const& minimum, Point2D const& maximum, bool invert) noexcept
        : _minimum(minimum), _maximum(maximum) {
    for (int n = 0; n < 2; ++n) {
        if (_minimum[n] == _maximum[n]) {
            *this = Box2D();
            return;
        } else if (_minimum[n] > _maximum[n]) {
            if (invert) {
                std::swap(_minimum[n], _maximum[n]);
            } else {
                *this = Box2D();
                return;
            }
        }
    }
}

Box2D::Box2D(Point2D const& corner, Extent2D const& dimensions, bool invert) noexcept
        : _minimum(corner), _maximum(corner + dimensions) {
    for (int n = 0; n < 2; ++n) {
        if (_minimum[n] == _maximum[n]) {
            *this = Box2D();
            return;
        } else if (_minimum[n] > _maximum[n]) {
            if (invert) {
                std::swap(_minimum[n], _maximum[n]);
            } else {
                *this = Box2D();
                return;
            }
        }
    }
}

Box2D::Box2D(Box2I const& other) noexcept
        : _minimum(Point2D(other.getMin()) - Extent2D(0.5)),
          _maximum(Point2D(other.getMax()) + Extent2D(0.5)) {
    if (other.isEmpty()) *this = Box2D();
}

Box2D Box2D::makeCenteredBox(Point2D const& center, Box2D::Extent const& size) noexcept {
    lsst::geom::Point2D corner(center);
    corner.shift(-0.5 * size);
    return lsst::geom::Box2D(corner, size, false);
}

bool Box2D::contains(Point2D const& point) const noexcept {
    // Can't delegate to IntervalD here because IntervalID is closed while
    // Box2D is half-open.
    return all(point.ge(this->getMin())) && all(point.lt(this->getMax()));
}

bool Box2D::contains(Box2D const& other) const noexcept {
    return getX().contains(other.getX()) && getY().contains(other.getY());
}

bool Box2D::overlaps(Box2D const& other) const noexcept {
    // Can't delegate to IntervalD here because IntervalID is closed while
    // Box2D is half-open.
    return !(other.isEmpty() || this->isEmpty() || any(other.getMax().le(this->getMin())) ||
             any(other.getMin().ge(this->getMax())));
}

bool Box2D::isDisjointFrom(Box2D const& other) const noexcept {
    return !overlaps(other);
}

void Box2D::grow(Extent2D const& buffer) {
    if (isEmpty()) return;  // should we throw an exception here instead of a no-op?
    _minimum -= buffer;
    _maximum += buffer;
    if (any(_minimum.ge(_maximum))) *this = Box2D();
}

void Box2D::shift(Extent2D const& offset) {
    if (isEmpty()) return;  // should we throw an exception here instead of a no-op?
    _minimum += offset;
    _maximum += offset;
}

void Box2D::flipLR(float xextent) {
    reflectAboutX(xextent);
}

void Box2D::flipTB(float yextent) {
    reflectAboutY(yextent);
}

void Box2D::include(Point2D const& point) noexcept {
    if (isEmpty()) {
        _minimum = point;
        _maximum = point;
        _tweakMax(0);
        _tweakMax(1);
        return;
    }
    for (int n = 0; n < 2; ++n) {
        if (point[n] < _minimum[n]) {
            _minimum[n] = point[n];
        } else if (point[n] >= _maximum[n]) {
            _maximum[n] = point[n];
            _tweakMax(n);
        }
    }
}

void Box2D::include(Box2D const& other) noexcept {
    if (other.isEmpty()) return;
    if (this->isEmpty()) {
        *this = other;
        return;
    }
    Point2D const& otherMin = other.getMin();
    Point2D const& otherMax = other.getMax();
    for (int n = 0; n < 2; ++n) {
        if (otherMin[n] < _minimum[n]) {
            _minimum[n] = otherMin[n];
        }
        if (otherMax[n] > _maximum[n]) {
            _maximum[n] = otherMax[n];
        }
    }
}

void Box2D::clip(Box2D const& other) noexcept {
    if (isEmpty()) return;
    if (other.isEmpty()) {
        *this = Box2D();
        return;
    }
    Point2D const& otherMin = other.getMin();
    Point2D const& otherMax = other.getMax();
    for (int n = 0; n < 2; ++n) {
        if (otherMin[n] > _minimum[n]) {
            _minimum[n] = otherMin[n];
        }
        if (otherMax[n] < _maximum[n]) {
            _maximum[n] = otherMax[n];
        }
    }
    if (any(_maximum.le(_minimum))) {
        *this = Box2D();
        return;
    }
}

Box2D & Box2D::dilateBy(Extent const & buffer) noexcept {
    // No advantage to delegating to IntervalD until Box2D is represented
    // in terms of x and y intervals.
    grow(buffer);
    return *this;
}

Box2D & Box2D::shiftBy(Extent const & offset) noexcept {
    // No advantage to delegating to IntervalD until Box2D is represented
    // in terms of x and y intervals.
    shift(offset);
    return *this;
}

Box2D & Box2D::reflectAboutX(Element x) noexcept {
    *this = Box2D(getX().reflectedAbout(x), getY());
    return *this;
}

Box2D & Box2D::reflectAboutY(Element y) noexcept {
    *this = Box2D(getX(), getY().reflectedAbout(y));
    return *this;
}

Box2D & Box2D::expandTo(Point const & other) noexcept {
    // Can't delegate to IntervalD here because IntervalID is closed while
    // Box2D is half-open.
    include(other);
    return *this;
}

Box2D & Box2D::expandTo(Box2D const & other) noexcept {
    // Can't delegate to IntervalD here because IntervalID is closed while
    // Box2D is half-open.
    include(other);
    return *this;
}

Box2D & Box2D::clipTo(Box2D const& other) noexcept {
    // No advantage to delegating to IntervalD until Box2D is represented
    // in terms of x and y intervals.
    clip(other);
    return *this;
}

bool Box2D::operator==(Box2D const& other) const noexcept {
    return (other.isEmpty() && this->isEmpty()) ||
           (other._minimum == this->_minimum && other._maximum == this->_maximum);
}

bool Box2D::operator!=(Box2D const& other) const noexcept {
    return !(other.isEmpty() && other.isEmpty()) &&
           (other._minimum != this->_minimum || other._maximum != this->_maximum);
}

std::size_t Box2D::hash_value() const noexcept {
    if (isEmpty()) {
        // All empty boxes are equal and must have equal hashes
        return 179;
    } else {
        // Completely arbitrary seed
        return utils::hashCombine(17, _minimum, _maximum);
    }
}

std::vector<Point2D> Box2D::getCorners() const {
    std::vector<Point2D> retVec;
    retVec.push_back(getMin());
    retVec.push_back(Point2D(getMaxX(), getMinY()));
    retVec.push_back(getMax());
    retVec.push_back(Point2D(getMinX(), getMaxY()));
    return retVec;
}

std::ostream& operator<<(std::ostream& os, Box2I const& box) {
    if (box.isEmpty()) return os << "Box2I()";
    return os << "Box2I(Point2I" << box.getMin() << ", Extent2I" << box.getDimensions() << ")";
}

std::ostream& operator<<(std::ostream& os, Box2D const& box) {
    if (box.isEmpty()) return os << "Box2D()";
    return os << "Box2D(Point2D" << box.getMin() << ", Extent2D" << box.getDimensions() << ")";
}

}  // namespace geom
}  // namespace lsst
