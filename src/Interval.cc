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
#include "lsst/geom/Interval.h"

namespace lsst {
namespace geom {

IntervalI IntervalI::fromMinMax(Element min, Element max) {
    BigElement size = 1 + static_cast<BigElement>(max) - static_cast<BigElement>(min);
    if (size <= 0) {
        return IntervalI();
    }
    return IntervalI(min, size);
}

IntervalI IntervalI::fromMinSize(Element min, Element size) {
    if (size == 0) {
        return IntervalI();
    }
    if (size < 0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("size (%s) < 0; use IntervalI() or size=0 "
                           "to construct an empty interval") % size).str()
        );
    }
    return IntervalI(min, size);
}

IntervalI IntervalI::fromMaxSize(Element max, Element size) {
    if (size == 0) {
        return IntervalI();
    }
    if (size < 0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("size (%s) < 0; use IntervalI() or size=0 "
                           "to construct an empty interval") % size).str()
        );
    }
    BigElement min = static_cast<BigElement>(max) - static_cast<BigElement>(size) + 1;
    return IntervalI(min, size);
}

IntervalI IntervalI::fromCenterSize(double center, Element size) {
    if (size == 0) {
        return IntervalI();
    }
    if (size < 0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("size (%s) < 0; use IntervalI() or size=0 "
                           "to construct an empty interval") % size).str()
        );
    }
    if (!std::isfinite(center)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Cannot make IntervalI with non-finite center");
    }
    double min = center + -0.5*size;
    // compensate for IntervalI's coordinate conventions (where max = min + size - 1)
    min += 0.5;
    return IntervalI(min, size);
}

IntervalI::IntervalI(IntervalD const& other, EdgeHandlingEnum edgeHandling) : _min(), _size() {
    if (other.isEmpty()) {
        *this = IntervalI();
        return;
    }
    if (!std::isfinite(other.getMin()) || !std::isfinite(other.getMax())) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Cannot convert non-finite IntervalD to IntervalI");
    }
    double fpMin(other.getMin() + 0.5);
    double fpMax(other.getMax() - 0.5);
    switch (edgeHandling) {
        case EdgeHandlingEnum::EXPAND:
            _min = static_cast<int>(std::floor(fpMin));
            _size = static_cast<int>(std::ceil(fpMax)) + 1 - _min;
            break;
        case EdgeHandlingEnum::SHRINK:
            _min = static_cast<int>(std::ceil(fpMin));
            _size = static_cast<int>(std::floor(fpMax)) + 1 - _min;
            break;
    }
}

ndarray::View<boost::fusion::vector1<ndarray::index::Range>> IntervalI::getSlice() const {
    return ndarray::view(getBegin(), getEnd());
}

bool IntervalI::contains(Element point) const noexcept {
    return point >= this->getMin() && point <= this->getMax();
}

bool IntervalI::contains(IntervalI const& other) const noexcept {
    return other.isEmpty() ||
           (other.getMin() >= this->getMin() && other.getMax() <= this->getMax());
}

bool IntervalI::overlaps(IntervalI const& other) const noexcept {
    return !isDisjointFrom(other);
}

bool IntervalI::isDisjointFrom(IntervalI const & other) const noexcept {
    if (isEmpty() || other.isEmpty()) {
        return true;
    }
    return getMin() > other.getMax() || getMax() < other.getMin();
}

IntervalI IntervalI::dilatedBy(Element buffer) const {
    if (isEmpty()) {
        return IntervalI();
    }
    BigElement min = static_cast<BigElement>(getMin()) - buffer;
    BigElement max = static_cast<BigElement>(getMax()) + buffer;
    if (max < min) {
        return IntervalI();
    }
    return fromMinMax(min, max);
}

IntervalI IntervalI::shiftedBy(Element offset) const{
    if (isEmpty()) {
        return IntervalI();
    }
    BigElement min = static_cast<BigElement>(getMin()) + offset;
    return fromMinSize(min, getSize());
}

IntervalI IntervalI::reflectedAbout(Element point) const {
    if (isEmpty()) {
        return IntervalI();
    }
    BigElement max = static_cast<BigElement>(point) - getMin();
    BigElement min = static_cast<BigElement>(point) - getMax();
    return fromMinMax(min, max);
}

IntervalI IntervalI::expandedTo(Element point) const {
    if (isEmpty()) {
        return IntervalI(point, 1);
    }
    return fromMinMax(std::min(static_cast<BigElement>(point),
                               static_cast<BigElement>(getMin())),
                      std::max(static_cast<BigElement>(point),
                               static_cast<BigElement>(getMax())));
}

IntervalI IntervalI::expandedTo(IntervalI const& other) const {
    if (other.isEmpty()) {
        return *this;
    }
    if (this->isEmpty()) {
        return other;
    }
    return fromMinMax(std::min(static_cast<BigElement>(other.getMin()),
                               static_cast<BigElement>(getMin())),
                      std::max(static_cast<BigElement>(other.getMax()),
                               static_cast<BigElement>(getMax())));
}

IntervalI & IntervalI::clipTo(IntervalI const& other) noexcept {
    if (isEmpty()) return *this;
    if (other.isEmpty()) {
        *this = IntervalI();
        return *this;
    }
    Element max = getMax();
    Element otherMin = other.getMin();
    Element otherMax = other.getMax();
    if (otherMin > _min) {
        _min = otherMin;
    }
    if (otherMax < max) {
        max = otherMax;
    }
    if (max < _min) {
        *this = IntervalI();
        return *this;
    }
    _size = 1 + max - _min;
    return *this;
}

bool IntervalI::operator==(IntervalI const& other) const noexcept {
    return other._min == this->_min && other._size == this->_size;
}

bool IntervalI::operator!=(IntervalI const& other) const noexcept {
    return !(other == *this);
}

std::size_t IntervalI::hash_value() const noexcept {
    // Completely arbitrary seed
    return utils::hashCombine(17, _min, _size);
}

IntervalI::Element IntervalI::checkForOverflow(BigElement x, char const * where) {
    if (x < std::numeric_limits<Element>::min() || x > std::numeric_limits<Element>::max()) {
        throw LSST_EXCEPT(
            pex::exceptions::OverflowError,
            (boost::format("Integer overflow (%d) in interval %s.") % x % where).str()
        );
    }
    return static_cast<Element>(x);
}

IntervalI::IntervalI(BigElement min, BigElement size) :
    _min(checkForOverflow(min, "minimum")), _size(checkForOverflow(size, "size"))
{
    checkForOverflow(min + size - 1, "maximum");
}


IntervalD::Element const IntervalD::EPSILON = std::numeric_limits<Element>::epsilon() * 2;

IntervalD::Element const IntervalD::INVALID = std::numeric_limits<Element>::quiet_NaN();

IntervalD::IntervalD() noexcept : _min(INVALID), _max(INVALID) {}

IntervalD IntervalD::fromMinMax(Element min, Element max) noexcept {
    if (max >= min) {
        return IntervalD(min, max);
    } else {
        return IntervalD();
    }
}

IntervalD IntervalD::fromMinSize(Element min, Element size) noexcept {
    if (size >= 0 && min < std::numeric_limits<Element>::infinity()) {
        return IntervalD(min, min + size);
    } else {
        return IntervalD();
    }
}

IntervalD IntervalD::fromMaxSize(Element max, Element size) noexcept {
    if (size >= 0 && max > -std::numeric_limits<Element>::infinity()) {
        return IntervalD(max - size, max);
    } else {
        return IntervalD();
    }
}

IntervalD IntervalD::fromCenterSize(Element center, Element size) noexcept {
    Element min = center - 0.5*size;
    return lsst::geom::IntervalD::fromMinSize(min, size);
}

IntervalD::IntervalD(IntervalI const& other) noexcept
        : _min(other.getMin() - 0.5),
          _max(other.getMax() + 0.5) {
    if (other.isEmpty()) *this = IntervalD();
}

bool IntervalD::contains(Element point) const noexcept {
    return point >= this->getMin() && point <= this->getMax();
}

bool IntervalD::contains(IntervalD const& other) const noexcept {
    return other.isEmpty() ||
           (other.getMin() >= this->getMin() && other.getMax() <= this->getMax());
}

bool IntervalD::overlaps(IntervalD const& other) const noexcept {
    return !isDisjointFrom(other);
}

bool IntervalD::isDisjointFrom(IntervalD const & other) const noexcept {
    if (isEmpty() || other.isEmpty()) {
        return true;
    }
    return getMin() > other.getMax() || getMax() < other.getMin();
}

IntervalD & IntervalD::dilateBy(Element buffer) noexcept {
    if (!isEmpty()) {
        _min -= buffer;
        _max += buffer;
        if (_min >= _max) *this = IntervalD();
    }
    return *this;
}

IntervalD & IntervalD::shiftBy(Element offset) noexcept {
    if (!isEmpty()) {
        _min += offset;
        _max += offset;
    }
    return *this;
}

IntervalD & IntervalD::reflectAbout(Element point) noexcept {
    if (!isEmpty()) {
        std::swap(_min, _max);
        _min = point - _min;
        _max = point - _max;
        // _size should remain unchanged
    }
    return *this;
}

IntervalD & IntervalD::expandTo(Element point) noexcept {
    if (isEmpty()) {
        _min = point;
        _max = point;
        return *this;
    }
    if (point < _min) {
        _min = point;
    } else if (point >= _max) {
        _max = point;
    }
    return *this;
}

IntervalD & IntervalD::expandTo(IntervalD const& other) noexcept {
    if (other.isEmpty()) return *this;
    if (this->isEmpty()) {
        *this = other;
        return *this;
    }
    Element otherMin = other.getMin();
    Element otherMax = other.getMax();
    if (otherMin < _min) {
        _min = otherMin;
    }
    if (otherMax > _max) {
        _max = otherMax;
    }
    return *this;
}

IntervalD & IntervalD::clipTo(IntervalD const& other) noexcept {
    if (isEmpty()) return *this;
    if (other.isEmpty()) {
        *this = IntervalD();
        return *this;
    }
    Element otherMin = other.getMin();
    Element otherMax = other.getMax();
    if (otherMin > _min) {
        _min = otherMin;
    }
    if (otherMax < _max) {
        _max = otherMax;
    }
    if (_max < _min) {
        *this = IntervalD();
    }
    return *this;
}

bool IntervalD::operator==(IntervalD const& other) const noexcept {
    return (other.isEmpty() && this->isEmpty()) ||
           (other._min == this->_min && other._max == this->_max);
}

bool IntervalD::operator!=(IntervalD const& other) const noexcept {
    return !(other == *this);
}

std::size_t IntervalD::hash_value() const noexcept {
    if (isEmpty()) {
        // All empty intervals are equal and must have equal hashes
        return 179;
    } else {
        // Completely arbitrary seed
        return utils::hashCombine(17, _min, _max);
    }
}

std::ostream& operator<<(std::ostream& os, IntervalI const& interval) {
    if (interval.isEmpty()) return os << "IntervalI()";
    return os << "IntervalI(" << interval.getMin() << ", " << interval.getSize() << ")";
}

std::ostream& operator<<(std::ostream& os, IntervalD const& interval) {
    if (interval.isEmpty()) return os << "IntervalD()";
    return os << "IntervalD(" << interval.getMin() << ", " << interval.getSize() << ")";
}

}  // namespace geom
}  // namespace lsst
