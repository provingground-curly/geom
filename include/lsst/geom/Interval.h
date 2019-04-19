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

#ifndef LSST_GEOM_INTERVAL_H
#define LSST_GEOM_INTERVAL_H

#include <vector>

#include "boost/format.hpp"
#include "ndarray.h"

#include "lsst/geom/Point.h"
#include "lsst/geom/Extent.h"

namespace lsst {
namespace geom {

class IntervalD;

/**
 *  A 1-d integer coordinate range.
 *
 *  IntervalI is an inclusive range that represents region of pixels.  An
 *  IntervalI never has negative size; the empty interval is defined to
 *  have zero-size size, and is treated as though it does not have a
 *  well-defined position (regardless of the return value of getMin() or
 *  getMax() for an empty interval).
 *
 *  All mutating IntervalI methods have a counterpart that returns a new
 *  object instead, and only the latter are exposed to Python.  That
 *  (intentionally) makes mutation in Python impossible; while this is not
 *  quite the same as making it fully immutable (as it is still possible to
 *  pass a Python IntervalI by non-const reference to a different C++ function
 *  that modifies it), the distinction is sufficiently rare that it can almost
 *  always be ignored.  IntervalI cannot be fully immutable (i.e. in C++)
 *  without making it much less convenient to use as a value type.
 *
 *  All IntervalI methods that mutate self or return a new instance (and are
 *  not marked `noexcept`) throw OverflowError if either bound or the size
 *  would be too large to fit in `int`.
 *
 *  @internal
 *
 *  IntervalI internally stores its minimum point and size, because we expect
 *  these will be the most commonly accessed quantities.
 *
 *  IntervalI sets the minimum point to the origin for an empty interval, and
 *  returns -1 for both elements of the maximum point in that case.
 */
class IntervalI final {
public:

    using Element = int;

    enum class EdgeHandlingEnum { EXPAND, SHRINK };

    /// Construct an empty interval.
    IntervalI() noexcept : _min(0), _size(0) {}

    //@{
    /**
     *  Construct an interval that contains all of the given points.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     */
    template <typename Iter>
    static IntervalI fromHull(Iter first, Iter last) {
        IntervalI result;
        for (auto i = first; i != last; ++i) {
            result.expandTo(*i);
        }
        return result;
    }
    static IntervalI fromHull(std::vector<Element> const & elements) {
        return fromHull(elements.begin(), elements.end());
    }
    static IntervalI fromHull(ndarray::Array<Element const, 1> const & elements) {
        return fromHull(elements.begin(), elements.end());
    }
    //@}

    /**
     *  Construct an interval from its lower and upper bounds.
     *
     *  If min > max, the resulting interval is empty.
     *
     *  @param[in] min       Minimum coordinate (inclusive).
     *  @param[in] max       Maximum coordinate (inclusive).
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     */
    static IntervalI fromMinMax(Element min, Element max);

    /**
     *  Construct an interval from its lower bound and size.
     *
     *  @param[in] min       Minimum coordinate (inclusive).
     *  @param[in] size      Number of pixels in interval.
     *                       Must be nonnegative.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     * @throws lsst::pex::exceptions::InvalidParameterError Thrown if size < 0.
     */
    static IntervalI fromMinSize(Element min, Element size);

    /**
     *  Construct an interval from its upper bound and size.
     *
     *  @param[in] max       Maximum coordinate (inclusive).
     *  @param[in] size      Number of pixels in interval.
     *                       Must be nonnegative.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     * @throws lsst::pex::exceptions::InvalidParameterError Thrown if size < 0.
     */
    static IntervalI fromMaxSize(Element max, Element size);

    /**
     * Create an interval centered as closely as possible on a particular
     * point.
     *
     * @param center The desired center of the interval.
     * @param size   Number of pixels in interval.  Must be nonnegative.
     *
     * @returns if `size` is positive, an interval with size `size`; if `size`
     *          is zero, an empty interval. If the returned interval is not
     *          empty, its center shall be within half a pixel of `center` in
     *          either dimension.
     *
     * @throws lsst::pex::exceptions::OverflowError Thrown if the resulting
     *     interval would overflow.
     * @throws lsst::pex::exceptions::InvalidParameterError Thrown if `center`
     *     is not finite and size is nonzero, or if size is negative.
     */
    static IntervalI fromCenterSize(double center, Element size);

    /**
     *  Construct an integer interval from a floating-point interval.
     *
     *  Floating-point to integer interval conversion is based on the concept
     *  that a pixel is not an infinitesimal point but rather a square of unit
     *  size centered on integer-valued coordinates.  Converting a
     *  floating-point interval to an integer interval thus requires a choice
     *  on how to handle pixels which are only partially contained by the
     *  input floating-point interval.
     *
     *  @param[in] other          A floating-point interval to convert.
     *  @param[in] edgeHandling   If EXPAND, the integer interval will contain
     *                            any pixels that overlap the floating-point
     *                            interval.  If SHRINK, the integer interval
     *                            will contain only pixels completely
     *                            contained by the floating-point interval.
     *
     * @throws lsst::pex::exceptions::InvalidParameterError Thrown if `other`
     *     is not finite.
     */
    explicit IntervalI(IntervalD const& other, EdgeHandlingEnum edgeHandling = EdgeHandlingEnum::EXPAND);

    /// Standard copy constructor.
    IntervalI(IntervalI const&) noexcept = default;

    /// Standard move constructor.
    IntervalI(IntervalI&&) noexcept = default;

    ~IntervalI() noexcept = default;

    void swap(IntervalI& other) noexcept {
        std::swap(_min, other._min);
        std::swap(_size, other._size);
    }

    /// Standard copy assignment operator.
    IntervalI& operator=(IntervalI const&) noexcept = default;

    /// Standard move assignment operatior.
    IntervalI& operator=(IntervalI&&) noexcept = default;

    /**
     *  @name Min/Max Accessors
     *
     *  Return the minimum and maximum coordinates of the interval (inclusive).
     */
    //@{
    Element getMin() const noexcept { return _min; }
    Element getMax() const noexcept { return _min + _size - 1; }
    //@}

    /**
     *  @name Begin/End Accessors
     *
     *  Return STL-style begin (inclusive) and end (exclusive) coordinates for
     *  the interval.
     */
    //@{
    Element getBegin() const noexcept { return _min; }
    Element getEnd() const noexcept { return _min + _size; }
    //@}

    /**
     *  @name Size Accessors
     *
     *  Return the size of the interval in pixels.
     */
    //@{
    Element getSize() const noexcept { return _size; }
    //@}

    /// Return slice to extract the interval's region from an ndarray::Array.
    ndarray::View<boost::fusion::vector1<ndarray::index::Range> > getSlice() const;

    /// Return true if the interval contains no points.
    bool isEmpty() const noexcept { return _size == 0; }

    /// Return true if the interval contains the point.
    bool contains(Element point) const noexcept;

    /**
     *  Return true if all points contained by other are also contained by
     *  this.
     *
     *  An empty interval is contained by every other interval, including
     *  other empty intervals.
     */
    bool contains(IntervalI const& other) const noexcept;

    //@{
    /**
     *  Return true if any points in other are also in this.
     *
     *  Any overlap operation involving an empty interval returns false.
     */
    bool overlaps(IntervalI const& other) const noexcept;
    bool intersects(IntervalI const& other) const noexcept { return overlaps(other); }
    //@}

    /**
     *  Return true if there a no points in both `this` and `other`.
     */
    bool isDisjointFrom(IntervalI const & other) const noexcept;

    //@{
    /**
     *  Increase the size of the interval by the given amount in both
     *  directions.
     *
     *  If `buffer` is negative, this is equivalent to eroding by `-buffer`.
     *
     *  Empty intervals remain empty after dilation.
     */
    IntervalI dilatedBy(Element buffer) const;
    IntervalI & dilateBy(Element buffer) {
        return (*this = dilatedBy(buffer)); // delegate mutator to factory for exception safety.
    }
    //@}

    //@{
    /**
     *  Decrease the size of the interval by the given amount in both
     *  directions.
     *
     *  If `buffer` is negative, this is equivalent to dilating by `-buffer`.
     *
     *  If the final size of the interval is less than or equal to zero (which
     *  can happen if `buffer` is negative), the interval will be made empty.
     *
     *  Empty intervals remain empty after erosion.
     */
    IntervalI erodedBy(Element buffer) const { return dilatedBy(-buffer); }
    IntervalI & erodeBy(Element buffer) { return dilateBy(-buffer); }
    //@}

    //@{
    /**
     * Shift the position of the interval by the given offset.
     *
     * Empty intervals remain empty when shifted.
     */
    IntervalI shiftedBy(Element offset) const;
    IntervalI & shiftBy(Element offset) {
        return (*this = shiftedBy(offset)); // delegate mutator to factory for exception safety.
    }
    //@}

    //@{
    /**
     * Reflect an interval about a point.
     *
     * Empty intervals remain empty when reflected.
     */
    IntervalI reflectedAbout(Element point) const;
    IntervalI & reflectAbout(Element point) {
        return (*this = reflectedAbout(point)); // delegate mutator to factory for exception safety.
    }
    //@}

    //@{
    /**
     * Expand an interval to ensure that `contains(other)` is true.
     *
     * Expanding an empty interval with a single point yields an interval with
     * size=1 at that point; expanding an empty interval with a second interval
     * is equivalent to assignment.
     */
    IntervalI expandedTo(Element other) const;
    IntervalI expandedTo(IntervalI const & other) const;
    IntervalI & expandTo(Element other) {
        return (*this = expandedTo(other)); // delegate mutator to factory for exception safety.
    }
    IntervalI & expandTo(IntervalI const & other) {
        return (*this = expandedTo(other)); // delegate mutator to factory for exception safety.
    }
    //@}

    //@{
    /**
     * Shrink an interval to ensure that it is contained by other.
     *
     * In particular, if `other` and this interval do not overlap this
     * interval will become empty.
     */
    IntervalI & clipTo(IntervalI const& other) noexcept;
    IntervalI clippedTo(IntervalI const & other) const noexcept { return IntervalI(*this).clipTo(other); }
    //@}

    /**
     *  Compare two intervals for equality.
     *
     *  All empty intervals are equal.
     */
    bool operator==(IntervalI const& other) const noexcept;

    /**
     *  Compare two intervals for equality.
     *
     *  All empty intervals are equal.
     */
    bool operator!=(IntervalI const& other) const noexcept;

    /// Return a hash of this object.
    std::size_t hash_value() const noexcept;

    std::string toString() const {
        return (boost::format("IntervalI(min=%s, max=%s)") % getMin() % getMax()).str();
    }

private:

    using BigElement = long long;

    static Element checkForOverflow(BigElement x, char const * where);

    // Internal constructor taking a larger interger size so it can check for
    // overflow.
    IntervalI(BigElement min, BigElement size);

    Element _min;
    Element _size;
};

/**
 *  A floating-point coordinate rectangle geometry.
 *
 *  IntervalD is closed (its bounds are considered included in the interval).
 *  An interval never has negative size; the empty interval is defined to
 *  zero-size size and its minimum and maximum values set to NaN.  Non-empty
 *  intervals representing infinitesimal points may also have zero size, and
 *  but not considered empty.
 *
 *  @internal
 *
 *  IntervalD internally stores its minimum point and maximum point, instead
 *  of minimum point and size, to ensure roundoff error does not affect
 *  whether points are contained by the interval.
 *
 *  Despite some recommendations to the contrary, IntervalD sets the minimum
 *  and maximum points to NaN for an empty interval.  In almost every case,
 *  special checks for emptiness would have been necessary anyhow, so there
 *  was little to gain in using the minimum > maximum condition to denote an
 *  empty interval, as was used in IntervalI.
 */
class IntervalD final {
public:
    using Element = double;

    /**
     *  Value the maximum coordinate is multiplied by to increase it by the
     *  smallest possible amount.
     */
    static Element const EPSILON;

    /// Value used to specify undefined coordinate values.
    static Element const INVALID;

    /// Construct an empty interval.
    IntervalD() noexcept;

    //@{
    /**
     *  Construct an interval that contains all of the given points.
     */
    template <typename Iter>
    static IntervalD fromHull(Iter first, Iter last) {
        IntervalD result;
        for (auto i = first; i != last; ++i) {
            result.expandTo(*i);
        }
        return result;
    }
    static IntervalD fromHull(std::vector<Element> const & elements) {
        return fromHull(elements.begin(), elements.end());
    }
    static IntervalD fromHull(ndarray::Array<Element const, 1> const & elements) {
        return fromHull(elements.begin(), elements.end());
    }
    //@}

    /**
     *  Construct an interval from its lower and upper bounds.
     *
     *  @param[in] min       Minimum coordinate (inclusive).
     *  @param[in] max       Maximum coordinate (inclusive).
     *
     *  If `min > max` or either is NaN, an empty interval is returned.
     */
    static IntervalD fromMinMax(Element min, Element max) noexcept;

    /**
     *  Construct an interval from its lower bound and size.
     *
     *  @param[in] min       Minimum coordinate (inclusive).
     *  @param[in] size      Size of interval.
     *
     *  If `size < )`, `size` is NaN, or `min` is infinty or NaN, an empty
     *  interval is returned.
     */
    static IntervalD fromMinSize(Element min, Element size) noexcept;

    /**
     *  Construct an interval from its upper bound and size.
     *
     *  @param[in] max       Maximum coordinate (inclusive).
     *  @param[in] size      Size of interval.
     *
     *  If `size < 0`, `size` is NaN, or `max` is -infinity or NaN, an empty
     *  interval is returned.
     */
    static IntervalD fromMaxSize(Element max, Element size) noexcept;

    /**
     * Create an interval centered as closely as possible on a particular
     * point.
     *
     * @param center The desired center of the interval.
     * @param size   Number of pixels in interval.
     *
     * @returns if `size` is nonnegative, an interval with size `size`;
     *          otherwise, an empty interval. If the returned interval is not
     *          empty, it shall be centered on `center`.
     */
    static IntervalD fromCenterSize(double center, Element size) noexcept;

    /**
     *  Construct a floating-point interval from an integer interval.
     *
     *  Integer to floating-point interval conversion is based on the concept
     *  that a pixel is not an infinitesimal point but rather a square of unit
     *  size centered on integer-valued coordinates.  While the output
     *  floating-point interval thus has the same size as the input integer
     *  interval, its minimum/maximum coordinates are 0.5 smaller/greater.
     */
    explicit IntervalD(IntervalI const& other) noexcept;

    /// Standard copy constructor.
    IntervalD(IntervalD const&) noexcept = default;

    /// Standard move constructor.
    IntervalD(IntervalD&&) noexcept = default;

    ~IntervalD() noexcept = default;

    void swap(IntervalD& other) noexcept {
        std::swap(_min, other._min);
        std::swap(_max, other._max);
    }

    /// Standard copy assignment operator.
    IntervalD& operator=(IntervalD const&) noexcept = default;

    /// Standard move assignment operator.
    IntervalD& operator=(IntervalD&&) noexcept = default;

    /**
     *  @name Min/Max Accessors
     *
     *  Return the minimum (inclusive) and maximum (exclusive) coordinates of
     *  the interval.
     */
    //@{
    Element getMin() const noexcept { return _min; }
    Element getMax() const noexcept { return _max; }
    //@}

    /**
     *  @name Size Accessors
     *
     *  Return the size of the interval.
     */
    //@{
    Element getSize() const noexcept { return isEmpty() ? 0.0 : _max - _min; }
    //@}

    /**
     *  @name Center Accessors
     *
     *  Return the center coordinate of the interval.
     */
    //@{
    Element getCenter() const noexcept { return 0.5*(_min + _max); }
    //@}

    /// Return true if the interval contains no points.
    bool isEmpty() const noexcept { return _min != _min; }

    /// Return true if the interval contains the point.
    bool contains(Element point) const noexcept;

    /**
     *  Return true if all points contained by other are also contained by
     *  this.
     *
     *  An empty interval is contained by every other interval, including
     *  other empty intervals.
     */
    bool contains(IntervalD const& other) const noexcept;

    //@{
    /**
     *  Return true if any points in other are also in this.
     *
     *  Any overlap operation involving an empty interval returns false.
     */
    bool overlaps(IntervalD const& other) const noexcept;
    bool intersects(IntervalD const& other) const noexcept { return overlaps(other); }
    //@}

    /**
     *  Return true if there a no points in both `this` and `other`.
     */
    bool isDisjointFrom(IntervalD const & other) const noexcept;

    //@{
    /**
     *  Increase the size of the interval by the given amount in both
     *  directions.
     *
     *  If `buffer` is negative, this is equivalent to eroding by `-buffer`.
     *
     *  Empty intervals remain empty after dilation.
     *
     *  Result is unspecified if `buffer` is NaN.
     */
    IntervalD & dilateBy(Element buffer) noexcept;
    IntervalD dilatedBy(Element buffer) const noexcept { return IntervalD(*this).dilateBy(buffer); }
    //@}

    //@{
    /**
     *  Decrease the size of the interval by the given amount in both
     *  directions.
     *
     *  If `buffer` is negative, this is equivalent to dilating by `-buffer`.
     *
     *  If the final size of the interval is less than zero (which can happen
     *  if `buffer` is negative), the interval will be made empty.
     *
     *  Empty intervals remain empty after erosion.
     *
     *  Result is unspecified if `buffer` is NaN.
     */
    IntervalD & erodeBy(Element buffer) noexcept { return dilateBy(-buffer); }
    IntervalD erodedBy(Element buffer) const noexcept { return dilatedBy(-buffer); }
    //@}

    //@{
    /**
     *  Shift the position of the interval by the given offset.
     *
     *  Empty intervals remain empty when shifted.
     *
     *  Result is specified if `offset` is not finite.
     */
    IntervalD & shiftBy(Element offset) noexcept;
    IntervalD shiftedBy(Element offset) const noexcept { return IntervalD(*this).shiftBy(offset); }
    //@}

    //@{
    /**
     *  Reflect an interval about a point.
     *
     *  Empty intervals remain empty when reflected.
     *
     *  Result is unspecified if `point` is not finite.
     */
    IntervalD & reflectAbout(Element point) noexcept;
    IntervalD reflectedAbout(Element point) const noexcept { return IntervalD(*this).reflectAbout(point); }
    //@}

    //@{
    /**
     *  Expand an interval to ensure that `contains(other)` is true.
     *
     *  Expanding an empty interval with a single point yields an interval
     *  with `size == 0` at that point; expanding an empty interval with a
     *  second interval is equivalent to assignment.
     *
     *  Expanding by NaN or an empty interval yields the original interval.
     */
    IntervalD & expandTo(Element other) noexcept;
    IntervalD & expandTo(IntervalD const & other) noexcept;
    IntervalD expandedTo(Element other) const noexcept { return IntervalD(*this).expandTo(other); }
    IntervalD expandedTo(IntervalD const & other) const noexcept { return IntervalD(*this).expandTo(other); }
    //@}

    //@{
    /**
     * Shrink an interval to ensure that it is contained by other.
     *
     * In particular, if `other` and this interval do not overlap this
     * interval will become empty.
     */
    IntervalD & clipTo(IntervalD const& other) noexcept;
    IntervalD clippedTo(IntervalD const & other) const noexcept { return IntervalD(*this).clipTo(other); }
    //@}

    /**
     *  Compare two intervals for equality.
     *
     *  All empty intervals are equal.
     */
    bool operator==(IntervalD const& other) const noexcept;

    /**
     *  Compare two intervals for equality.
     *
     *  All empty intervals are equal.
     */
    bool operator!=(IntervalD const& other) const noexcept;

    /// Return a hash of this object.
    std::size_t hash_value() const noexcept;

    std::string toString() const {
        return (boost::format("IntervalD(min=%s, max=%s)") % getMin() % getMax()).str();
    }

private:

    void _tweakMax() noexcept {
        if (_max < 0.0) {
            _max *= (1.0 - EPSILON);
        } else if (_max > 0.0) {
            _max *= (1.0 + EPSILON);
        } else {
            _max = EPSILON;
        }
    }

    IntervalD(Element min, Element max) : _min(min), _max(max) {}

    Element _min;
    Element _max;
};

std::ostream& operator<<(std::ostream& os, IntervalI const& interval);

std::ostream& operator<<(std::ostream& os, IntervalD const& interval);

inline void swap(IntervalI & a, IntervalI & b) { a.swap(b); }

inline void swap(IntervalD & a, IntervalD & b) { a.swap(b); }

}  // namespace geom
}  // namespace lsst

namespace std {
template <>
struct hash<lsst::geom::IntervalI> {
    using argument_type = lsst::geom::IntervalI;
    using result_type = size_t;
    size_t operator()(argument_type const& x) const noexcept { return x.hash_value(); }
};

template <>
struct hash<lsst::geom::IntervalD> {
    using argument_type = lsst::geom::IntervalD;
    using result_type = size_t;
    size_t operator()(argument_type const& x) const noexcept { return x.hash_value(); }
};
}  // namespace std

#endif
