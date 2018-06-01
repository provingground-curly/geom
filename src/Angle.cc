
#include <iostream>

#include "lsst/geom/Angle.h"

namespace lsst {
namespace geom {

std::ostream& operator<<(std::ostream& s, Angle a) { return s << static_cast<double>(a) << " rad"; }

}
}
