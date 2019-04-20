#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

__all__ = ["Box2I"]

import numpy as np

from lsst.utils import continueClass
from ._geom import Box2I


@continueClass  # noqa: F811
class Box2I:

    def grid(self, dtype=np.int32):
        """Return pair of arrays with the centers of all pixels in the box.

        Parameters
        ----------
        dtype : convertible to `numpy.dtype`
            The data type of the returned arrays.

        Returns
        -------
        x : `numpy.ndarray`
            Array with shape `(self.height, self.width)` containing x
            coordinate values.
        y : `numpy.ndarray`
            Array with shape `(self.height, self.width)` containing x
            coordinate values.
        """
        return np.meshgrid(self.x.arange(dtype=dtype), self.y.arange(dtype=dtype))
