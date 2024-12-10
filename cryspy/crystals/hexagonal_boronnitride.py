from numpy import array, sin, cos, pi, stack

from ..structures.atom import Atom
from ..structures.unit_cell import UnitCell

A = 2.51
B = 2.51
C = 7.71

BORON = Atom("Boron", "B", 74, 1.660, 0.950, 0, 0)
NITROGREN = Atom("Nitrogen", "N", 34, 1.660, -0.980, 0, 0)

UNIT_CELL: UnitCell = UnitCell([BORON, NITROGREN])
LATTICE_VECTORS = stack(
    [
        array([sin(pi / 6), cos(pi / 6), 0]) * A / 2,
        array([-1 * sin(pi / 6), cos(pi / 6), 0]) * B / 2,
        array([0, 0, 1]) * C,
    ]
)
