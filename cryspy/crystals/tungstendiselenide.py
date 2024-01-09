from numpy import array, sin, cos, pi, stack

from ..structures.atom import Atom
from ..structures.unit_cell import UnitCell

A = 3.32
B = 3.32
C = 12.98 / 2

TUNGSTEN = Atom("Tungsten", "W", 74, 1.660, 0.950, 0, 74)
SELENIUM0 = Atom("Selenium", "Se", 34, 1.660, -0.980, 1.663, 0)
SELENIUM1 = Atom("Selenium", "Se", 34, 1.660, -0.980, -1.663, 0)

UNIT_CELL: UnitCell = UnitCell([TUNGSTEN, SELENIUM0, SELENIUM1])
LATTICE_VECTORS = stack(
    [
        array([sin(pi / 6), cos(pi / 6), 0]) * A / 2,
        array([-1 * sin(pi / 6), cos(pi / 6), 0]) * B / 2,
        array([0, 0, 1]) * C,
    ]
)
