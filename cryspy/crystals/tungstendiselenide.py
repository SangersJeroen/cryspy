from numpy import array, sin, cos, pi, stack
from ..structures.atom import Atom
from ..structures.unit_cell import UnitCell

A = 3.32
B = 3.32
C = 12.98

CELL_VECTOR_A = array([sin(30 / 180 * pi), cos(30 / 180 * pi), 0]) * 3.32
CELL_VECTOR_B = array([1, 0, 0]) * 3.32
CELL_VECTOR_C = array([0, 0, 12.98])

CELL_DIAG = CELL_VECTOR_A + CELL_VECTOR_B

LATTICE_VECTORS = stack([CELL_VECTOR_A, CELL_VECTOR_B, CELL_VECTOR_C])

TUNGSTEN_ATOM: Atom = Atom("Tungsten", "W", 74, *(2 / 3 * CELL_DIAG[:-1]), 0, 2)
SELENIUM_ATOM0: Atom = Atom("Selenium", "Se", 34, *(1 / 3 * CELL_DIAG[:-1]), 1.67, 2)
SELENIUM_ATOM1: Atom = Atom("Selenium", "Se", 34, *(1 / 3 * CELL_DIAG[:-1]), -1.67, 2)

UNIT_CELL: UnitCell = UnitCell([TUNGSTEN_ATOM, SELENIUM_ATOM0, SELENIUM_ATOM1])
