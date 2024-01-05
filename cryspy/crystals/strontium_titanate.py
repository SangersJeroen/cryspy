from numpy import array, sin, cos, pi, stack
from ..structures.atom import Atom
from ..structures.unit_cell import UnitCell

A = B = C = 3.91

STRONTIUM = Atom("Strontium", "Sr", 38, 0, 0, 0, 38)
TITANIUM = Atom("Titanium", "Ti", 22, 1 / 2 * A, 1 / 2 * B, 1 / 2 * C, 22)
OXYGEN0 = Atom("Oxygen", "O", 8, 1 / 2 * A, 1 / 2 * B, 0 * C, 3)
OXYGEN1 = Atom("Oxygen", "O", 8, 0 * A, 1 / 2 * B, 1 / 2 * C, 3)
OXYGEN2 = Atom("Oxygen", "O", 8, 1 / 2 * A, 0 * B, 1 / 2 * C, 3)

UNIT_CELL: UnitCell = UnitCell([STRONTIUM, TITANIUM, OXYGEN0, OXYGEN1, OXYGEN2])

CELL_VECTOR_A = array([A, 0, 0]) / 2
CELL_VECTOR_B = array([0, B, 0]) / 2
CELL_VECTOR_C = array([0, 0, C]) / 2
LATTICE_VECTORS = stack([CELL_VECTOR_A, CELL_VECTOR_B, CELL_VECTOR_C])
