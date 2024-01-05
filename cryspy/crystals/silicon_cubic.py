from numpy import array, sin, cos, pi, stack
from ..structures.atom import Atom
from ..structures.unit_cell import UnitCell

A = B = C = 5.44

SILICON0 = Atom("Silicon", "Si", 14, 0, 0, 0)
SILICON2 = Atom("Silicon", "Si", 14, 1 / 4 * A, 1 / 4 * B, 1 / 4 * C)
SILICON1 = Atom("Silicon", "Si", 14, 1 / 2 * A, 0 * B, 1 / 2 * C)
SILICON3 = Atom("Silicon", "Si", 14, 0, 1 / 2 * B, 1 / 2 * C)
SILICON4 = Atom("Silicon", "Si", 14, 1 / 2 * A, 1 / 2 * B, 0 * C)
SILICON5 = Atom("Silicon", "Si", 14, 1 / 4 * A, 3 / 4 * B, 3 / 4 * C)
SILICON6 = Atom("Silicon", "Si", 14, 3 / 4 * A, 1 / 4 * B, 3 / 4 * C)
SILICON7 = Atom("Silicon", "Si", 14, 3 / 4 * A, 3 / 4 * B, 1 / 4 * C)

UNIT_CELL: UnitCell = UnitCell(
    [SILICON0, SILICON1, SILICON2, SILICON3, SILICON4, SILICON5, SILICON6, SILICON7]
)

CELL_VECTOR_A = array([A, 0, 0]) / 2
CELL_VECTOR_B = array([0, B, 0]) / 2
CELL_VECTOR_C = array([0, 0, C]) / 2
LATTICE_VECTORS = stack([CELL_VECTOR_A, CELL_VECTOR_B, CELL_VECTOR_C])
