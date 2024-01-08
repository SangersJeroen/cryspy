from cryspy.crystals.tungstendiselenide import (
    UNIT_CELL as unit_cell,
    LATTICE_VECTORS as vectors,
    A as A,
    B as B,
    C as C,
)

from cryspy.structures.lattice import CrystalLattice
from cryspy.structures.specimen import Specimen

if __name__ == "__main__":
    lattice1 = CrystalLattice(unit_cell)
    lattice1.construct_lattice(vectors, (40, 40, 5))

    specimen = Specimen([lattice1])
    # specimen.export_point_mass_dict("test.dat")
    # specimen.project_onto_plane((0, 0, 1))
    specimen.plot_3d()
