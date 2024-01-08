import numpy as np
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
    lattice1.construct_lattice(vectors, (10, 10, 1))

    specimen = Specimen([lattice1])
    # specimen.project_onto_plane((0, 0, 1))
    specimen.build_model()
    # specimen.plot_2d(view_axis=np.array([1, 1, 1]))
    arrs = specimen.point_mass_arrays()
    print(arrs.shape)
