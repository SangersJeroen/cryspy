#!/usr/bin/env python3
import numpy as np

from cryspy.crystals.hexagonal_boronnitride import (
    UNIT_CELL as unit_cell,
    LATTICE_VECTORS as vectors,
)

from cryspy.structures.lattice import CrystalLattice
from cryspy.structures.specimen import Specimen


if __name__ == "__main__":
    latticeA = CrystalLattice(unit_cell)
    latticeA.construct_lattice(vectors, (100, 100, 1))

    specimen = Specimen([latticeA])
    specimen.build_model()
    specimen.plot_3d()
    # specimen.project_onto_plane(np.array([0, 1, 1]))
    specimen.plot_2d(view_axis=np.array([1, 0, 0]))

    specimen.export_to_xyzfile_temsim(r'output/hbn')
