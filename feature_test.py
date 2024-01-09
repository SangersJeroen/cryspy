#!/usr/bin/env python3
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
    lattice = CrystalLattice(unit_cell)
    lattice.construct_lattice(vectors, (50, 50, 1))

    ranges: list[tuple[float, float]] = [(-10, 10), (-10, 10), (-10, 10)]

    latticeA = CrystalLattice(unit_cell)
    latticeA.set_origin(np.array([0, 0, C]))
    latticeA.construct_lattice(vectors, (50, 50, 1))
    latticeA.rotate_lattice(20)
    specimen = Specimen([lattice, latticeA])
    # specimen.project_onto_plane((0, 0, 1))
    specimen.build_model()
    specimen.plot_3d()
    specimen.plot_2d()
    arrs = specimen.point_mass_arrays()
    specimen.export_point_mass_dict(
        "wse",
        extent=ranges,
    )
