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
    lattice.set_origin(np.array([0, 0, 2 * C]))
    lattice.construct_lattice(vectors, (50, 50, 1))

    ranges: list[tuple[float, float]] = [(0, 20), (0, 20), (0, 20)]

    latticeA = CrystalLattice(unit_cell)
    latticeA.set_origin(np.array([0, 0.5 * B, C]))
    latticeA.construct_lattice(vectors, (50, 50, 1))
    latticeA.rotate_lattice(20)
    specimen = Specimen([lattice, latticeA])
    # specimen.project_onto_plane((0, 0, 1))
    specimen.build_model()
    specimen.plot_3d()
    specimen.plot_2d(ranges=ranges)
    arrs = specimen.point_mass_arrays()
    specimen.export_to_xyzfile(
        "output/wse",
        extent=ranges,
    )
