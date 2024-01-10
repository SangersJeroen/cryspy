#!/usr/bin/env python3
import numpy as np

from cryspy.crystals.tungstendiselenide import (
    A as a,
    B as b,
    C as c,
    UNIT_CELL as unit_cell,
    LATTICE_VECTORS as vectors,
)

from cryspy.structures.lattice import CrystalLattice
from cryspy.structures.specimen import Specimen

ranges = [(0, 80), (0, 80), (0, 80)]
reps = 110


if __name__ == "__main__":
    latticeA = CrystalLattice(unit_cell)
    latticeA.set_origin(np.array([1.27, 0.13, 2 * c]))
    latticeA.construct_lattice(vectors, (reps, reps, 2))
    latticeA.rotate_lattice(123.7 / 180 * np.pi)

    latticeB = CrystalLattice(unit_cell)
    latticeB.set_origin(np.array([0.91, 0.51, 4 * c]))
    latticeB.construct_lattice(vectors, (reps, reps, 2))
    latticeB.rotate_lattice(125.6 / 180 * np.pi)

    latticeC = CrystalLattice(unit_cell)
    latticeC.set_origin(np.array([0.837981, 0.71, 6 * c]))
    latticeC.construct_lattice(vectors, (reps, reps, 2))
    latticeC.rotate_lattice(124.4 / 180 * np.pi)

    specimen = Specimen([latticeA, latticeB, latticeC])
    specimen.build_model()
    # specimen.plot_3d()
    specimen.plot_2d(ranges=ranges)
    specimen.export_to_xyzfile("output/e2_heterostructure", extent=ranges)
