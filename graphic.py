#!/usr/bin/env python3
import numpy as np
import cProfile
from time import time_ns

from cryspy.crystals.tungstendiselenide import (
    A as a,
    B as b,
    C as c,
    UNIT_CELL as unit_cell,
    LATTICE_VECTORS as vectors,
)

from cryspy.structures.lattice import CrystalLattice
from cryspy.structures.specimen import Specimen

ranges = [(-100, 100), (-100, 100), (0, 80)]

if __name__ == "__main__":
    latticeA = CrystalLattice(unit_cell)
    latticeA.set_origin(np.array([0, 0, 0]))
    latticeA.construct_lattice(vectors, (30, 30, 1))
    latticeA.rotate_lattice(0)

    latticeB = CrystalLattice(unit_cell)
    latticeB.set_origin(np.array([2*a, 3*b, c]))
    latticeB.construct_lattice(vectors, (25, 25, 1))
    latticeB.rotate_lattice(5 / 180 * np.pi)

    # latticeC = CrystalLattice(unit_cell)
    # latticeC.set_origin(np.array([5*a, 5*a, c]))
    # latticeC.construct_lattice(vectors, (25, 25, 1))
    # latticeC.rotate_lattice(2 / 180 * np.pi)

    specimen = Specimen([latticeA, latticeB])
    specimen.build_model()

    specimen.plot_2d()
    specimen.plot_3d()
    # cProfile.run("specimen.plot_2d(ranges=ranges)", sort="tottime")

    specimen.export_to_xyzfile_blender(r"./output/small_graphic")
