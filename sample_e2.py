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

sreps = 70
breps = sreps

if __name__ == "__main__":
    print(f"big layer reps: {breps}")
    print(f"small layer reps: {sreps}")

    start = time_ns()
    latticeA = CrystalLattice(unit_cell)
    latticeA.set_origin(np.array([1.27, 0.13, 0 * c]))
    latticeA.construct_lattice(vectors, (breps, breps, 1))
    latticeA.rotate_lattice(123.7 / 180 * np.pi)

    latticeB = CrystalLattice(unit_cell)
    latticeB.set_origin(np.array([0.91, 0.51, 1 * c]))
    latticeB.construct_lattice(vectors, (breps, breps, 1))
    latticeB.rotate_lattice(125.6 / 180 * np.pi)

    latticeC = CrystalLattice(unit_cell)
    latticeC.set_origin(np.array([0.1, 0.2, 2 * c]))
    latticeC.construct_lattice(vectors, (sreps, breps, 1))
    latticeC.rotate_lattice(124.4 / 180 * np.pi)
    lap0 = time_ns() - start
    start = time_ns()
    print(f"Lattices initiated in {lap0/1e9:.2f} seconds")

    specimen = Specimen([latticeA, latticeB, latticeC])
    specimen.build_model()

    lap1 = time_ns() - start
    print(f"Specimen build in {lap1/1e9:.2f} seconds")

    # specimen.plot_3d()
    specimen.plot_2d()
    # cProfile.run("specimen.plot_2d(ranges=ranges)", sort="tottime")

    specimen.export_to_xyzfile_blender("output/e2_heterostructure")
