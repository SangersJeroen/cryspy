import numpy as np

from cryspy.crystals.tungstendiselenide import (
    TUNGSTEN_DISELENIDE_UNIT_CELL as wse2_unit_cell,
    LATTICE_VECTORS as wse2_vectors,
)

from cryspy.structures.lattice import CrystalLattice
from glumpy.api.matplotlib import Figure, PanZoom, LinearScale, PointCollection, Trackball, OrthographicProjection, IdentityProjection

if __name__ == "__main__":
    colors = {
        "W": (1, 0, 0, 1),
        "Se": (0, 1, 0, 1),
    }
    offset = np.array([np.random.rand(), np.random.rand(), 12.98 * 2])
    reps = 100

    lattice1 = CrystalLattice(wse2_unit_cell)
    lattice1.construct_lattice(wse2_vectors, (reps, reps, 2))
    lattice1.rotate_lattice(123.7 / 180 * np.pi)

    lattice2 = CrystalLattice(wse2_unit_cell)
    lattice2.set_origin(offset)
    lattice2.construct_lattice(wse2_vectors, (reps, reps, 2))
    lattice2.rotate_lattice(125.6 / 180 * np.pi)
    # lattice2.rotate_lattice(50 / 180 * np.pi)

    lattice3 = CrystalLattice(wse2_unit_cell)
    lattice3.set_origin(offset * 2.0)
    lattice3.construct_lattice(wse2_vectors, (reps, reps, 2))
    lattice3.rotate_lattice(124.4 / 180 * np.pi)


    figure = Figure((24, 12))
    axis = figure.add_axes(
        aspect=1,
        rect=[0, 0, 1, 1],
        interface=PanZoom(name="panzoom", aspect=1),
        # interface=Trackball(name="trackball"),
        # zscale=LinearScale(domain=[0.0, 1.0], range=[-1 / reps / 12.98, 1 / reps / 12.98]),
        xscale=LinearScale(
            domain=[0.0, 1.0], range=[-1 / reps / 3.32, 1 / reps / 3.32]
        ),
        yscale=LinearScale(
            domain=[0.0, 1.0], range=[-1 / reps / 3.32, 1 / reps / 3.32]
        ),
    )
    collection = PointCollection("agg", color="local", size="local")
    axis.add_drawable(collection)

    # lattice1.plot_lattice_3d(collection, colors)
    # lattice2.plot_lattice_3d(collection, colors)
    # lattice3.plot_lattice_3d(collection, colors)
    lattice1.plot_lattice(collection, colors)
    lattice2.plot_lattice(collection, colors)
    lattice3.plot_lattice(collection, colors)
    figure.show()
