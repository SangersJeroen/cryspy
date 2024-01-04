import numpy as np

from cryspy.crystals.tungstendiselenide import (
    TUNGSTEN_DISELENIDE_UNIT_CELL as wse2_unit_cell,
    LATTICE_VECTORS as wse2_vectors,
)

from cryspy.structures.lattice import CrystalLattice
from cryspy.structures.specimen import Specimen
from glumpy.api.matplotlib import (
    Figure,
    PanZoom,
    LinearScale,
    PointCollection,
    Trackball,
)

if __name__ == "__main__":
    colors = {
        "W": (1, 0, 0, 1),
        "Se": (0, 1, 0, 1),
    }
    reps = 3

    lattice1 = CrystalLattice(wse2_unit_cell)
    lattice1.construct_lattice(wse2_vectors, (reps, reps, reps))

    # figure = Figure((24, 12))
    # axis = figure.add_axes(
    #     aspect=1,
    #     rect=[0, 0, 1, 1],
    #     # interface=PanZoom(name="panzoom", aspect=1),
    #     interface=Trackball(name="trackball"),
    #     zscale=LinearScale(
    #         domain=[0.0, 1.0], range=[-1 / reps / 12.98, 1 / reps / 12.98]
    #     ),
    #     xscale=LinearScale(
    #         domain=[0.0, 1.0], range=[-1 / reps / 3.32, 1 / reps / 3.32]
    #     ),
    #     yscale=LinearScale(
    #         domain=[0.0, 1.0], range=[-1 / reps / 3.32, 1 / reps / 3.32]
    #     ),
    # )
    # collection = PointCollection("agg", color="local", size="local")
    # axis.add_drawable(collection)
    #
    # lattice1.plot_lattice_3d(collection, colors)
    # figure.show()

    figure = Figure((24, 12))
    axis = figure.add_axes(
        aspect=1,
        rect=[0, 0, 1, 1],
        # interface=PanZoom(name="panzoom", aspect=1),
        interface=Trackball(name="trackball"),
        zscale=LinearScale(
            domain=[0.0, 1.0], range=[-1 / reps / 12.98, 1 / reps / 12.98]
        ),
        xscale=LinearScale(
            domain=[0.0, 1.0], range=[-1 / reps / 3.32, 1 / reps / 3.32]
        ),
        yscale=LinearScale(
            domain=[0.0, 1.0], range=[-1 / reps / 3.32, 1 / reps / 3.32]
        ),
    )
    collection = PointCollection("agg", color="local", size="local")
    axis.add_drawable(collection)

    specimen = Specimen([lattice1])
    specimen.project_onto_plane((0, 0, 1))
    specimen._constituent_structures[0].plot_lattice_3d(collection, colors)

    figure.show()
