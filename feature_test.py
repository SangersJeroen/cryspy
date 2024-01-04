import numpy as np

from cryspy.crystals.silicon_cubic import (
    SILICON_UNIT_CELL as unit_cell,
    LATTICE_VECTORS as vectors,
    A as A,
    B as B,
    C as C,
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
    colors = {"W": (1, 0, 0, 1), "Se": (0, 1, 0, 1), "Si": (1, 0, 1, 1)}
    reps = 10

    lattice1 = CrystalLattice(unit_cell)
    lattice1.construct_lattice(vectors, (reps, reps, reps))

    figure = Figure((24, 12))
    axis = figure.add_axes(
        aspect=1,
        rect=[0, 0, 1, 1],
        # interface=PanZoom(name="panzoom", aspect=1),
        interface=Trackball(name="trackball"),
        zscale=LinearScale(domain=[0.0, 1.0], range=[-1 / reps / C, 1 / reps / C]),
        xscale=LinearScale(domain=[0.0, 1.0], range=[-1 / reps / A, 1 / reps / A]),
        yscale=LinearScale(domain=[0.0, 1.0], range=[-1 / reps / B, 1 / reps / B]),
    )
    collection = PointCollection("agg", color="local", size="local")
    axis.add_drawable(collection)

    specimen = Specimen([lattice1])
    specimen.project_onto_plane(np.array([1, 2, 1]))
    lattice1.plot_lattice_3d(collection, colors)
    figure.show()
