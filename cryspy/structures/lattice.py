from .unit_cell import UnitCell
from numpy._typing import NDArray
from numpy import array, prod, dot, sin, cos, arange, meshgrid, double
from copy import deepcopy
from glumpy.graphics.collections.agg_point_collection import AggPointCollection


class CrystalLattice:
    def __init__(self, unit_cell: UnitCell) -> None:
        self.unit_cell = unit_cell
        self.origin = None

    def construct_lattice(
        self, vectors: NDArray[double], repetitions: tuple[int, int, int]
    ) -> None:
        self.lattice = [deepcopy(self.unit_cell) for _ in range(prod(repetitions))]

        if vectors.shape[0] != len(repetitions):
            raise ValueError("Number of Vectors and repetitions must have same length")

        xreps, yreps, zreps = repetitions
        xrep_range = arange(xreps) - (xreps - 1) / 2
        yrep_range = arange(yreps) - (yreps - 1) / 2
        zrep_range = arange(zreps) - (zreps - 1) / 2
        xidx, yidx, zidx = meshgrid(xrep_range, yrep_range, zrep_range, indexing="ij")

        xidx = xidx.flatten()
        yidx = yidx.flatten()
        zidx = zidx.flatten()

        for index, cell in enumerate(self.lattice):
            if self.origin is not None:
                origin = self.origin
            else:
                origin = array([0, 0, 0])
            position = (
                dot(vectors.T, array([xidx[index], yidx[index], zidx[index]])) + origin
            )
            cell.set_origin(*position)
            cell.place_atoms()

    def set_origin(self, vector: NDArray[double]) -> None:
        self.origin = vector

    def rotate_lattice(self, angle: float) -> None:
        self.angle = angle
        rotation_matrix = array(
            [
                [cos(angle), -sin(angle), 0],
                [sin(angle), cos(angle), 0],
                [0, 0, 1],
            ]
        )
        for cell in self.lattice:
            cell.origin = rotation_matrix.dot(cell.origin)
            cell.angle = angle
            cell.place_atoms()

    def plot_lattice(
        self,
        collection: AggPointCollection,
        colors: dict[str, tuple[float, float, float, float]],
    ) -> None:
        for cell in self.lattice:
            for atom in cell:
                pos = atom.get_position()
                collection.append(
                    (*pos[:-1], 0), color=colors[atom.symbol], s=atom.size
                )

    def plot_lattice_3d(
        self,
        collection: AggPointCollection,
        colors: dict[str, tuple[float, float, float, float]],
    ) -> None:
        for cell in self.lattice:
            for atom in cell:
                pos = atom.get_position()
                print(pos)
                collection.append(pos, color=colors[atom.symbol], s=atom.size)
