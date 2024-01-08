from .unit_cell import UnitCell
from .atom import Atom
from numpy._typing import NDArray
from numpy import array, prod, dot, sin, cos, arange, meshgrid, double
from copy import deepcopy
from glumpy.graphics.collections.agg_point_collection import AggPointCollection


Vector = tuple[float, float, float] | list[float] | NDArray[double]
Position = tuple[float, float, float] | list[float] | NDArray[double]
PointMassDict = dict[int, list[Position]]


class CrystalLattice:
    def __init__(self, unit_cell: UnitCell) -> None:
        self.unit_cell = unit_cell
        self.origin = None
        self._constructed: bool = False

    def construct_lattice(
        self,
        vectors: NDArray[double],
        repetitions: tuple[int, int, int],
        construct: bool = False,
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
        if construct:
            self._construct()

    def point_mass_dict(self) -> PointMassDict:
        self._atoms: list[Atom] = []
        for unit_cell in self.lattice:
            self._atoms.extend(unit_cell.atoms)
        point_mass_positions: PointMassDict = {}
        for atom in self._atoms:
            if atom.number not in point_mass_positions:
                point_mass_positions[atom.number] = []
            point_mass_positions[atom.number].append(atom.position)

        self._point_mass_dict = point_mass_positions
        return point_mass_positions

    def point_mass_arrays(self) -> dict[int, NDArray[double]]:
        point_mass_arrays: dict[int, NDArray[double]] = {}
        for number, positions in self.point_mass_dict().items():
            point_mass_arrays[number] = array(positions)
        self._point_mass_arrays = point_mass_arrays
        return point_mass_arrays

    def _construct(self):
        if self._constructed is False:
            for unit_cell in self.lattice:
                unit_cell.place_atoms()
                self._constructed = True
                # print("Lattice constructed")
        else:
            print("Lattice already constructed")

    def set_origin(self, vector: NDArray[double]) -> None:
        self.origin = vector

    def rotate_lattice(self, angle: float) -> None:
        self.angle = angle
        for cell in self.lattice:
            cell.angle = angle

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
                collection.append(pos, color=colors[atom.symbol], s=atom.size)
