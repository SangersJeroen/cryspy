from .lattice import CrystalLattice
from .atom import Atom
from numpy._typing import NDArray
from numpy import array, double, dot, matrix, identity
import numpy.linalg as linalg
from ..crystals.colors import colors


from glumpy.graphics.text import FontManager
from glumpy.api.matplotlib import (
    Figure,
    PanZoom,
    PointCollection,
    LinearScale,
    PointCollection,
    Trackball,
)


Vector = tuple[float, float, float] | list[float] | NDArray[double]
Position = tuple[float, float, float] | list[float] | NDArray[double]
PointMassDict = dict[int, list[Position]]

regular = FontManager.get("OpenSans-Regular.ttf")
bold = FontManager.get("OpenSans-Bold.ttf")


class Specimen:
    def __init__(self, lattices: list[CrystalLattice]) -> None:
        self._constituent_structures = lattices

    def point_mass_dict(self) -> PointMassDict:
        self._atoms: list[Atom] = []
        for lattice in self._constituent_structures:
            for unit_cell in lattice.lattice:
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

    def repeating_feature(self):
        pass

    def build_model(self):
        for lattice in self._constituent_structures:
            lattice._construct()

    def find_ranges(self) -> tuple[float, float, float]:
        xmax: float = 0.0
        ymax: float = 0.0
        zmax: float = 0.0

        for atom in self.point_mass_arrays().values():
            xmax = max(xmax, max(atom[:, 0]))
            ymax = max(ymax, max(atom[:, 1]))
            zmax = max(zmax, max(atom[:, 2]))
        self._extremes = (xmax, ymax, zmax)
        return (xmax, ymax, zmax)

    def project_onto_plane(self, plane_vector_normal: Vector):
        unit_normal: Vector = array(plane_vector_normal)

        for lattice in self._constituent_structures:
            for unit_cell in lattice.lattice:
                unit_cell.place_atoms()
                for atom in unit_cell.atoms:
                    old_pos = atom.get_position()
                    new_pos = (
                        old_pos
                        - dot(old_pos, unit_normal)
                        / dot(unit_normal, unit_normal)
                        * unit_normal
                    )
                    atom.update_position(new_pos)

    def export_point_mass_dict(self, filename_stem: str):
        """export_point_mass_dict
        Function exports the model layers to a series of plain text .dat
        files that can be used as direct inputs to C. Kittel's* soft-
        ware.

        It does this by slicing the model along the specified axis where
        all atoms in the same layer are written to one layer file.

        for multiple layers the files are named as:
                `filename_stem[a-Z].dat`
        """
        extremes = self.find_ranges()

        with open(filename, "w") as file:
            file.write(
                f"  {extremes[0]:.4f}  {extremes[1]:.4f}  {layer_thickness:.4f}\n"
            )
            file.write(" 0\n")
            file.write("\n")
            for atom_number in self._point_mass_arrays:
                positions = self._point_mass_arrays[atom_number]
                file.write(f"{atom_number}\n")
                for position in positions:
                    file.write(
                        f"  {position[0]/extremes[0]:.4f}  {position[1]/extremes[1]:.4f}  {position[2]/extremes[2]:.4f}\n"
                    )
                file.write("\n")
            file.write("\n")

    def plot_3d(self):
        figure = Figure()
        axis = figure.add_axes(
            aspect=1,
            rect=(0, 0, 1, 1),
            interface=Trackball(name="trackball"),
            zscale=LinearScale(domain=[0, 1], range=[-1 / 50, 1 / 50]),
            xscale=LinearScale(domain=[0, 1], range=[-1 / 50, 1 / 50]),
            yscale=LinearScale(domain=[0, 1], range=[-1 / 50, 1 / 50]),
            facecolor=(0, 0, 0, 1),
        )
        atoms = PointCollection(
            mode="agg",
            color="local",
            size="local",
        )
        axis.add_drawable(atoms)

        self.build_model()
        for lattice in self._constituent_structures:
            for unit_cell in lattice.lattice:
                for atom in unit_cell.atoms:
                    pos = atom.get_position()
                    size = atom.size
                    color = colors[atom.symbol]
                    atoms.append(pos, color=color, s=size)
        figure.show()

    def plot_2d(self, view_axis: Vector = None):
        figure = Figure()
        axis = figure.add_axes(
            aspect=1,
            rect=(0, 0, 1, 1),
            interface=PanZoom(name="panzoom", aspect=1),
            xscale=LinearScale(domain=[0, 1], range=[-1 / 50, 1 / 50]),
            yscale=LinearScale(domain=[0, 1], range=[-1 / 50, 1 / 50]),
            facecolor=(0, 0, 0, 1),
        )
        atoms = PointCollection(
            mode="agg",
            color="local",
            size="local",
        )
        axis.add_drawable(atoms)

        rotation = identity(3)
        if view_axis is not None:
            vec1, vec2 = array([-1, 2, 1]), array([2, -1, -1])
            vectors = [view_axis, vec1, vec2]
            orth_basis = matrix(gram_schmidt(vectors))
            print(orth_basis)
            rotation = linalg.inv(orth_basis)

        self.build_model()
        for lattice in self._constituent_structures:
            for unit_cell in lattice.lattice:
                for atom in unit_cell.atoms:
                    pos = atom.get_position()
                    pos = dot(rotation, pos)
                    size = atom.size
                    color = colors[atom.symbol]
                    atoms.append(pos, color=color, s=size)
        figure.show()


def gram_schmidt(vectors):
    basis = []
    for v in vectors:
        w = v - sum(dot(v, b) * b for b in basis)
        if (w > 1e-10).any():
            basis.append(w / linalg.norm(w))
    return array(basis)
