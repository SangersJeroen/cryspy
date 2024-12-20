from .lattice import CrystalLattice
import string
from .atom import Atom
from numpy._typing import NDArray
from numpy import (
    vstack,
    asarray,
    empty,
    unique,
    identity,
    array,
    double,
    dot,
    min as npmin,
    max as npmax,
)
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
Extent = list[tuple[float, float]]


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

    def point_mass_arrays(self) -> NDArray[double]:
        rows = 0
        for crystal_lattice in self._constituent_structures:
            rows += len(crystal_lattice.unit_cell.atoms) * len(crystal_lattice.lattice)
        datablock = empty((rows, 5))
        s = 0
        for crystal_lattice in self._constituent_structures:
            for unit_cell in crystal_lattice.lattice:
                for atom in unit_cell:
                    datablock[s, :] = array([*atom.position, atom.number, atom.size])
                    s += 1
        self._datablock = datablock
        return datablock

    def repeating_feature(self):
        pass

    def build_model(self):
        for lattice in self._constituent_structures:
            lattice._construct()

    def find_ranges(self) -> list[tuple[float, float]]:
        minpos = npmin(self._datablock, axis=0)
        maxpos = npmax(self._datablock, axis=0)
        self._extremes = [
            (minpos[0], maxpos[0]),
            (minpos[1], maxpos[1]),
            (minpos[2], maxpos[2]),
        ]
        return self._extremes

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
        self.point_mass_arrays()

    def _limit_to_extend(self, extent: Extent) -> NDArray[double]:
        if self._datablock is None:
            data = self.point_mass_arrays()

        data = self._datablock
        xmin, xmax = extent[0]
        ymin, ymax = extent[1]

        xmask: NDArray[bool] = (data[:, 0] >= xmin) & (data[:, 0] <= xmax)
        ymask: NDArray[bool] = (data[:, 1] >= ymin) & (data[:, 1] <= ymax)
        mask: NDArray[bool] = xmask & ymask
        return data[mask]

    def export_to_datfiles(self, filename_stem: str, extent: Extent | None = None):
        """export_point_mass_dict
        Function exports the model layers to a series of plain text .dat
        files that can be used as direct inputs to C. Kittel's* soft-
        ware.

        It does this by slicing the model along the z-axis where all
        atoms in the same layer are written to one layer file.

        for multiple layers the files are named as:
                `filename_stem[a-Z].dat`

        eventually any axis
        """
        if extent is None:
            data: NDArray[double] = self.point_mass_arrays()
            self.find_ranges()
            _, xmax = self._extremes[0]
            _, ymax = self._extremes[1]
        else:
            data = self._limit_to_extend(extent)
            _, xmax = extent[0]
            _, ymax = extent[1]

        unique_z = unique(data[:, 2])
        unique_z.sort()

        unique_w = unique(data[:, -2])
        unique_w.sort()

        slice_thicknesses: list[float] = []
        for i, z in enumerate(unique_z):
            if i + 1 == len(unique_z):
                slice_thicknesses.append(z - unique_z[i - 1])
            else:
                slice_thicknesses.append(unique_z[i + 1] - z)

        for i, (z, t) in enumerate(zip(unique_z, slice_thicknesses)):
            letter = (list(string.ascii_lowercase) + list(string.ascii_uppercase))[i]
            filename = filename_stem + letter + ".dat"

            with open(filename, "w") as file:
                file.write(f"  {xmax:.4f}  {ymax:.4f}  {t:.4f}\n")
                file.write("0\n")

                at_z = data[data[:, 2] == z]
                for w in unique_w:
                    at_w = at_z[at_z[:, -2] == w]
                    if at_w.shape[0] != 0:
                        file.write(
                            f"{int(w)}\n"
                        )  # TODO: Fix this, this should be atomic number...
                        for row_idx in range(at_w.shape[0]):
                            row_dat = at_w[row_idx, :]
                            file.write(
                                f"  {1.0000}  {row_dat[0]/xmax/2:.4f}  {row_dat[1]/ymax/2:.4f}\n"
                            )
                        file.write("\n")
                file.write("\n\n")

    def export_to_xyzfile_temsim(self, filename: str, extent: Extent | None = None):
        """export_point_mass_dict
        Function exports the model layers to a series of plain text .dat
        files that can be used as direct inputs to C. Kittel's* soft-
        ware.

        It does this by slicing the model along the z-axis where all
        atoms in the same layer are written to one layer file.

        for multiple layers the files are named as:
                `filename_stem[a-Z].dat`

        eventually any axis
        """
        data = self.point_mass_arrays()

        if extent is None:
            self.find_ranges()
            xmin, xmax = self._extremes[0]
            ymin, ymax = self._extremes[1]
            zmin, zmax = self._extremes[2]
        else:
            data = self._limit_to_extend(extent)
            xmin, xmax = extent[0]
            ymin, ymax = extent[1]
            zmin, zmax = extent[2]

        print("exporting specimen")
        print(data.shape)
        print(f"xmin, xmax: ({xmin:.3f},{xmax:.3f}")
        print(f"ymin, ymax: ({ymin:.3f},{ymax:.3f}")
        print(f"zmin, zmax: ({zmin:.3f},{zmax:.3f}")

        with open(filename + ".xyz", "w") as file:
            file.write("-- boop boop beep ---\n")
            file.write(f"\t{xmax:.4f}\t{ymax:.4f}\t{zmax:.4f}\n")
            for row_idx in range(data.shape[0]):
                posx, posy, posz, znum, _ = data[row_idx]
                file.write(
                    f"{int(znum)}\t{posx:.4f}\t{posy:.4f}\t{posz:.4f}\t1.000\t0.000\n"
                )
            file.write("-1")
        print(f"{data.shape[0]} atoms exported")

    def export_to_xyzfile_blender(self, filename: str, extent: Extent | None = None):
        """export_point_mass_dict
        Function exports the model layers to a series of plain text .dat
        files that can be used as direct inputs to C. Kittel's* soft-
        ware.

        It does this by slicing the model along the z-axis where all
        atoms in the same layer are written to one layer file.

        for multiple layers the files are named as:
                `filename_stem[a-Z].dat`

        eventually any axis
        """
        data = self.point_mass_arrays()
        num_atoms = data.shape[0]

        with open(filename + ".xyz", "w") as file:
            file.write(f"{num_atoms}\n\n")
            for lattice in self._constituent_structures:
                for ucell in lattice.lattice:
                    for atom in ucell.atoms:
                        file.write(
                            f"{atom.symbol}\t{atom.get_position()[0]:.4f}\t{atom.get_position()[1]:.4f}\t{atom.get_position()[2]:.4f}\n"
                        )
        print(f"{data.shape[0]} atoms exported")

    def plot_3d(self, ranges: Extent | None = None):
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

        if hasattr(self, "_datablock"):
            datablock = self._datablock
        else:
            datablock = self.point_mass_arrays()

        if ranges is not None:
            datablock = self._limit_to_extend(extent=ranges)

        for idx in range(datablock.shape[0]):
            pos = datablock[idx, :3]
            size = datablock[idx, -1]
            color = colors[str(int(datablock[idx, -2]))]
            atoms.append(pos, color=color, s=size)
        figure.show()

    def plot_2d(self, view_axis: Vector = None, ranges: Extent | None = None):
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
            vec1, vec2 = array([17, 2, 5]), array([7, 11, 13])
            vectors = [view_axis, vec1, vec2]
            orth_basis = asarray(gram_schmidt(vectors))
            rotation = linalg.inv(orth_basis)

        if hasattr(self, "_datablock"):
            datablock = self._datablock
        else:
            datablock = self.point_mass_arrays()

        if ranges is not None:
            datablock = self._limit_to_extend(extent=ranges)

        for idx in range(datablock.shape[0]):
            pos = datablock[idx, :3]
            pos = asarray(dot(rotation, pos))
            size = datablock[idx, -1]
            color = colors[str(int(datablock[idx, -2]))]
            atoms.append((pos[0], pos[1], 0), color=color, s=size)
        figure.show()


def gram_schmidt(vectors):
    basis = []
    for v in vectors:
        w = v - sum(dot(v, b) * b for b in basis)
        if (w > 1e-10).any():
            basis.append(w / linalg.norm(w))
    return array(basis)
