import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from typing import Iterator
from copy import deepcopy


HEX_VEC = np.array(
    [
        [0, 1, 0],
        [np.cos(30 / 180 * np.pi), np.sin(30 / 180 * np.pi), 0],
        [np.cos(-60 / 180 * np.pi), np.sin(-60 / 180 * np.pi), 0],
    ]
)

CUBIC_VEC = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


class Atom:
    def __init__(
        self, name: str, symbol: str, number: int, x: float, y: float, z: float
    ) -> None:
        self.name = name
        self.symbol = symbol
        self.number = number
        self.position = np.array([x, y, z])

    def __repr__(self) -> str:
        return f"Atom({self.name}, {self.symbol}, {self.number}, {self.position})"

    def update_position(self, new_position: np.ndarray[float]) -> None:
        self.position = new_position

    def get_position(self) -> np.ndarray[float]:
        return self.position


class UnitCell:
    def __init__(self, atoms: list[Atom] | None) -> None:
        self.atoms = atoms

    def __repr__(self) -> str:
        return f"Origin: {self.origin}, Atoms: {self.atoms}"

    def __iter__(self) -> Iterator[Atom]:
        if self.atoms is not None:
            return iter(self.atoms)
        else:
            raise ValueError("No atoms in unit cell")

    def add_atoms(self, atoms: list[Atom]) -> None:
        if self.atoms is not None:
            self.atoms.extend(atoms)
        else:
            self.atoms = atoms

    def set_origin(self, x: float, y: float, z: float) -> None:
        self.origin = np.array([x, y, z])

    def place_atoms(self) -> None:
        if self.atoms is None:
            raise ValueError("No atoms in unit cell")
        for atom in self.atoms:
            atom.update_position(atom.position + self.origin)


class CrystalLattice:
    def __init__(self, unit_cell: UnitCell) -> None:
        self.unit_cell = unit_cell

    def construct_lattice(
        self, vectors: np.ndarray[float], repetitions: tuple[int, int, int]
    ) -> None:
        self.lattice = [deepcopy(self.unit_cell) for _ in range(np.prod(repetitions))]

        if vectors.shape[0] != len(repetitions):
            raise ValueError("Number of Vectors and repetitions must have same length")

        xreps, yreps, zreps = repetitions
        xrep_range = np.arange(xreps) - (xreps - 1) / 2
        yrep_range = np.arange(yreps) - (yreps - 1) / 2
        zrep_range = np.arange(zreps) - (zreps - 1) / 2
        xidx, yidx, zidx = np.meshgrid(
            xrep_range, yrep_range, zrep_range, indexing="ij"
        )

        xidx = xidx.flatten()
        yidx = yidx.flatten()
        zidx = zidx.flatten()

        for index, cell in enumerate(self.lattice):
            position = np.dot(
                vectors.T, np.array([xidx[index], yidx[index], zidx[index]])
            )
            cell.set_origin(*position)
            cell.place_atoms()

    def rotate_lattice(self, angle: float) -> None:
        rotation_matrix = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1],
            ]
        )
        for cell in self.lattice:
            cell.origin = rotation_matrix.dot(cell.origin)
            cell.place_atoms()

    def plot_lattice_3d(self) -> None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        for cell in self.lattice:
            for atom in cell:
                ax.scatter(*atom.get_position())

        # equal axis
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        plt.show()

    def plot_lattice(self, ax: plt.Axes, colors: dict[str, str]) -> None:
        for cell in self.lattice:
            for atom in cell:
                pos = atom.get_position()
                ax.scatter(pos[0], pos[1], color=colors[atom.symbol], s=5)

        # equal axis
        ax.set_aspect("equal")


if __name__ == "__main__":
    tungsten = Atom("W", "W", 74, 2 / 3, 1 / 3, 0)
    selenium0 = Atom("Se", "Se", 34, 1 / 3, 2 / 3, 0.87)
    selenium1 = Atom("Se", "Se", 34, 1 / 3, 2 / 3, -0.87)
    unit_cell = UnitCell([tungsten, selenium0, selenium1])

    lattice0 = CrystalLattice(unit_cell)
    lattice0.construct_lattice(HEX_VEC, (15, 15, 1))

    tungsten_ = Atom("W", "W", 74, 2 / 3, 1 / 3, 0)
    selenium0_ = Atom("Se", "Se", 34, 1 / 3, 2 / 3, 0.87)
    selenium1_ = Atom("Se", "Se", 34, 1 / 3, 2 / 3, -0.87)
    unit_cell_ = UnitCell([tungsten_, selenium0_, selenium1_])

    lattice1 = CrystalLattice(unit_cell_)
    lattice1.construct_lattice(HEX_VEC, (15, 15, 1))
    lattice1.rotate_lattice(30 / 180 * np.pi)

    colors = {
        "W": "black",
        "Se": "Yellow",
    }
    fig, ax = plt.subplots()
    lattice0.plot_lattice(ax, colors)
    lattice1.plot_lattice(ax, colors)
    plt.show()
