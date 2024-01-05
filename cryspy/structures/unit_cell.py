from .atom import Atom
from numpy import array, sin, cos
from typing import Iterator


class UnitCell:
    def __init__(self, atoms: list[Atom] | None) -> None:
        self.atoms = atoms
        self.angle = 0
        self.origin = array([0, 0, 0])
        self.atoms_placed = False

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
        self.origin = array([x, y, z])

    def place_atoms(self) -> None:
        if self.atoms is None:
            raise ValueError("No atoms in unit cell")
        for atom in self.atoms:
            rotation_matrix = array(
                [
                    [cos(self.angle), -sin(self.angle), 0],
                    [sin(self.angle), cos(self.angle), 0],
                    [0, 0, 1],
                ]
            )
            atom.update_position(rotation_matrix.dot(atom.position) + self.origin)
