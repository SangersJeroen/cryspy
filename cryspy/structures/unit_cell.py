from .atom import Atom
from numpy import array
from typing import Iterator


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
        self.origin = array([x, y, z])

    def place_atoms(self) -> None:
        if self.atoms is None:
            raise ValueError("No atoms in unit cell")
        for atom in self.atoms:
            atom.update_position(atom.position + self.origin)
