from .lattice import CrystalLattice
from .atom import Atom
from numpy._typing import NDArray
from numpy import array, double, dot
from numpy.linalg import norm

Vector = tuple[float, float, float] | list[float] | NDArray[double]
Position = tuple[float, float, float] | list[float] | NDArray[double]
PointMassDict = dict[int, list[Position]]


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

        return point_mass_positions

    def point_mass_arrays(self) -> dict[int, NDArray[double]]:
        point_mass_arrays: dict[int, NDArray[double]] = {}
        for number, positions in self.point_mass_dict().items():
            point_mass_arrays[number] = array(positions)
        self._point_mass_arrays = point_mass_arrays
        return point_mass_arrays

    def repeating_feature(self):
        pass

    def project_onto_plane(self, plane_vector_normal: Vector):
        unit_normal: Vector = array(plane_vector_normal)

        for lattice in self._constituent_structures:
            for unit_cell in lattice.lattice:
                for atom in unit_cell.atoms:
                    old_pos = atom.get_position()
                    new_pos = (
                        old_pos
                        - dot(old_pos, unit_normal)
                        / dot(unit_normal, unit_normal)
                        * unit_normal
                    )
                    atom.update_position(new_pos)
                    print(old_pos, new_pos)

    def export_point_mass_dict(self, filename: str):
        with open(filename, "w") as file:
            for number, positions in self.point_mass_dict().items():
                for position in positions:
                    file.write(f"{number} {position[0]} {position[1]} {position[2]}\n")
