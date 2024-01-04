from .lattice import CrystalLattice

class Specimen:

    def __init__(self, lattices: Iterable[CrystalLattice]) -> None:
        self._constituent_structures = lattices

    def repeating_feature(self):
        pass

    def point_mass_coord_list(self):
        pass

    def project_pmcl(self, plane):
        pass

