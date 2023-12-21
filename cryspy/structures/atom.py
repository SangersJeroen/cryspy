from numpy import array, double
from numpy._typing import NDArray


class Atom:
    def __init__(
        self,
        name: str,
        symbol: str,
        number: int,
        x: float,
        y: float,
        z: float,
        size: float | None = None,
    ) -> None:
        self.name = name
        self.symbol = symbol
        self.number = number
        self.position = array([x, y, z])
        self.size = size if not None else 5

    def __repr__(self) -> str:
        return f"Atom({self.name}, {self.symbol}, {self.number}, {self.position})"

    def update_position(self, new_position: NDArray[double]) -> None:
        self.position = new_position

    def get_position(self) -> NDArray[double]:
        return self.position
