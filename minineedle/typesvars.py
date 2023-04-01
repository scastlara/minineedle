from __future__ import annotations

from typing import Any, Protocol, TypeVar


class Aligneable(Protocol):
    def __str__(self) -> str:
        ...

    def __eq__(self, obj: Any) -> bool:
        ...


ItemToAlign = TypeVar("ItemToAlign", bound=Aligneable)
