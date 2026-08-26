from typing import Any, Literal, assert_type

from minineedle import needle
from minineedle.core import Gap


def assert_lists_equal(list1: list, list2: list) -> None:
    assert len(list1) == len(list2)
    assert all(map(lambda x: x[0] == x[1], zip(list1, list2)))


def test_custom_comparator_is_applied():
    class Item:
        def __init__(self, wild_type: Literal["A", "C", "T", "G"], mutant: Literal["A", "C", "T", "G"]) -> None:
            self.wild_type = wild_type
            self.mutant = mutant

        def __str__(self) -> str:
            if self.wild_type == self.mutant:
                return self.wild_type
            return f"{self.wild_type}->{self.mutant}"

        def __eq__(self, other: Any) -> bool:
            if isinstance(other, Item):
                return self.wild_type == other.wild_type and self.mutant == other.mutant
            else:
                return False

        @staticmethod
        def match_wild_type(left: Item, right: Item) -> bool:
            return left.wild_type == right.wild_type

    seq1 = [
        Item("A", "A"),
        Item("C", "C"),
        Item("G", "T"),
        Item("G", "G"),
        Item("T", "A"),
    ]
    seq2 = [
        Item("A", "A"),
        Item("C", "T"),
        Item("G", "G"),
        Item("T", "T"),
    ]
    expected_nmatrix = [
        [0, -1, -2, -3, -4, -5],
        [-1, 1, 0, -1, -2, -3],
        [-2, 0, 0, -1, -2, -3],
        [-3, -1, -1, -1, 0, -1],
        [-4, -2, -2, -2, -1, -1],
    ]

    alignment = needle.NeedlemanWunsch(seq1, seq2)
    alignment.align()

    assert expected_nmatrix == alignment.get_almatrix()

    custom_alignment = needle.NeedlemanWunsch(seq1, seq2, comparison_function=Item.match_wild_type)
    custom_alignment.align()

    expected_nmatrix = [
        [0, -1, -2, -3, -4, -5],
        [-1, 1, 0, -1, -2, -3],
        [-2, 0, 2, 1, 0, -1],
        [-3, -1, 1, 3, 2, 1],
        [-4, -2, 0, 2, 2, 3],
    ]
    assert expected_nmatrix == custom_alignment.get_almatrix()
