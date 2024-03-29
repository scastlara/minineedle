import pytest
from minineedle import core, needle, smith
from typing_extensions import assert_type


def test_iterable_true() -> None:
    """
    Tests if iterable objects are accepted by NeedlemanWunsch.
    """
    seq1 = "ACTG"
    seq2 = "ACTG"
    alignment = needle.NeedlemanWunsch(seq1, seq2)
    assert_type(alignment, needle.NeedlemanWunsch[str])


def test_iterable_true_list() -> None:
    """
    Tests if iterable (lists) objects are accepted by NeedlemanWunsch.
    """
    seq1 = ["A", "C", "T", "G"]
    seq2 = ["A", "C", "T", "G"]
    alignment = needle.NeedlemanWunsch(seq1, seq2)
    assert_type(alignment, needle.NeedlemanWunsch[str])


def test_iterable_different() -> None:
    """
    Tests if iterable (lists) objects are accepted by NeedlemanWunsch.
    """
    seq1 = ["A", "C", "T", "G"]
    seq2 = "ACTG"
    alignment = needle.NeedlemanWunsch(seq1, seq2)
    assert_type(alignment, needle.NeedlemanWunsch[str])


def test_change_matrix_true() -> None:
    """
    Tests if ScoreMatrix object is accepted by NeedlemanWunsch's change_matrix
    method.
    """
    seq1 = "ACTG"
    seq2 = "ACTG"
    matrix = core.ScoreMatrix(1, -1, -1)
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)

    needle_alignment.change_matrix(matrix)


def test_needleman_dynamic_matrix() -> None:
    """
    Checks if alignment score matrix is correct.
    """
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()
    expected_matrix = [
        [0, -1, -2, -3, -4, -5, -6, -7],
        [-1, 1, 0, -1, -2, -3, -4, -5],
        [-2, 0, 0, 1, 0, -1, -2, -3],
        [-3, -1, -1, 0, 2, 1, 0, -1],
        [-4, -2, -2, -1, 1, 1, 0, -1],
        [-5, -3, -3, -1, 0, 0, 0, -1],
        [-6, -4, -2, -2, -1, -1, 1, 0],
        [-7, -5, -3, -1, -2, -2, 0, 0],
    ]
    assert needle_alignment._nmatrix == expected_matrix


def test_needleman_pointer_matrix() -> None:
    """
    Tests if pointer matrix is correct.
    """
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()
    expected_matrix = [
        [None, "left", "left", "left", "left", "left", "left", "left"],
        ["up", "diag", "left", "left", "left", "diag", "left", "left"],
        ["up", "up", "diag", "diag", "left", "left", "left", "left"],
        ["up", "up", "diag", "up", "diag", "left", "left", "left"],
        ["up", "up", "diag", "up", "diag", "diag", "diag", "diag"],
        ["up", "up", "diag", "diag", "up", "diag", "diag", "diag"],
        ["up", "up", "diag", "up", "up", "diag", "diag", "left"],
        ["up", "up", "up", "diag", "left", "diag", "up", "diag"],
    ]
    assert needle_alignment._pmatrix == expected_matrix


def test_needleman_alignment() -> None:
    """
    Checks if align method works
    """
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()

    assert needle_alignment._alseq1 == ["G", "C", "A", core.Gap(), "T", "G", "C", "U"]
    assert needle_alignment._alseq2 == ["G", core.Gap(), "A", "T", "T", "A", "C", "A"]


def test_needleman_alignment_score() -> None:
    """
    Tests if score is correctly computed.
    """
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()
    assert needle_alignment.get_score() == 0


def test_needleman_alignment_score_list() -> None:
    """
    Tests if score is correctly computed.
    """
    seq1 = ["G", "C", "A", "T", "G", "C", "U"]
    seq2 = ["G", "A", "T", "T", "A", "C", "A"]
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()

    assert needle_alignment.get_score() == 0


def test_needleman_alignment_score_different() -> None:
    """
    Tests if score is correctly computed.
    """
    seq1 = ["G", "C", "A", "T", "G", "C", "U"]
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()

    assert needle_alignment.get_score() == 0


def test_needleman_alignment_integers() -> None:
    """
    Tests if algorithm can handle integers
    """
    seq1 = [1, 2, 3, 5, 1]
    seq2 = [1, 2, 9, 9, 9, 3, 5, 1]
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()

    assert needle_alignment.get_score() == 2


def test_needleman_alignment_residue_mix() -> None:
    """
    Tests if algorithm can handle different types in list
    """
    seq1 = [1, "C", "A", 2, "G", "C", "U"]
    seq2 = ["G", "A", "T", "T", "A", "C", "A"]
    # mypy will correctly complain about this
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))

    needle_alignment.align()


def test_needleman_alignment_identity() -> None:
    """
    Tests if score is correctly computed.
    """
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()

    assert needle_alignment.get_identity() == 50.0


def test_needleman_to_string() -> None:
    """
    Tests string formatting of alignment object.
    """
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle_alignment.align()

    assert str(needle_alignment) == """Alignment of SEQUENCE 1 and SEQUENCE 2:\n\tGCA-TGCU\n\tG-ATTACA\n"""


def test_needleman_comparison_equal() -> None:
    """
    Tests score comparison of equal alignments
    """
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle1 = needle.NeedlemanWunsch(seq1, seq2)
    needle1.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle1.align()
    needle2 = needle.NeedlemanWunsch(seq1, seq2)
    needle2.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle2.align()

    assert needle1 == needle2


def test_needleman_comparison_lower() -> None:
    """
    Tests score comparison of different alignments
    """
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    seq3 = "LLLLLLL"
    needle1 = needle.NeedlemanWunsch(seq1, seq2)
    needle1.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle1.align()
    needle2 = needle.NeedlemanWunsch(seq1, seq3)
    needle2.change_matrix(core.ScoreMatrix(1, -1, -1))
    needle2.align()

    assert needle1 > needle2


def test_needleman_get_aligned_sequences_list() -> None:
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.align()
    seq_a, seq_b = needle_alignment.get_aligned_sequences()

    assert isinstance(seq_a, list)
    assert isinstance(seq_b, list)


def test_needleman_get_aligned_sequences_string() -> None:
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.align()
    seq_a, seq_b = needle_alignment.get_aligned_sequences("str")

    assert isinstance(seq_a, str)
    assert isinstance(seq_b, str)


def test_needleman_get_aligned_sequences_wrong() -> None:
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.align()

    with pytest.raises(ValueError):
        needle_alignment.get_aligned_sequences("wrong")  # type: ignore[call-overload]


def test_smith_dynamic_matrix() -> None:
    """
    Checks if alignment score matrix is correct.
    """
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    smith_alignment = smith.SmithWaterman(seq1, seq2)
    smith_alignment.change_matrix(core.ScoreMatrix(3, -3, -2))
    smith_alignment.align()
    expected_matrix = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 3, 1, 0, 0, 0, 3, 3],
        [0, 0, 3, 1, 0, 0, 0, 3, 6],
        [0, 3, 1, 6, 4, 2, 0, 1, 4],
        [0, 3, 1, 4, 9, 7, 5, 3, 2],
        [0, 1, 6, 4, 7, 6, 4, 8, 6],
        [0, 0, 4, 3, 5, 10, 8, 6, 5],
        [0, 0, 2, 1, 3, 8, 13, 11, 9],
        [0, 3, 1, 5, 4, 6, 11, 10, 8],
        [0, 1, 0, 3, 2, 7, 9, 8, 7],
    ]
    assert smith_alignment._nmatrix == expected_matrix


def test_smith_alignment() -> None:
    """
    Checks if align method works
    """
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    smith_alignment = smith.SmithWaterman(seq1, seq2)
    smith_alignment.change_matrix(core.ScoreMatrix(3, -3, -2))
    smith_alignment.align()
    assert smith_alignment.get_aligned_sequences()[0] == [
        "G",
        "T",
        "T",
        core.Gap(),
        "A",
        "C",
    ]
    assert smith_alignment.get_aligned_sequences()[1] == ["G", "T", "T", "G", "A", "C"]


def test_smith_alignment_score() -> None:
    """
    Checks if align method works
    """
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    smith_alignment = smith.SmithWaterman(seq1, seq2)
    smith_alignment.change_matrix(core.ScoreMatrix(3, -3, -2))
    smith_alignment.align()

    assert smith_alignment.get_score() == 13


def test_gap_char_change() -> None:
    """
    Checks change of Gap character
    """
    gap = core.Gap()
    gap.character = "imagap"

    assert str(gap) == "imagap"


def test_gap_change_alignment() -> None:
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(3, -3, -2))
    needle_alignment.align()
    needle_alignment.gap_character = "-gap-"
    aligned_seq1, _ = needle_alignment.get_aligned_sequences(sequence_format="str")

    assert aligned_seq1 == "TGTT-gap-ACGG"


def test_gap_change_alignment_before() -> None:
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(3, -3, -2))
    needle_alignment.gap_character = "-gap-"
    needle_alignment.align()
    aligned_seq1, _ = needle_alignment.get_aligned_sequences(sequence_format="str")

    assert aligned_seq1 == "TGTT-gap-ACGG"


def test_strange_alignment() -> None:
    seq1 = "TG--TA--CTA"
    seq2 = "GG--TGA--CTA"
    needle_alignment = needle.NeedlemanWunsch(seq1, seq2)
    needle_alignment.change_matrix(core.ScoreMatrix(2, -2, -3))
    needle_alignment.align()
    needle_alignment.gap_character = "+"

    assert needle_alignment.get_aligned_sequences("str")[0] == "TG--T+A--CTA"


def test_get_gap_character() -> None:
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    alignment = needle.NeedlemanWunsch(seq1, seq2)

    assert alignment.gap_character == "-"


def test_gap_equal() -> None:
    g1 = core.Gap()
    g2 = core.Gap()

    assert g1 == g2


def test_gap_char() -> None:
    gap = core.Gap("a")

    assert str(gap) == "a"
