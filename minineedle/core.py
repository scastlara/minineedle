from __future__ import annotations

from enum import Enum
from typing import Any, Generic, Literal, Optional, Sequence, overload

from minineedle.typesvars import ItemToAlign


class ScoreMatrix:
    def __init__(self, match: int, miss: int, gap: int) -> None:
        self.match = match
        self.miss = miss
        self.gap = gap

    def __str__(self) -> str:
        return f"Match:{self.match} Missmatch:{self.miss} Gap:{self.gap}"


class AlignmentFormat(str, Enum):
    list = "list"
    str = "str"


class OptimalAlignment(Generic[ItemToAlign]):
    def __init__(self, seq1: Sequence[ItemToAlign], seq2: Sequence[ItemToAlign]) -> None:
        self.seq1 = seq1
        self.seq2 = seq2
        self._alseq1: list[ItemToAlign | Gap] = []
        self._alseq2: list[ItemToAlign | Gap] = []
        self.smatrix = ScoreMatrix(match=1, miss=-1, gap=-1)
        self._score = int()
        self._identity = float()
        self._nmatrix = self._initialize_number_matrix()
        self._pmatrix = self._initialize_pointers_matrix()
        self._gap_character = "-"

        self._is_iterable(self.seq1)
        self._is_iterable(self.seq2)

    def __str__(self) -> str:
        if not self._alseq1:
            self.align()
        return "Alignment of {} and {}:\n\t{}\n\t{}\n".format(
            "SEQUENCE 1",
            "SEQUENCE 2",
            "".join([str(x) for x in self._alseq1]),
            "".join([str(x) for x in self._alseq2]),
        )

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, OptimalAlignment):
            raise ValueError(f"{other} must be instance of {self.__class__.__name__}")

        return self.get_score() < other.get_score()

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, OptimalAlignment):
            raise ValueError(f"{other} must be instance of {self.__class__.__name__}")

        return self.get_score() == other.get_score()

    @property
    def gap_character(self) -> str:
        return self._gap_character

    @gap_character.setter
    def gap_character(self, var: str) -> None:
        self._gap_character = var
        if self._alseq1:
            self._alseq1 = self._change_gap_char(self._alseq1)
            self._alseq2 = self._change_gap_char(self._alseq2)

    def _is_iterable(self, iterable: Sequence[ItemToAlign]) -> None:
        iter(iterable)

    def get_score(self) -> int | float:
        if not self._alseq1:
            self.align()
        return self._score

    def change_matrix(self, newmatrix: ScoreMatrix) -> None:
        """
        Changes ScoreMatrix

        Args:
            newmatrix (ScoreMatrix): Matrix containing match, miss, and gap penalties.
        """
        if isinstance(newmatrix, ScoreMatrix):
            self.smatrix = newmatrix
        else:
            raise ValueError("New matrix should be a ScoreMatrix object.")

    def align(self) -> None:
        """
        Performs a Needleman-Wunsch or Smith-Waterman alignment with the given sequences and the
        corresponding ScoreMatrix.
        """
        self._add_initial_pointers()
        self._add_gap_penalties()
        self._fill_matrices()

        imax, jmax = self._get_last_cell_position()
        self._get_alignment_score(imax, jmax)
        self._trace_back_alignment(imax, jmax)

    def get_almatrix(self) -> list[list[int]]:
        """
        Returns the alignment matrix (list of lists)
        """
        if not self._alseq1:
            self.align()
        return self._nmatrix

    def get_identity(self) -> float:
        """
        Returns the % of identity of the alignment
        """
        if not self._alseq1:
            self.align()
        return round(self._identity, 2)  # Two decimal points

    @overload
    def get_aligned_sequences(self, sequence_format: Literal[AlignmentFormat.str] | Literal["str"]) -> tuple[str, str]:
        ...

    @overload
    def get_aligned_sequences(
        self, sequence_format: Literal[AlignmentFormat.list] | Literal["list"] = "list"
    ) -> tuple[list[ItemToAlign | Gap], list[ItemToAlign | Gap]]:
        ...

    def get_aligned_sequences(
        self, sequence_format: Literal["str"] | AlignmentFormat | Literal["list"] = "list"
    ) -> tuple[str, str] | tuple[list[ItemToAlign | Gap], list[ItemToAlign | Gap]]:
        """
        Returns tuple with both aligned sequences as lists or as strings.
        """
        if sequence_format == AlignmentFormat.list:
            return self._change_gap_char(self._alseq1), self._change_gap_char(self._alseq2)
        elif sequence_format == AlignmentFormat.str:
            return "".join([str(x) for x in self._change_gap_char(self._alseq1)]), "".join(
                [str(x) for x in self._change_gap_char(self._alseq2)]
            )
        else:
            raise ValueError("Sequence_format has to be either 'list' or 'str'!")

    def _change_gap_char(self, iterable: list[ItemToAlign | Gap]) -> list[ItemToAlign | Gap]:
        new_sequence: list[ItemToAlign | Gap] = []
        for it in iterable:
            if isinstance(it, Gap):
                it.character = self._gap_character
            new_sequence.append(it)
        return new_sequence

    def _add_initial_pointers(self) -> None:
        """
        Fills the pointers matrix first row with "left" pointer and
        the initial column with "up" pointers.
        """
        for row in self._pmatrix:
            row[0] = "up"

        for jcol in range(0, len(self._pmatrix[0])):
            self._pmatrix[0][jcol] = "left"

        self._pmatrix[0][0] = None

    def _add_gap_penalties(self) -> None:
        """
        Fills number matrix first row and first column with the gap penalties.
        Will change depending if using NeedlemanWunsch or SmithWaterman.
        """
        raise NotImplementedError("NeedlemanWunsch or SmithWaterman should be used instead!")

    def _get_alignment_score(self, imax: int, jmax: int) -> None:
        """
        Stores the alignment score value in the _score attribute.
        """
        self._score = self._nmatrix[imax][jmax]

    def _get_last_cell_position(self) -> tuple[int, int]:
        """
        Returns the cell row and column of the last cell in the matrix in which
        the alignment ends. For Needleman-Wunsch this will be the last cell of the matrix,
        for Smith-Waterman will be the cell with the highest score.
        """
        raise NotImplementedError("NeedlemanWunsch or SmithWaterman should be used instead!")

    def _initialize_number_matrix(self) -> list[list[int]]:
        """
        Initializes the matrix where the computed scores are stored.
        """
        return [[0 for x in range(len(self.seq1) + 1)] for x in range(len(self.seq2) + 1)]

    def _initialize_pointers_matrix(self) -> list[list[Optional[str]]]:
        """
        Initializes the matrix where the "up", "down", "diag" pointers are stored.
        """
        return [[None for x in range(len(self.seq1) + 1)] for x in range(len(self.seq2) + 1)]

    def _fill_matrices(self) -> None:
        for irow in range(0, len(self.seq2)):
            for jcol in range(0, len(self.seq1)):
                # Scores
                topscore = self._nmatrix[irow][jcol + 1] + self.smatrix.gap
                leftscore = self._nmatrix[irow + 1][jcol] + self.smatrix.gap
                diagscore = self._nmatrix[irow][jcol]
                if self.seq1[jcol] == self.seq2[irow]:
                    diagscore += self.smatrix.match
                else:
                    diagscore += self.smatrix.miss

                self._check_best_score(diagscore, topscore, leftscore, irow, jcol)

    def _trace_back_alignment(self, irow: int, jcol: int) -> None:
        while True:
            if self._pmatrix[irow][jcol] == "diag":
                self._alseq1.append(self.seq1[jcol - 1])
                self._alseq2.append(self.seq2[irow - 1])
                if self.seq1[jcol - 1] == self.seq2[irow - 1]:
                    self._identity += 1
                irow -= 1
                jcol -= 1
            elif self._pmatrix[irow][jcol] == "up":
                self._alseq1.append(Gap(self._gap_character))
                self._alseq2.append(self.seq2[irow - 1])
                irow -= 1
            elif self._pmatrix[irow][jcol] == "left":
                self._alseq1.append(self.seq1[jcol - 1])
                self._alseq2.append(Gap(self._gap_character))
                jcol -= 1
            else:
                break
        self._alseq1 = list(reversed(self._alseq1))
        self._alseq2 = list(reversed(self._alseq2))
        self._identity = (self._identity / len(self._alseq1)) * 100

    def _check_best_score(self, diagscore: int, topscore: int, leftscore: int, irow: int, jcol: int) -> None:
        """
        Decides best score for a given cell. Will change depending if using
        Needleman-Wunsch or Smith-Waterman
        """
        raise NotImplementedError


class Gap(object):
    def __init__(self, character: str = "-") -> None:
        self.character = character

    def __str__(self) -> str:
        return str(self.character)

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, Gap)
