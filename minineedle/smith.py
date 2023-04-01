from typing import Optional, Sequence

from minineedle.core import OptimalAlignment
from minineedle.typesvars import ItemToAlign


class SmithWaterman(OptimalAlignment[ItemToAlign]):
    """
    Smith-Waterman algorithm
    """

    def __init__(self, seq1: Sequence[ItemToAlign], seq2: Sequence[ItemToAlign]) -> None:
        super().__init__(seq1, seq2)

    def _add_gap_penalties(self) -> None:
        """
        Fills number matrix first row and first column with the gap penalties.
        """
        for i in range(1, len(self.seq1) + 1):
            self._nmatrix[0][i] = 0

        for j in range(1, len(self.seq2) + 1):
            self._nmatrix[j][0] = 0

    def _get_last_cell_position(self) -> tuple[int, int]:
        """
        Returns the cell row and column of the last cell in the matrix in which
        the alignment ends. For Needleman-Wunsch this will be the last cell of the matrix,
        for Smith-Waterman will be the cell with the highest score.
        """
        imax, jmax = 0, 0
        max_score = 0
        for irow in range(0, len(self._nmatrix)):
            for jcol in range(0, len(self._nmatrix[0])):
                score = self._nmatrix[irow][jcol]
                if score > max_score:
                    imax = irow
                    jmax = jcol
                    max_score = score
        return imax, jmax

    def _check_best_score(self, diagscore: int, topscore: int, leftscore: int, irow: int, jcol: int) -> None:
        best_pointer: Optional[str] = ""
        best_score = 0

        if diagscore >= topscore:
            if diagscore >= leftscore:
                best_pointer, best_score = ("diag", diagscore)
            else:
                best_pointer, best_score = ("left", leftscore)
        else:
            if topscore > leftscore:
                best_pointer, best_score = ("up", topscore)
            else:
                best_pointer, best_score = ("left", leftscore)

        if best_score < 0:
            best_pointer = None
            best_score = 0

        self._pmatrix[irow + 1][jcol + 1] = best_pointer
        self._nmatrix[irow + 1][jcol + 1] = best_score
