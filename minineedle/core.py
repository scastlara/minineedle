import sys

class OptimalAlignment(object):
    def __init__(self, seq1, seq2):
        self.seq1     = seq1
        self.seq2     = seq2
        self._alseq1   = []
        self._alseq2   = []
        self.smatrix  = ScoreMatrix(match=1, miss=-1, gap=-1)
        self._score    = int()
        self._identity = float()
        self._nmatrix  = self._initialize_number_matrix()
        self._pmatrix  = self._initialize_pointers_matrix()
        self._gap_character = "-"

        self._is_iterable(self.seq1)
        self._is_iterable(self.seq2)

    
    @property
    def gap_character(self):
        return self._gap_character
    
    @gap_character.setter
    def gap_character(self, var):
        self._gap_character = var
        if self._alseq1:
            self._alseq1 = self._change_gap_char(self._alseq1)
            self._alseq2 = self._change_gap_char(self._alseq2)

    def _is_iterable(self, iterable):
        iter(iterable)

    def __str__(self):
        if not self._alseq1:
            self.align()
        return "Alignment of {} and {}:\n\t{}\n\t{}\n".format("SEQUENCE 1",
                                                              "SEQUENCE 2",
                                                              "".join([ str(x) for x in self._alseq1 ]),
                                                              "".join([ str(x) for x in self._alseq2 ])
                                                              )

    def __lt__(self, other):
        return self.get_score() < other.get_score()

    def __eq__(self, other):
        return self.get_score() == other.get_score()

    def get_score(self):
        if not self._alseq1:
            self.align()
        return self._score

    def change_matrix(self, newmatrix):
        """
        Changes ScoreMatrix

        Args:
            newmatrix (ScoreMatrix): Matrix containing match, miss, and gap penalties.
        """
        if isinstance(newmatrix, ScoreMatrix):
            self.smatrix = newmatrix
        else:
            raise ValueError("New matrix should be a ScoreMatrix object.")

    def align(self):
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

    def get_almatrix(self):
        """
        Returns the alignment matrix (list of lists)
        """
        if not self._alseq1:
            self.align()
        return self._nmatrix
    
    def get_identity(self):
        """
        Returns the % of identity of the alignment
        """
        if not self._alseq1:
            self.align()
        return round(self._identity, 2) # Two decimal points

    def get_aligned_sequences(self, sequence_format="list"):
        """
        Returns tuple with both aligned sequences as lists or as strings.
        """
        if sequence_format == "list":
            seq_a = self._change_gap_char(self._alseq1)
            seq_b = self._change_gap_char(self._alseq2)
        elif sequence_format == "string" or sequence_format == "str":
            seq_a = "".join([ str(x) for x in self._change_gap_char(self._alseq1) ])
            seq_b = "".join([ str(x) for x in self._change_gap_char(self._alseq2) ])
        else:
            raise ValueError("Sequence_format has to be either 'list' or 'str'!")
        return seq_a, seq_b

    def _change_gap_char(self, iterable):
        new_sequence = []
        for it in iterable:
            if isinstance(it, Gap):
                it.character = self._gap_character
            new_sequence.append(it)
        return new_sequence

    def _add_initial_pointers(self):
        """
        Fills the pointers matrix first row with "left" pointer and 
        the initial column with "up" pointers.
        """
        for row in self._pmatrix:
            row[0] = "up"

        for jcol in range(0, len(self._pmatrix[0])):
            self._pmatrix[0][jcol] = "left"

        self._pmatrix[0][0] = None

    def _add_gap_penalties(self):
        """
        Fills number matrix first row and first column with the gap penalties.
        Will change depending if using NeedlemanWunsch or SmithWaterman.
        """
        raise NotImplementedError("NeedlemanWunsch or SmithWaterman should be used instead!")
        

    def _get_alignment_score(self, imax, jmax):
        """
        Stores the alignment score value in the _score attribute.
        """
        self._score = self._nmatrix[imax][jmax]
    
    def _get_last_cell_position(self):
        """
        Returns the cell row and column of the last cell in the matrix in which 
        the alignment ends. For Needleman-Wunsch this will be the last cell of the matrix,
        for Smith-Waterman will be the cell with the highest score.
        """
        raise NotImplementedError("NeedlemanWunsch or SmithWaterman should be used instead!")

    def _initialize_number_matrix(self):
        """
        Initializes the matrix where the computed scores are stored.
        """
        return [ [ 0 for x in range(len(self.seq1) + 1)] for x in range(len(self.seq2) +1)]

    def _initialize_pointers_matrix(self):
        """
        Initializes the matrix where the "up", "down", "diag" pointers are stored.
        """
        return [ [ None for x in range(len(self.seq1) + 1)] for x in range(len(self.seq2) +1)]

    def _fill_matrices(self):
        for irow in range(0, len(self.seq2)):
            for jcol in range(0, len(self.seq1)):
                # Scores
                topscore  = self._nmatrix[irow][jcol + 1] + self.smatrix.gap
                leftscore = self._nmatrix[irow + 1][jcol] + self.smatrix.gap
                diagscore = self._nmatrix[irow][jcol]
                if self.seq1[jcol] == self.seq2[irow]:
                    diagscore += self.smatrix.match
                else:
                    diagscore += self.smatrix.miss

                self._check_best_score(diagscore, topscore, leftscore, irow, jcol)
        return

    def _trace_back_alignment(self, irow, jcol):
        alseq1 = list()
        alseq2 = list()

        while True:
            if self._pmatrix[irow][jcol] == "diag":
                alseq1.append(self.seq1[jcol - 1])
                alseq2.append(self.seq2[irow - 1])
                if self.seq1[jcol - 1] == self.seq2[irow - 1]:
                    self._identity += 1
                irow -= 1
                jcol -= 1
            elif self._pmatrix[irow][jcol] == "up":
                alseq1.append(Gap(self._gap_character))
                alseq2.append(self.seq2[irow - 1] )
                irow -= 1
            elif self._pmatrix[irow][jcol] == "left":
                alseq1.append(self.seq1[jcol - 1] )
                alseq2.append(Gap(self._gap_character))
                jcol -= 1
            else:
                break
        self._alseq1 = list(reversed(alseq1))
        self._alseq2 = list(reversed(alseq2))
        self._identity = (self._identity / len(self._alseq1)) * 100
        return

    def _check_best_score(self, diagscore, topscore, leftscore, irow, jcol):
        """
        Decides best score for a given cell. Will change depending if using
        Needleman-Wunsch or Smith-Waterman
        """
        return

class ScoreMatrix(object):
    def __init__(self, match, miss, gap):
        self.match = match
        self.miss  = miss
        self.gap   = gap

        try:
            assert(isinstance(self.match, int) is True)
            assert(isinstance(self.miss, int) is True)
            assert(isinstance(self.gap, int) is True)
        except AssertionError:
            raise ValueError("match, miss, and gap should be integers.")

    def __str__(self):
        print_str = "Match:%s Missmatch:%s Gap:%s" % (self.match, self.miss, self.gap)
        return print_str

class Gap(object):
    def __init__(self, character="-"):
        self.character = character

    def __str__(self):
        return str(self.character)

    def __eq__(self, other):
        return isinstance(other, Gap)