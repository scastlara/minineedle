import sys

class Needleman(object):
    '''
    Needleman Wunsch Alignment object. Takes two sequence objects (seq1 and seq2) and aligns them with the method align.
    '''
    def __init__(self, seq1, seq2):
        self.seq1     = seq1
        self.seq2     = seq2
        self.alseq1   = list()
        self.alseq2   = list()
        self.smatrix  = ScoreMatrix(match=1, miss=-1, gap=-1)
        self.__score    = int()
        self.__identity = float()
        self.__nmatrix  = [ [ 0 for x in range(len(self.seq1) + 1)] for x in range(len(self.seq2) +1)]
        self.__pmatrix  = [ [ None for x in range(len(self.seq1) + 1)] for x in range(len(self.seq2) +1)]

        try:
            self.seq1.get_sequence()
            self.seq2.get_sequence()
            self.seq1.get_identifier()
            self.seq2.get_identifier()
        except NoSequenceObject as error:
            sys.stderr.write(str(error))
            sys.exit(1)


    def __str__(self):
        if not self.alseq1:
            self.align()
        return str("Alignment of %s and %s:\n\t%s\n\t%s\n" % (self.seq1.get_identifier(),
                                                              self.seq2.get_identifier(),
                                                              "".join(self.alseq1),
                                                              "".join(self.alseq2)
                                                              ))

    def __lt__(self, other):
        return self.get_score() < other.get_score()

    def __eq__(self, other):
        return self.get_score() == other.get_score()

    def get_score(self):
        if not self.alseq1:
            self.align()
        return self.__score

    def change_matrix(self, newmatrix):
        if isinstance(newmatrix, ScoreMatrix):
            self.smatrix = newmatrix
        else:
            raise ValueError("New matrix should be a ScoreMatrix object.")

    def align(self):
        # Init results and pointers matrix
        for row in self.__pmatrix:
            row[0] = "up"

        for jcol in range(0, len(self.__pmatrix[0])):
            self.__pmatrix[0][jcol] = "left"

        self.__pmatrix[0][0] = None

        # Add Gap penalties
        for i in range(1, len(self.seq1) + 1 ):
            self.__nmatrix[0][i] = self.__nmatrix[0][i-1] + self.smatrix.gap

        for j in range(1, len(self.seq2) +1 ):
            self.__nmatrix[j][0] = self.__nmatrix[j-1][0] + self.smatrix.gap

        # Iterate through the sequences and nmatrix
        self.__fill_matrices()

        # Get score and alignment
        imax = len(self.__nmatrix) - 1
        jmax = len(self.__nmatrix[0]) -1
        self.__score = self.__nmatrix[imax][jmax]
        self.__get_alignment(imax, jmax)

    def get_almatrix(self):
        """
        Returns the alignment matrix (list of lists)
        """
        if not self.alseq1:
            self.align()
        return self.__nmatrix

    def get_identity(self):
        """
        Returns the % of identity of the alignment
        """
        if not self.alseq1:
            self.align()
        return round(self.__identity, 2) # Two decimal points

    def __fill_matrices(self):
        for irow in range(0, len(self.seq2)):
            for jcol in range(0, len(self.seq1)):
                # Scores
                topscore  = self.__nmatrix[irow][jcol + 1] + self.smatrix.gap
                leftscore = self.__nmatrix[irow + 1][jcol] + self.smatrix.gap
                diagscore = self.__nmatrix[irow][jcol]
                if self.seq1.get_sequence()[jcol] == self.seq2.get_sequence()[irow]:
                    diagscore += self.smatrix.match
                else:
                    diagscore += self.smatrix.miss

                # Get best score
                self.__check_best(diagscore, topscore, leftscore, irow, jcol)
        return

    def __get_alignment(self, irow, jcol):
        alseq1 = list()
        alseq2 = list()

        while True:
            if self.__pmatrix[irow][jcol] == "diag":
                alseq1.append(self.seq1.get_sequence()[jcol - 1])
                alseq2.append(self.seq2.get_sequence()[irow - 1])
                if self.seq1.get_sequence()[jcol - 1] == self.seq2.get_sequence()[irow - 1]:
                    self.__identity += 1
                irow -= 1
                jcol -= 1
            elif self.__pmatrix[irow][jcol] == "up":
                alseq1.append("-")
                alseq2.append(self.seq2.get_sequence()[irow - 1] )
                irow -= 1
            elif self.__pmatrix[irow][jcol] == "left":
                alseq1.append(self.seq1.get_sequence()[jcol - 1] )
                alseq2.append("-")
                jcol -= 1
            else:
                break
        self.alseq1 = list(reversed(alseq1))
        self.alseq2 = list(reversed(alseq2))
        self.__identity = (self.__identity / len(self.alseq1)) * 100
        return

    def __check_best(self, diagscore, topscore, leftscore, irow, jcol):
        best_pointer = str()
        best_score   = int()
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

        self.__pmatrix[irow + 1][jcol +1] = best_pointer
        self.__nmatrix[irow + 1][jcol +1] = best_score
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


class NoSequenceObject(Exception):
    '''
    Exception to handle no sequence objects passed to Needleman
    '''
    def __str__(self):
        return("ERROR: Class Needleman takes two sequence objects")
