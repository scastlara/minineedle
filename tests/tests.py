import minineedle

def test_iterable_false():
    '''
    Tests if non-iterable objects are rejected by NeedlemanWunsch.
    '''
    seq1 = "ACTG"
    seq2 = 1
    try:
        minineedle.NeedlemanWunsch(seq1, seq2)
        assert(1 == 0)
    except TypeError:
        assert(1 == 1)
    
def test_iterable_true():
    '''
    Tests if iterable objects are accepted by NeedlemanWunsch.
    '''
    seq1 = "ACTG"
    seq2 = "ACTG"
    try:
        minineedle.NeedlemanWunsch(seq1, seq2)
        assert(1 == 1)
    except Exception:
        assert(1 == 0)

def test_iterable_true_list():
    '''
    Tests if iterable (lists) objects are accepted by NeedlemanWunsch.
    '''
    seq1 = ['A', 'C', 'T', 'G']
    seq2 = ['A', 'C', 'T', 'G']
    try:
        minineedle.NeedlemanWunsch(seq1, seq2)
        assert(1 == 1)
    except Exception:
        assert(1 == 0)

def test_iterable_different():
    '''
    Tests if iterable (lists) objects are accepted by NeedlemanWunsch.
    '''
    seq1 = ['A', 'C', 'T', 'G']
    seq2 = "ACTG"
    try:
        minineedle.NeedlemanWunsch(seq1, seq2)
        assert(1 == 1)
    except Exception:
        assert(1 == 0)

def test_change_matrix_false():
    '''
    Tests if non-ScoreMatrix object is rejected by NeedlemanWunsch's change_matrix
    method. 
    '''
    seq1 = "ACTG"
    seq2 = "ACTG"
    matrix = "NOT A MATRIX"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    
    try:
        needle.change_matrix(matrix)
        assert(1 == 0)
    except ValueError:
        assert(1 == 1)

def test_change_matrix_true():
    '''
    Tests if ScoreMatrix object is accepted by NeedlemanWunsch's change_matrix
    method. 
    '''
    seq1 = "ACTG"
    seq2 = "ACTG"
    matrix = minineedle.ScoreMatrix(1, -1, -1)
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    
    try:
        needle.change_matrix(matrix)
        assert(1 == 1)
    except ValueError:
        assert(1 == 0)

def test_matrix_parameters():
    '''
    Checks if ScoreMatrix object rejects non-integer parameters.
    '''
    try:
        matrix = minineedle.ScoreMatrix(1, 2, "b")
        assert(1 == 0)
    except ValueError:
        assert(1 == 1)

def test_needleman_dynamic_matrix():
    '''
    Checks if alignment score matrix is correct.
    '''
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    expected_matrix = [
        [0,-1,-2,-3,-4,-5,-6,-7],
        [-1,1,0,-1,-2,-3,-4,-5],
        [-2,0,0,1,0,-1,-2,-3],
        [-3,-1,-1,0,2,1,0,-1],
        [-4,-2,-2,-1,1,1,0,-1],
        [-5,-3,-3,-1,0,0,0,-1],
        [-6,-4,-2,-2,-1,-1,1,0],
        [-7,-5,-3,-1,-2,-2,0,0]
        
    ]
    assert(needle._nmatrix == expected_matrix)

def test_needleman_pointer_matrix():
    '''
    Tests if pointer matrix is correct.
    '''
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    expected_matrix = [
        [None, 'left', 'left', 'left', 'left', 'left', 'left', 'left'],
        ['up', 'diag', 'left', 'left', 'left', 'diag', 'left', 'left'],
        ['up', 'up', 'diag', 'diag', 'left', 'left', 'left', 'left'],
        ['up', 'up', 'diag', 'up', 'diag', 'left', 'left', 'left'],
        ['up', 'up', 'diag', 'up', 'diag', 'diag', 'diag', 'diag'],
        ['up', 'up', 'diag', 'diag', 'up', 'diag', 'diag', 'diag'],
        ['up', 'up', 'diag', 'up', 'up', 'diag', 'diag', 'left'],
        ['up', 'up', 'up', 'diag', 'left', 'diag', 'up', 'diag']
    ]
    assert(needle._pmatrix == expected_matrix)

def test_needleman_alignment():
    '''
    Checks if align method works
    '''
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle._alseq1 == ["G", "C", "A", "-", "T", "G", "C", "U"])
    assert(needle._alseq2 == ["G", "-", "A", "T", "T", "A", "C", "A" ])

def test_needleman_alignment_score():
    '''
    Tests if score is correctly computed.
    '''
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle.get_score() == 0)

def test_needleman_alignment_score_list():
    '''
    Tests if score is correctly computed.
    '''
    seq1 = ["G","C","A","T","G","C","U"]
    seq2 = ["G","A","T","T","A","C","A"]
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle.get_score() == 0)

def test_needleman_alignment_score_different():
    '''
    Tests if score is correctly computed.
    '''
    seq1 = ["G","C","A","T","G","C","U"]
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle.get_score() == 0)

def test_needleman_alignment_integers():
    '''
    Tests if algorithm can handle integers
    '''
    seq1 = [1,2,3,5,1]
    seq2 = [1,2,9,9,9,3,5,1]
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle.get_score() == 2)


def test_needleman_alignment_residue_mix():
    '''
    Tests if algorithm can handle different types in list
    '''
    seq1 = [1,"C","A",2,"G","C","U"]
    seq2 = ["G","A","T","T","A","C","A"]
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    try:
        assert(1==1)
    except ValueError:
        assert(1==0)

def test_needleman_alignment_identity():
    '''
    Tests if score is correctly computed.
    '''
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle.get_identity() == 50.0)

def test_needleman_to_string():
    '''
    Tests string formatting of alignment object.
    '''
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(str(needle) == """Alignment of SEQUENCE 1 and SEQUENCE 2:\n\tGCA-TGCU\n\tG-ATTACA\n""")

def test_needleman_comparison_equal():
    '''
    Tests score comparison of equal alignments
    '''
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle1 = minineedle.NeedlemanWunsch(seq1, seq2)
    needle1.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle1.align()
    needle2 = minineedle.NeedlemanWunsch(seq1, seq2)
    needle2.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle2.align()
    assert(needle1 == needle2)

def test_needleman_comparison_lower():
    '''
    Tests score comparison of different alignments
    '''
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    seq3 = "LLLLLLL"
    needle1 = minineedle.NeedlemanWunsch(seq1, seq2)
    needle1.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle1.align()
    needle2 = minineedle.NeedlemanWunsch(seq1, seq3)
    needle2.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle2.align()
    assert(needle1 > needle2)

def test_needleman_get_aligned_squences_list():
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.align()
    seq_a, seq_b = needle.get_aligned_sequences()
    assert(isinstance(seq_a, list))
    assert(isinstance(seq_b, list))

def test_needleman_get_aligned_squences_string():
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.align()
    seq_a, seq_b = needle.get_aligned_sequences("str")
    assert(isinstance(seq_a, str))
    assert(isinstance(seq_b, str))

def test_needleman_get_aligned_squences_wrong():
    seq1 = "GCATGCU"
    seq2 = "GATTACA"
    needle = minineedle.NeedlemanWunsch(seq1, seq2)
    needle.align()
    try:
        seq_a, seq_b = needle.get_aligned_sequences("wrong")
        assert(0 == 1)
    except ValueError:
        assert(1 == 1)
    else:
        assert(0 == 1)

def test_smith_dynamic_matrix():
    '''
    Checks if alignment score matrix is correct.
    '''
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    smith = minineedle.SmithWaterman(seq1, seq2)
    smith.change_matrix(minineedle.ScoreMatrix(3, -3, -2))
    smith.align()
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
        [0, 1, 0, 3, 2, 7, 9, 8, 7]
    ]
    assert(smith._nmatrix == expected_matrix)

def test_smith_alignment():
    '''
    Checks if align method works
    '''
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    smith = minineedle.SmithWaterman(seq1, seq2)
    smith.change_matrix(minineedle.ScoreMatrix(3, -3, -2))
    smith.align()
    assert(smith._alseq1 == ["G", "T", "T", "-", "A", "C"])
    assert(smith._alseq2 == ["G", "T", "T", "G", "A", "C"])

def test_smith_alignment_score():
    '''
    Checks if align method works
    '''
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    smith = minineedle.SmithWaterman(seq1, seq2)
    smith.change_matrix(minineedle.ScoreMatrix(3, -3, -2))
    smith.align()
    assert(smith.get_score() == 13)