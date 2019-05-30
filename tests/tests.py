import miniseq
import minineedle

def test_sequence_false():
    '''
    Tests if non-sequence objects are rejected by Needleman.
    '''
    seq1 = miniseq.DNASequence("seq1", "ACTG")
    seq2 = "NOT A SEQUENCE"
    try:
        minineedle.Needleman(seq1, seq2)
        assert(1 == 0)
    except Exception:
        assert(1 == 1)
    
def test_sequence_true():
    '''
    Tests if sequence objects are accepted by Needleman.
    '''
    seq1 = miniseq.DNASequence("seq1", "ACTG")
    seq2 = miniseq.DNASequence("seq2", "ACTG")
    try:
        minineedle.Needleman(seq1, seq2)
        assert(1 == 1)
    except Exception:
        assert(1 == 0)

def test_change_matrix_false():
    '''
    Tests if non-ScoreMatrix object is rejected by Needleman's change_matrix
    method. 
    '''
    seq1 = miniseq.DNASequence("seq1", "ACTG")
    seq2 = miniseq.DNASequence("seq2", "ACTG")
    matrix = "NOT A MATRIX"
    needle = minineedle.Needleman(seq1, seq2)
    
    try:
        needle.change_matrix(matrix)
        assert(1 == 0)
    except ValueError:
        assert(1 == 1)

def test_change_matrix_true():
    '''
    Tests if ScoreMatrix object is accepted by Needleman's change_matrix
    method. 
    '''
    seq1 = miniseq.DNASequence("seq1", "ACTG")
    seq2 = miniseq.DNASequence("seq2", "ACTG")
    matrix = minineedle.ScoreMatrix(1, -1, -1)
    needle = minineedle.Needleman(seq1, seq2)
    
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

def test_dynamic_matrix():
    '''
    Checks if alignment score matrix is correct.
    '''
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
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
    assert(needle._Needleman__nmatrix == expected_matrix)

def test_pointer_matrix():
    '''
    Tests if pointer matrix is correct.
    '''
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
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
    assert(needle._Needleman__pmatrix == expected_matrix)


def test_alignment():
    '''
    Checks if align method works
    '''
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle.alseq1 == ["G", "C", "A", "-", "T", "G", "C", "U"])
    assert(needle.alseq2 == ["G", "-", "A", "T", "T", "A", "C", "A" ])

def test_alignment_score():
    '''
    Tests if score is correctly computed.
    '''
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle.get_score() == 0)

def test_alignment_identity():
    '''
    Tests if score is correctly computed.
    '''
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(needle.get_identity() == 50.0)

def test_to_string():
    '''
    Tests string formatting of alignment object.
    '''
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
    needle.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle.align()
    assert(str(needle) == """Alignment of seq1 and seq2:\n\tGCA-TGCU\n\tG-ATTACA\n""")

def test_comparison_equal():
    '''
    Tests score comparison of equal alignments
    '''
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle1 = minineedle.Needleman(seq1, seq2)
    needle1.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle1.align()
    needle2 = minineedle.Needleman(seq1, seq2)
    needle2.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle2.align()
    assert(needle1 == needle2)

def test_comparison_lower():
    '''
    Tests score comparison of different alignments
    '''
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    seq3 = miniseq.Sequence("seq3", "LLLLLLL")
    needle1 = minineedle.Needleman(seq1, seq2)
    needle1.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle1.align()
    needle2 = minineedle.Needleman(seq1, seq3)
    needle2.change_matrix(minineedle.ScoreMatrix(1, -1, -1))
    needle2.align()
    assert(needle1 > needle2)

def test_get_aligned_squences_list():
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
    needle.align()
    seq_a, seq_b = needle.get_aligned_sequences()
    assert(isinstance(seq_a, list))
    assert(isinstance(seq_b, list))

def test_get_aligned_squences_string():
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
    needle.align()
    seq_a, seq_b = needle.get_aligned_sequences("str")
    assert(isinstance(seq_a, str))
    assert(isinstance(seq_b, str))

def test_get_aligned_squences_wrong():
    seq1 = miniseq.Sequence("seq1", "GCATGCU")
    seq2 = miniseq.Sequence("seq2", "GATTACA")
    needle = minineedle.Needleman(seq1, seq2)
    needle.align()
    try:
        seq_a, seq_b = needle.get_aligned_sequences("wrong")
        assert(0 == 1)
    except ValueError:
        assert(1 == 1)
    else:
        assert(0 == 1)
    