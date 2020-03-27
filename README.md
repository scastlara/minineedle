[![Build Status](https://travis-ci.org/scastlara/minineedle.svg?branch=master)](https://travis-ci.org/scastlara/minineedle) [![Coverage Status](https://coveralls.io/repos/github/scastlara/minineedle/badge.svg?branch=master)](https://coveralls.io/github/scastlara/minineedle?branch=master) [![PyPI version](https://badge.fury.io/py/minineedle.svg)](https://badge.fury.io/py/minineedle)

<img width="250" src="https://github.com/scastlara/minineedle/blob/master/minineedle/logo.png"/>

Needleman-Wunsch and Smith-Waterman algorithms in python for any iterable objects.

## Version
v2.1.0

## Algorithms

### Needleman-Wunsch
<img src="https://upload.wikimedia.org/wikipedia/commons/3/3f/Needleman-Wunsch_pairwise_sequence_alignment.png" width="300px">

> The Needleman–Wunsch algorithm is an algorithm used in bioinformatics to align protein or nucleotide sequences. It was one of the first applications of dynamic programming to compare biological sequences. The algorithm was developed by Saul B. Needleman and Christian D. Wunsch and published in 1970. The algorithm essentially divides a large problem (e.g. the full sequence) into a series of smaller problems and uses the solutions to the smaller problems to reconstruct a solution to the larger problem. It is also sometimes referred to as the optimal matching algorithm and the global alignment technique. The Needleman–Wunsch algorithm is still widely used for optimal global alignment, particularly when the quality of the global alignment is of the utmost importance. 
>
> -- From the <cite>[Wikipedia article](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)</cite>

### Smith-Waterman
<img src="https://upload.wikimedia.org/wikipedia/commons/9/92/Smith-Waterman-Algorithm-Example-En.gif" width="300px">

> The Smith–Waterman algorithm performs local sequence alignment; that is, for determining similar regions between two strings of nucleic acid sequences or protein sequences. Instead of looking at the entire sequence, the Smith–Waterman algorithm compares segments of all possible lengths and optimizes the similarity measure. 
>
> -- From the <cite>[Wikipedia article](https://en.wikipedia.org/wiki/Smith–Waterman_algorithm)</cite>


## Usage

```python
import minineedle


# Use miniseq objects
# Load sequences as miniseq FASTA object
import miniseq
fasta = miniseq.FASTA(filename="myfasta.fa")
seq1, seq2 = fasta[0], fasta[1]

# Or use strings, lists, etc
# seq1, seq2 = "ACTG", "ATCTG"
# seq1, seq2 = ["A","C","T","G"], ["A","T","C","T","G"]

# Create the instance
alignment = minineedle.NeedlemanWunsch(seq1, seq2)
# or
# alignment = minineedle.SmithWaterman(seq1, seq2)

# Make the alignment
alignment.align()

# Get the score
alignment.get_score()

# Get the sequences aligned as lists
al1, al2 = alignment.get_aligned_sequences("list")

# Get the sequences as strings
al1, al2 = alignment.get_aligned_sequences("str")

# Change the matrix and run again
alignment.change_matrix(minineedle.ScoreMatrix(match=4, miss=-4, gap=-2))
alignment.align()

# Print the sequences aligned
print(alignment)

# Change gap character
alignment.gap_character = "-gap-"
print(alignment)

# Sort a list of alignments by score
first_al  = NeedlemanWunsch(seq1, seq2)
second_al = NeedlemanWunsch(seq3, seq4)

for align in sorted([first_al, second_al], reverse=True):
    print(align)

```



## Install
```bash
pip install minineedle
```


## Classes

### NeedlemanWunsch
Needleman-Wunsch alignment class. It has the following attributes:
- seq1
- seq2     
- alseq1   
- alseq2
- nmatrix   
- pmatrix   
- smatrix  
- score    
- identity
- gap_character

To create the instance you have to provide two iterable objects with elements that can be compared with "==".

### SmithWaterman
Smith-Waterman alignment class. It has the following attributes:
- seq1
- seq2     
- alseq1   
- alseq2
- nmatrix   
- pmatrix   
- smatrix  
- score    
- identity

To create the instance you have to provide two iterable objects with elements that can be compared with "==".

### ScoreMatrix
With this class you can define your own score matrices. It has three attributes:
- match
- miss
- gap


## Methods
### align()
Performs the alignment.

### get_score()
Returns the score of the alignment. It runs align() if it has not been done yet.

### change_matrix(newmatrix)
Takes a ScoreMatrix object and updates the matrix for the alignment. You still have to run it calling `align()`.

### get identity()
Returns the % of identity (rounded with 2 decimal points).

### get_almatrix()
Return the alignment matrix as a list of lists.
