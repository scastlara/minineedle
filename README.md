<img width="250" src="https://github.com/scastlara/minineedle/blob/master/minineedle/logo.png"/>

Needleman Wunsch algorithm in python using [miniseq](https://github.com/scastlara/miniseq) objects

## Version
v1.0.0

## Algorithm
<img src="https://upload.wikimedia.org/wikipedia/commons/3/3f/Needleman-Wunsch_pairwise_sequence_alignment.png" width="300px">

> The Needleman–Wunsch algorithm is an algorithm used in bioinformatics to align protein or nucleotide sequences. It was one of the first applications of dynamic programming to compare biological sequences. The algorithm was developed by Saul B. Needleman and Christian D. Wunsch and published in 1970. The algorithm essentially divides a large problem (e.g. the full sequence) into a series of smaller problems and uses the solutions to the smaller problems to reconstruct a solution to the larger problem. It is also sometimes referred to as the optimal matching algorithm and the global alignment technique. The Needleman–Wunsch algorithm is still widely used for optimal global alignment, particularly when the quality of the global alignment is of the utmost importance. 
>
> -- From the <cite>[Wikipedia article](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)</cite>

## Usage

```python
import miniseq
import minineedle

# Load sequences as miniseq FASTA object
fasta = miniseq.FASTA(filename="myfasta.fa")
seq1, seq2 = fasta[0], fasta[1]

# Create the instance
alignment = minineedle.Needleman(seq1, seq2)

# Make the alignment
alignment.align()

# Get the score
alignment.get_score()

# Get the sequences aligned as lists
al1 = alignment.alseq1
al2 = alignment.alseq2

# Change the matrix and run again
alignment.change_matrix(minineedle.ScoreMatrix(match=4, miss=-4, gap=-2))
alignment.align()

# Print the sequences aligned
print(alignment)

# Sort a list of alignments by score
first_al  = Needleman(seq1, seq2)
second_al = Needleman(seq3, seq4)

for align in sorted([first_al, second_al], reverse=True):
    print(align)

```

## Install
```bash
git clone https://github.com/scastlara/minineedle.git
sudo python3 setup.py install
```


## Classes

### Needleman
Alignment class. It has the following attributes:
- seq1
- seq2     
- alseq1   
- alseq2
- nmatrix   
- pmatrix   
- smatrix  
- score    
- identity

To create the instance you have to provide seq1 and seq2 as two strings.

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
