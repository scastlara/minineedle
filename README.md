<img width="250" src="https://github.com/scastlara/minineedle/blob/master/minineedle/logo.png"/>
Needleman Wunsch algorithm in python using [miniseq](https://github.com/scastlara/miniseq) objects

## Version
v0.1.0

## Usage

```python
import miniseq
import minineedle

# Load sequences as miniseq FASTA object
fasta = miniseq.FASTA(filename="myfasta.fa")
seq1, seq2 = fasta.sequences

fasta
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
