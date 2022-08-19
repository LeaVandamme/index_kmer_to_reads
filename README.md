# Indexing K-mers to long reads

## Compilation

git clone https://gitlab.cristal.univ-lille.fr/lvandamm/indexkmertoread
git submodule init
git submodule update

make -j

## Launch 

#### Arguments :

- Fasta file (reads) : path to the fasta file containing reads
- K-mer length
- BitVector Size
- Compressed or not (-c or -nc)

Used to create the index.

#### Options :

- Fasta file (containing 1 kmer to search it in the index)
- Fasta file (containing 1 sequence, to search and count kmers in the index)
- Fasta file (containing a list of sequences, to search and count all the kmers in the index)

<code>
./index [Fasta file (reads)] [Kmer length] [Bitvector Size] [Compressed (-c or -nc)] -k [Fasta file (1 kmer)] -s [Fasta file (1 sequence)] -f [Fasta file (several sequences)]
</code>

<br/>Example to create the index :

<code>
make clean && make -j && ./index input_files/phage/phage_hifi.fa 15 26 -c
</code>

<br/> 
