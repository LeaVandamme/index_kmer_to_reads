# Indexing K-mers to long reads

## Compilation

git clone --recurse-submodules https://github.com/LeaVandamme/index_kmer_to_reads.git

cd SIMDCompressionAndIntersection/ ; make

cd .. ; cd TurboPFor-Integer-Compression/ ; make

cd .. ; make -j

## Launch 

#### Arguments :

- Fasta file (reads) : path to the fasta file containing reads
- K-mer length
- BitVector Size

<code>
./index [Fasta file (reads)] [Kmer length] [Bitvector Size]
</code>

<br/>Example to create the index :

<code>
make clean && make -j && ./index input_files/phage/phage_1000.fa 15 32
</code>

<br/> 
