 #ifndef UTILS
 #define UTILS

#include <iostream>
#include <string>
#include "index.h"

using namespace std;

string seq_from_fasta(const string& filename);
uint64_t getMemorySelfMaxUsed ();
uint32_t str2numstrand(const string& str);

void delta_encoding(const rh_umap& hmap);
void delta_decoding(const vector<uint32_t> reads);

void vector_compression(rh_umap& hmap);
void vector_decompression(const vector<uint32_t> reads);

double koToBit(double ko);
string intToString(uint64_t num);
uint32_t uniqueKmers(const sparse_vector_u32& index);
void write_values(string output_file, float val1, const string& val2);
kmer hash64shift(kmer key);
uint32_t revhash(uint32_t x);
uint32_t unrevhash(uint32_t x);
uint64_t nuc2int(char c); 
char revCompChar(const char& c);
uint64_t shiftnucl (uint64_t kmerint, char c, uint16_t k);
uint64_t shiftnucl_revcomp(uint64_t val, const char& c, uint16_t k);
string revcomp(const string& s);
uint64_t min_kmer(uint64_t kmer1, uint64_t kmer2);
uint64_t find_canonical(uint64_t kmer_int, uint64_t revComp_int, uint64_t bitvector_size);
uint64_t modulo(uint64_t nb, uint32_t m);
void transform_file(const string& fich1, const string& fich2);
string compare_files(const string& file1, const string& file2);

#endif