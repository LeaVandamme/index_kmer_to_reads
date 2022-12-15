 #ifndef UTILS
 #define UTILS

#include <iostream>
#include <string>
#include "index.h"

using namespace std;

string seq_from_fasta(const string& filename);
uint64_t getMemorySelfMaxUsed ();
uint32_t str2numstrand(const string& str);
string num2strstrand(uint32_t integer);

void delta_decoding(const vector<uint32_t> reads);

double koToBit(double ko);
string intToString(uint64_t num);
uint32_t uniqueKmers(const sparse_vector_u32& index);
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

vector<uint32_t> read_line_position(const string& filename);
string get_read_sequence(const string& filename, uint32_t read_id);
void get_reads(const string& read_file_out, const vector<uint32_t>& id_reads);


#endif