#include <iostream>
#include <string>
#include <fstream>
#include <sys/resource.h>
#include <algorithm>
#include <bitset>
#include "../include/fastDelta.h"
#include "../headers/index.h"
#include "../headers/bloomfilter.h"
#include "../SIMDCompressionAndIntersection/include/codecfactory.h"
#include "../SIMDCompressionAndIntersection/include/intersection.h"




using namespace SIMDCompressionLib;
using namespace std;




string seq_from_fasta(const string& filename){

    ifstream fichier(filename, ios::in);
    string seq;
    if(fichier){   
        string ligne;
        while(getline(fichier,ligne)){
            if (ligne[0] != '>'){
                seq = ligne;
            }
        }
        fichier.close();
    }
    
    else cerr << "Error opening the file." << endl;
    return seq;
}




uint64_t getMemorySelfMaxUsed (){
	uint64_t result = 0;
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage)==0){  
        result = usage.ru_maxrss;  
    }
	return result;
}




uint32_t str2numstrand(const string& str) {
  uint32_t res(0);
  for(uint i(0);i<str.size();i++) {
    res<<=2;
    switch (str[i]){
      case 'A':res+=0;break;
      case 'C':res+=1;break;
      case 'G':res+=2;break;
      case 'T':res+=3;break;

      case 'a':res+=0;break;
      case 'c':res+=1;break;
      case 'g':res+=2;break;
      case 't':res+=3;break;
      default: return 0 ;break;
    }
  }
  return (res);
}



// void delta_encoding(const rh_umap& index){
//   for (auto const &pair: index) {
//     int N = pair.second.size();
//     compute_deltas_inplace((uint32_t *) pair.second.data(),N,0);
//   }
// }



void delta_decoding(const vector<uint32_t> reads){
  int N = reads.size();
  compute_prefix_sum_inplace((uint32_t *) reads.data(),N,0);
}



// void vector_compression(rh_umap& index){
//   IntegerCODEC &codec = *CODECFactory::getFromName("s4-bp128-d4");
//   for (auto &pair: index) {
//     if(pair.second->size() > 5 or true){
//        int N = index[pair.first]->size();
//         vector<uint32_t> compressed_output=vector<uint32_t>(N+8);
//         size_t compressedsize = compressed_output.size();
//         codec.encodeArray((uint32_t *) pair.second->data(), N, compressed_output.data(), compressedsize);
//         delete(index[pair.first]);
//         compressed_output.resize(compressedsize);
//         index[pair.first]  = new vector<uint32_t>;
//         for(uint i(0);i<compressedsize;++i){
//           (index[pair.first] )->push_back(compressed_output[i]);
//         }
//           vector<uint32_t> mydataback(N);
//           size_t recoveredsize = mydataback.size();
//           //
//           codec.decodeArray(compressed_output.data(), compressedsize,
//                             mydataback.data(), recoveredsize);
//             if (mydataback != *pair.second)
//               throw runtime_error("bug!"); 
//          //*light_vector=compressed_output;
//          cout << static_cast<double>(compressed_output.size()) << " " << static_cast<double>(pair.second->size()) << endl;
//          cout<<N<<" "<<index[pair.first] ->size()<<" "<<index[pair.first] ->capacity()<<endl;
//     }
//   }
// }




void vector_decompression(const vector<uint32_t> reads){
  IntegerCODEC &codec = *CODECFactory::getFromName("s4-bp128-d4");
  int N = reads.size();
  vector<uint32_t> mydataback(N);
  size_t recoveredsize = mydataback.size();
  codec.decodeArray(reads.data(), N, mydataback.data(), recoveredsize);
  //mydataback.shrink_to_fit();
}




double koToBit(double ko){
  return ko * 8000;
}




string intToString(uint64_t num) {
  string s = std::to_string(num);

  int n = s.length() - 3;
  int end = (num >= 0) ? 0 : 1;
  while (n > end) {
    s.insert(n, ",");
    n -= 3;
  }
  return s;
}




void write_values(string output_file, float val1, const string& val2){
  ofstream outp(output_file, ios::app);
    if (outp){
        outp << val1 << " " << val2 << endl;
    }
  outp.close();
}




uint32_t revhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x2c1b3c6d;
	x = ((x >> 16) ^ x) * 0x297a2d39;
	x = ((x >> 16) ^ x);
	return x;
}




uint32_t unrevhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ((x >> 16) ^ x) * 0x64ea2d65;
	x = ((x >> 16) ^ x);
	return x;
}



uint64_t shiftnucl(uint64_t val, char c, uint16_t k) {
	val = val << 2;
  val = val + str2numstrand(string(1,c));
  uint mask = pow(2,(k*2))-1;
  return val & mask;
}


char revCompChar(const char& c){
  switch (c){
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
    }
    return 'A';
}



string revcomp(const string& s){
  string res(s.size(), 0);
  for(int i((int)s.length() - 1); i >= 0; i--){
    res[s.size()-1-i] = revCompChar(s[i]);
  }
  return res;
}




uint64_t shiftnucl_revcomp(uint64_t val, const char& c, uint16_t k) {
	val = val >> 2;
  return val + (str2numstrand(string(1,revCompChar(c))) << ((k-1)*2));
}




uint64_t min_kmer(uint64_t kmer1, uint64_t kmer2){
  if(kmer1<kmer2){
    return kmer1;
  }
  return kmer2;
}



uint64_t modulo(uint64_t nb, uint32_t m){
  return nb & (m-1);
}



uint64_t find_canonical(uint64_t kmer_int, uint64_t revComp_int, uint64_t bitvector_size){
  return modulo(revhash(min_kmer(kmer_int, revComp_int)), bitvector_size);
}



void transform_file(const string& fich1, const string& fich2){
  
  ifstream fichier1(fich1, ios::in);
  ofstream outp(fich2);
  if(outp){
    string ligne;
    while(!fichier1.eof()){
      getline(fichier1,ligne);
      ligne.erase(0, ligne.find(":")+1);
      ligne.erase(std::remove_if(ligne.begin(), ligne.end(), ::isspace), ligne.end());
      outp << ligne << "\n";
    }
  }
  fichier1.close();
  outp.close();
}



string compare_files(const string& file1, const string& file2){
  string res = "";
  int cpt = 1;
  ifstream fichier1(file1, ios::in);
  ifstream fichier2(file2, ios::in);

  if(fichier1 && fichier2){
    string ligne1, ligne2;
    while(!fichier1.eof() && !fichier2.eof()){
      getline(fichier1,ligne1);
      getline(fichier2,ligne2);

      if(ligne1 != ligne2){
        res += std::to_string(cpt);
        res += " ";
      }
      cpt++;
    }
  }
  return res;
}



