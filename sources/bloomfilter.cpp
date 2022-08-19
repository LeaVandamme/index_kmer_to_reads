#include <iostream>
#include <vector>
#include <climits>
#include <assert.h>
#include "bool_iterator.cpp"
#include "../headers/bloomfilter.h"
#include "../headers/utils.h"

using namespace std;

BloomFilter::BloomFilter(uint8_t i_nb_hash_function, uint32_t i_size){
   nb_hash_function = i_nb_hash_function;
   size = i_size;
   filter.resize(size,false);
   it_begin=bool_iterator(&filter);
   it_end=bool_iterator(true);
}
  

vector<bool>* BloomFilter::get_filter(){
   return &(filter);
}



bool_iterator* BloomFilter::getBegin(){
   return &it_begin;
}



bool_iterator* BloomFilter::getEnd(){
   return &it_end;
}



// FUNCTIONS



void BloomFilter::insert(uint32_t kmer){
   size_t hashed = revhash(kmer);
   filter[hashed] = true;
}



bool BloomFilter::contains(uint32_t kmer){
   size_t hashed = revhash(kmer);
   return filter[hashed];
}



vector<uint32_t> BloomFilter::filter_to_array(const vector<bool>& filter){
   vector<uint32_t> res;
   int j = 0;
   for(uint i=0; i<filter.size(); i++){
      if(filter[i] == true){
         res.push_back(i);
         j++;
      }
   }
   return res;
}