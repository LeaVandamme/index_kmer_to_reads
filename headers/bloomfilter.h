 #ifndef BLOOMFILTER
 #define BLOOMFILTER

#include <iostream>
#include <vector>

//#include "utils.h"
#include "../sources/bool_iterator.cpp"

using namespace std;

class BloomFilter{


    private:

        uint8_t nb_hash_function;
        uint32_t size;
        vector<bool> filter;
        bool_iterator it_begin,it_end;


    public:

        BloomFilter(uint8_t nb_hash_function, uint32_t size);

        vector<bool>* get_filter();
        bool_iterator* getBegin();
        bool_iterator* getEnd();
        
        void insert(uint32_t kmer);
        bool contains(uint32_t kmer);
        void clear_filter();
        vector<uint32_t> filter_to_array(const vector<bool>& filter);
        bool_iterator begin()const{return it_begin;}  ;
        bool_iterator end()const{return it_end;} ;
};

#endif