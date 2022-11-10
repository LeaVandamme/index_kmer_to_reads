 #ifndef INDEX
 #define INDEX

#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include "../include/robin_hood.h"
#include "bloomfilter.h"
#include "../BitMagic/src/bm.h"
#include "../BitMagic/src/bmundef.h"

using namespace std;

struct occ_in_read {
            uint32_t total_count;
            uint32_t read_id;
        };



struct rank_kmer {
            uint64_t rank;
            uint32_t nb_apparition;
        };

typedef robin_hood::unordered_map<uint32_t, vector<uint32_t>*> rh_umap;
typedef vector<occ_in_read> vect_occ_read;
typedef uint64_t kmer;
typedef vector<vector<uint32_t>*>  index_vector;
typedef unique_ptr<bm::bvector<>::rs_index_type> bvector_rankselect;
typedef bm::sparse_vector<uint32_t, bm::bvector<> > sparse_vector_u32;



class Index{

    private:
        string filename;
        uint16_t k;
        uint32_t nb_kmers;
        
    public:
        Index(string& filename, uint16_t k);

        void set_nb_kmers(uint32_t nb_kmers);
        uint32_t get_nb_kmers();

        vector<uint32_t> index_fasta_rankselect(const string& read_file, uint16_t k, uint16_t bitvector_size);
        vector<uint8_t> index_fasta_rankselect_compressed(const string& read_file, uint16_t k, uint16_t bitvector_size);
        vect_occ_read query_sequence_rankselect(const vector<uint8_t>& index, const string& sequence);
        vector<uint32_t> query_kmer_rankselect(const vector<uint8_t>& index_vect, const string& kmer);
        vect_occ_read query_fasta_rankselect(const vector<uint8_t>& index_vect, const string& filename);
        void write_index_rankselect(const vector<uint32_t>& read_vector, const string& output_file);
        uint32_t uniqueKmers();
};

#endif