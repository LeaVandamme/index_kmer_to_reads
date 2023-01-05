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
typedef bm::bvector<>::rs_index_type bvector_rankselect;
typedef bm::sparse_vector<uint32_t, bm::bvector<> > sparse_vector_u32;



class Index{

    public:
        string filename;
        uint16_t k;
        uint32_t nb_kmers;
        vector<uint8_t> compressed_index;
        sparse_vector_u32* counting_bf;
        bm::bvector<> bf;
        bvector_rankselect* rs_idx;
        vector<uint32_t> vect_pos;
        vector<uint32_t> read_line_pos;
        
        Index(string& filename, uint16_t k);

        void set_nb_kmers(uint32_t nb_kmers);
        void set_compressed_index(vector<uint8_t>& compressed_index);
        void set_counting_bf(sparse_vector_u32* counting_bf);
        void set_bf(bm::bvector<>& bf);
        void set_rs_idx(bvector_rankselect* rs_idx);
        void set_vect_pos(vector<uint32_t>& vect_pos);
        void set_read_line_pos(const vector<uint32_t>& read_line_pos);

        void index_compressed(const string& read_file, uint16_t k, uint16_t bitvector_size, const string& memory_file);
        vector<uint32_t> query_kmer(const string& kmer) const;
        vect_occ_read query_sequence(const string& sequence) const;
        void query_sequence_fp(const string& sequence, const string& file_out, uint16_t threshold);
        vect_occ_read query_fasta(const string& filename) const;
        uint32_t uniqueKmers();

        vector<uint32_t> kmers_array(const string& sequence) const;
        unordered_map<uint32_t, uint16_t> get_possible_reads_threshold(const string& sequence, uint threshold) const;
        bool textual_verif(const string& kmer, const string& read);
        vector<uint32_t> verif_fp(const unordered_map<uint32_t, uint16_t>& reads_to_verify, const vector<uint32_t>& kmers, uint16_t threshold) const;

        void memory_file(const string& output_file, uint16_t k, uint16_t bv_size, size_t uncomp_index_size) const;

};

#endif