#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <math.h>
#include <thread>
#include <unistd.h>
#include <bits/stdc++.h>
#include "../BitMagic/src/bm.h"
#include "../BitMagic/src/bmundef.h"
#include "../BitMagic/src/bmsparsevec.h"
#include "../headers/utils.h"
#include "../include/robin_hood.h"
#include "../TurboPFor-Integer-Compression/vp4.h"



using namespace std;
using namespace chrono;

using std::find;


Index::Index(string& fastaFilename, uint16_t kmerLength){
    filename = fastaFilename;
    k = kmerLength;
}

void Index::set_nb_kmers(uint32_t nb_kmers){
    this-> nb_kmers = nb_kmers;
}
void Index::set_compressed_index(vector<uint8_t>& compressed_index){
    this-> compressed_index = compressed_index;
}
void Index::set_counting_bf(sparse_vector_u32* counting_bf){
    this-> counting_bf = counting_bf;
}
void Index::set_bf(bm::bvector<>& bf){
    this-> bf = bf;
}
void Index::set_rs_idx(bvector_rankselect* rs_idx){
    this-> rs_idx = rs_idx;
}
void Index::set_vect_pos(vector<uint32_t>& vect_pos){
    this-> vect_pos = vect_pos;
}
void Index::set_read_line_pos(const vector<uint32_t>& read_line_pos){
    this-> read_line_pos = read_line_pos;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CREATION



void Index::index_compressed(const string& read_file, uint16_t k, uint16_t bitvector_size, const string& memory_output){
    bm::bvector<> bitvector;
    bm::bvector<>::rs_index_type* rs_idx(new bm::bvector<>::rs_index_type());

    vector<uint32_t> vect_pos = {0};
    sparse_vector_u32* occ_of_kmer = new sparse_vector_u32; // Calculer sa taille en bits
    bitvector.resize(pow(2,bitvector_size));
    occ_of_kmer->resize(pow(2,bitvector_size));
    uint32_t bvsize = bitvector.size();
    vector<uint32_t> read_ids;
    uint32_t num_kmers = 1;
    uint32_t kmer_int, revComp_int, kmer_hash, revComp_hash;
    auto start_bitvector = high_resolution_clock::now();
    auto start_total = high_resolution_clock::now();
    using value_type = bm::sparse_vector<int, bm::bvector<> >::value_type;
    ifstream fichier(read_file, ios::in);
    if(fichier){

        vector<uint32_t> tmp_vect_kmers; // BITVECTOR

        // BLOOM FILTER
        
        while(!fichier.eof()){
            string ligne, kmer; 
            getline(fichier,ligne);
            bool is_first = true;
            if (ligne[0] != '>' && !ligne.empty()){
                uint32_t i = k;
                while(i<ligne.length()){
                    if(is_first){
                        kmer = ligne.substr(0, k);
                        kmer_int = str2numstrand(kmer);
                        revComp_int = str2numstrand(revcomp(kmer));
                        kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                        bitvector.set(kmer_hash);
                        is_first = false;
                    }
                    else{
                        kmer_int = shiftnucl(kmer_int, ligne[i], k);
                        revComp_int = shiftnucl_revcomp(revComp_int, ligne[i], k);
                        kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                        i++;            
                        bitvector.set(kmer_hash);
                        num_kmers++;
                    }

                    // IF KMER NOT ALREADY SEEN FOR THIS READ
                    if (find(tmp_vect_kmers.begin(), tmp_vect_kmers.end(), kmer_int) == tmp_vect_kmers.end()) {
                        occ_of_kmer->inc(kmer_hash);
                        tmp_vect_kmers.push_back(kmer_int);
                    }
                }   
            }
            tmp_vect_kmers.clear();
        }
        set_nb_kmers(num_kmers);
        bitvector.build_rs_index(rs_idx);
        vector<uint32_t>().swap(tmp_vect_kmers);
    }
    else{
        cerr << "Error opening the file." << endl;
    }

    fichier.close();
    auto end_bitvector = high_resolution_clock::now();
    auto querying_bitvector = duration_cast<std::chrono::seconds>(end_bitvector- start_bitvector);
    cout << "ok bitvector : " << querying_bitvector.count() << " seconds."  << endl;

    // VECTOR WITH START POSITIONS FOR READS

    uint64_t total_count = 0;

    for(auto it = occ_of_kmer->begin(); it != occ_of_kmer->end(); it++){
		if(*it != 0){
            total_count += *it;
            vect_pos.push_back(total_count);
        }
	}

    cout << "ok vectpos" << endl;
    read_ids.resize(total_count);

    // INDEX
        
    auto start_index = high_resolution_clock::now();
    uint32_t num_read_index = 0;
    ifstream fichier_index(read_file, ios::in);
    if(fichier_index){
        vector<uint32_t> tmp_vect_kmers;
        string ligne_index, kmer_index;
            while(!fichier_index.eof()){
                bool is_first = true;
                getline(fichier_index,ligne_index);
                if (ligne_index[0] != '>' && !ligne_index.empty()){
                    uint32_t i = k;

                    while(i<ligne_index.length()){
                        if(is_first){
                            kmer_index = ligne_index.substr(0, k);
                            kmer_int = str2numstrand(kmer_index);
                            revComp_int = str2numstrand(revcomp(kmer_index));
                            kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                            is_first = false;
                        }
                        else{
                            kmer_int = shiftnucl(kmer_int, ligne_index[i], k);
                            revComp_int = shiftnucl_revcomp(revComp_int, ligne_index[i], k);
                            kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                            i++;
                        }
                        if (find(tmp_vect_kmers.begin(), tmp_vect_kmers.end(), kmer_int) == tmp_vect_kmers.end()) {
                            auto rank_hash = bitvector.rank(kmer_hash, *rs_idx);
                            uint32_t start_pos = vect_pos[rank_hash-1];
                            uint32_t nb_apparition = vect_pos[rank_hash] - vect_pos[rank_hash-1];

                            (read_ids)[start_pos]++;

                            uint32_t indice_to_insert = start_pos + (read_ids)[start_pos];
                            if(indice_to_insert < (start_pos + nb_apparition)){
                                (read_ids)[indice_to_insert] = num_read_index;
                            }
                            else{
                                (read_ids)[start_pos] = num_read_index;
                            }

                            tmp_vect_kmers.push_back(kmer_int);
                        }
                    }
                    num_read_index ++;
                }
            tmp_vect_kmers.clear();
        }
        vector<uint32_t>().swap(tmp_vect_kmers);
    }
    fichier_index.close();
    cout << "ok index non compressé" << endl;

    // COMPRESSION
    
    vector<uint32_t>::const_iterator first, last;
    vector<unsigned char> compressed_vector;
    uint scomp;
    uint tmp_indice_vect_pos = 0;
    vector<uint8_t> compressed_id;
    vect_pos.insert(vect_pos.begin(), 0);
    vector<uint32_t>::const_iterator begin = read_ids.begin();
    for(uint32_t i = 2; i<vect_pos.size(); i++){
        first = begin + vect_pos[i-1];
        last = begin + vect_pos[i];
        
        vector<uint32_t> sub_index(first, last);
        sort(sub_index.begin(), sub_index.end());

        // compression du bloc
        compressed_vector = vector<unsigned char>(sub_index.size()*16 + 100); 
        scomp=p4ndenc32(sub_index.data(), sub_index.size() , compressed_vector.data());
        compressed_vector.resize(scomp);
        // ajout de chaque élément compressé dans l'index
        for(auto elt : compressed_vector){
            compressed_id.push_back(elt);
        }

        // modif du vecteur avec les positions d'insertion des numéros de read pour chaque bloc
        tmp_indice_vect_pos += scomp;
        vect_pos[i-1] = tmp_indice_vect_pos;
    }
    cout << "ok index compressé" << endl;

    auto end_index = high_resolution_clock::now();
    auto end_total = high_resolution_clock::now();
    auto querying_index = duration_cast<std::chrono::seconds>(end_index - start_index);
    auto total = duration_cast<std::chrono::seconds>(end_total - start_total);
    cout << "ok index : " << querying_index.count() << " seconds."  << endl;
    cout << "Total time : " << total.count() << " seconds."  << endl;
    bm::bvector<>::statistics st;
    bitvector.calc_stat(&st);
    cout << "Taille du Bloom filter : " << st.memory_used << endl;
    cout << "Taille du vecteur positions : " << vect_pos.size() * 4 <<endl;
    cout << "Taille du vecteur identifiants non compressé : " << read_ids.size() * 4<< endl;
    cout << "Taille du vecteur identifiants compressé : " << compressed_id.size() << endl;

    this->set_counting_bf(occ_of_kmer);
    this->set_compressed_index(compressed_id);
    this->set_vect_pos(vect_pos);
    this->set_rs_idx(rs_idx);
    this->set_bf(bitvector);
    vector<uint32_t> read_line_pos = read_line_position(read_file);
    this->set_read_line_pos(read_line_pos);

    memory_file(memory_output, k, bitvector_size, read_ids.size()*4);
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// QUERY




vector<uint32_t> Index::query_kmer(const string& kmer) const {
    bm::bvector<> bitvector = this->bf;
    bvector_rankselect* rs_idx = this->rs_idx;
    vector<uint32_t> vect_pos = this->vect_pos;
    uint32_t kmer_int = str2numstrand(kmer);
    uint32_t revComp_int = str2numstrand(revcomp(kmer));
    uint64_t kmer_hash = find_canonical(kmer_int, revComp_int, bitvector.size());
    uint length = 0;
    uint start = 0;
    if (bitvector[kmer_hash] == 1){ // if the kmer exists
    
        auto rank_hash = bitvector.rank(kmer_hash, *rs_idx);

        vector<uint8_t>::const_iterator first, last;
        vector<uint8_t>::const_iterator begin = this->compressed_index.begin();
        
        first = begin + vect_pos[rank_hash-1];
        last = begin + vect_pos[rank_hash];
        vector<unsigned char> to_decompress(first, last);
        
        vector<uint32_t> uncompressed_vector;
        uint scomp2;

        uint32_t size_uncompressed_index = this->counting_bf->at(kmer_hash);

        uncompressed_vector = vector<uint32_t>(size_uncompressed_index*2+100);
        scomp2=p4nddec32(to_decompress.data(), size_uncompressed_index, uncompressed_vector.data());

        uncompressed_vector.resize(size_uncompressed_index);
        return uncompressed_vector;
    }
    else{
        cout << "K-mer not found." << '\n';
        vector<uint32_t> v;
        return v;
    }
}



vect_occ_read Index::query_sequence(const string& sequence) const{

    vect_occ_read res;
    unordered_map<uint32_t, uint32_t> tmp_reads;
    bool is_first = true;
    for (long unsigned int i=0; i<=sequence.length()-this->k; i++){
        string kmer = sequence.substr(i, this->k);
        vector<uint32_t> reads_of_kmer = query_kmer(kmer);
        for (auto& idRead : reads_of_kmer){
            if(tmp_reads.find(idRead) != tmp_reads.end()){
                tmp_reads[idRead] ++;
            }
            else{
                tmp_reads[idRead] = 1;
            }
        }
    }

    for (auto const &pair: tmp_reads){
        occ_in_read structRead;
        structRead.read_id = pair.first;
        structRead.total_count = pair.second;
        res.push_back(structRead);
    }

    return res;
}



void Index::query_sequence_fp(const string& sequence, const string& file_out, uint16_t threshold){
    vector<uint32_t> kmers = kmers_array(sequence);
    unordered_map<uint32_t, uint16_t> poss_reads = get_possible_reads_threshold(sequence, threshold);
    vector<uint32_t> res = verif_fp(poss_reads, kmers, threshold);
    get_reads(filename, read_line_pos, file_out, res);
}




vect_occ_read Index::query_fasta(const string& filename) const{
    
    vect_occ_read res;
    unordered_map<uint32_t, uint32_t> tmp_reads;

    ifstream fichier(filename, ios::in);

    if(fichier){   
        string ligne;
        while( !fichier.eof()){
            getline(fichier,ligne);
            if (ligne[0] != '>' && !ligne.empty()){
                vect_occ_read occ_read = query_sequence(ligne);
                for (auto& struct_read : occ_read){
                    int read_id = struct_read.read_id;
                    int total_count = struct_read.total_count;
                    if(tmp_reads.find(read_id) != tmp_reads.end()){
                        tmp_reads[read_id] += total_count;
                    }
                    else{
                        tmp_reads[read_id] = total_count;
                    }
                }
            }
        }
        for (auto const &pair: tmp_reads){
            occ_in_read structRead;
            structRead.read_id = pair.first;
            structRead.total_count = pair.second;
            res.push_back(structRead);
        }
        fichier.close();
        return res;
    }
    else{
        cerr << "Error opening the file (query_fasta)." << endl;
        return res;
    }
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OTHER



uint32_t Index::uniqueKmers(){
    vector<uint32_t> vect_pos = this->vect_pos;
    uint32_t nb_unique = 0;
    for (uint i=0; i<vect_pos.size()-1; i++) {
        if((vect_pos[i+1] - vect_pos[i]) == 1){
            nb_unique ++;
        }
    }
    return nb_unique;
}


vector<uint32_t> Index::kmers_array(const string& sequence) const{
    vector<uint32_t> kmers;
    string curr_kmer;
    uint32_t kmer_int;
    for (long unsigned int i=0; i<=sequence.length()-this->k+1; i++){
        curr_kmer = sequence.substr(i, this->k);
        kmer_int = str2numstrand(curr_kmer);
        kmers.push_back(kmer_int);
    }
    sort(kmers.begin(), kmers.end());
    return kmers;
}


unordered_map<uint32_t, uint16_t> Index::get_possible_reads_threshold(const string& sequence, uint threshold) const{
    unordered_map<uint32_t, uint16_t> id_to_count = {};

    string curr_kmer;
    vector<uint32_t> curr_ids_read;
    uint i = 0;

    // INSERTION HASHMAP DU NOMBRE D'APPARITION POUR CHAQUE READ
    for (long unsigned int i=0; i<=sequence.length()-this->k; i++){

        curr_kmer = sequence.substr(i, this->k);

        // INSERTION DANS ARRAY DES KMERS (INT)
        curr_ids_read = query_kmer(curr_kmer);
        
        for(auto id_read : curr_ids_read){
            if(id_to_count.find(id_read) == id_to_count.end()){
                id_to_count[id_read] = 1;
            }
            else{
                for (auto i = id_to_count.begin(); i != id_to_count.end(); i++){
                    if(i->first == id_read){
                        i->second++; // UN MEME READ PEUT ETRE INCRÉMENTÉ PLUSIEURS FOIS POUR LE MEME KMER
                    }
                }
            }
        }
    }

    // VERIF SEUIL
    unordered_map<uint32_t, uint16_t>::iterator it = id_to_count.begin();
    while (it != id_to_count.end()) {
        if ((it->second) < threshold) {
            it = id_to_count.erase(it);
        } else
            it++;
    }
    return id_to_count;
}



vector<uint32_t> Index::verif_fp(const unordered_map<uint32_t, uint16_t>& reads_to_verify, const vector<uint32_t>& kmers, uint16_t threshold) const{

    string read_seq, curr_kmer;
    uint32_t kmer_int, revComp_int;
    uint64_t kmer_hash;
    vector<uint32_t> reads_to_return = {};
    for(auto read = reads_to_verify.begin(); read != reads_to_verify.end(); read++){
        uint cpt = 0;
        read_seq = get_read_sequence(read_line_pos, filename, read->first);
        bool is_first = true;
        string ligne, kmer;
        uint i = k;
        while(i<read_seq.length()){
            if(is_first){
                kmer = read_seq.substr(0, k);
                kmer_int = str2numstrand(kmer);
                is_first = false;
            }
            else{
                kmer_int = shiftnucl(kmer_int, read_seq[i], k);
                i++;
            }

            if (binary_search(kmers.begin(), kmers.end(), kmer_int)) {
                cpt++;
            }
            if(cpt == threshold){
                reads_to_return.push_back(read->first);
                break;
            }
        }
    }
    return reads_to_return;
}



void Index::memory_file(const string& output_file, uint16_t k, uint16_t bv_size, size_t uncomp_index_size) const{
  ofstream file(output_file, ios::out | ios_base::app);
  bm::bvector<>::statistics st;
  this->bf.calc_stat(&st);
  if(file){
    file << k << ";" << bv_size << endl;
    file << "K-mers:" << this->nb_kmers << endl;
    file << "Bf:" << st.memory_used << endl;
    file << "Counting bf:" << this->counting_bf->size() << endl;
    file << "Position vector:" << this->vect_pos.size() * 4 << endl;
    file << "Uncompressed index:" << uncomp_index_size << endl;
    file << "Compressed index:" << this->compressed_index.size() << endl;
  }
  else{
    cerr << "Error opening the file (memory_file)." << endl;
  }
}