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



uint32_t Index::get_nb_kmers(){
    return this->nb_kmers;
}



vector<uint32_t> vect_pos = {0};
bm::bvector<> bitvector;
bvector_rankselect rs_idx(new bm::bvector<>::rs_index_type());
vector<uint32_t> Index::index_fasta_rankselect(const string& read_file, uint16_t k, uint16_t bitvector_size = 32){

    bitvector.resize(pow(2,bitvector_size));
    uint32_t bvsize = bitvector.size();
    vector<uint32_t> index;
    uint32_t num_kmers = 1;
    uint32_t kmer_int, revComp_int, kmer_hash, revComp_hash;
    unordered_map<uint32_t, uint32_t> hmap;

    auto start_bitvector = high_resolution_clock::now();
    ifstream fichier(read_file, ios::in);
    if(fichier){

        // BLOOM FILTER
        
        while(!fichier.eof()){
            string ligne, kmer;
            getline(fichier,ligne);
            bool is_first = true;
            if (ligne[0] != '>' && !ligne.empty()){
                uint16_t i = k;
                while(i<ligne.length()){
                    if(is_first){
                        kmer = ligne.substr(0, k);
                        kmer_int = str2numstrand(kmer);
                        revComp_int = str2numstrand(revcomp(kmer));
                        kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                        is_first = false;
                    }
                    else{
                        kmer_int = shiftnucl(kmer_int, ligne[i], k);
                        revComp_int = shiftnucl_revcomp(revComp_int, ligne[i], k);
                        kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                        i++;
                    }
                    if(bitvector[kmer_hash] != 1){
                        bitvector.set(kmer_hash);
                        num_kmers++;
                    }
                }
            }
        }
        set_nb_kmers(num_kmers);
        bitvector.build_rs_index(rs_idx.get());
    }
    else{
        cerr << "Error opening the file." << endl;
    }
    cout << endl;
    fichier.close();
    auto end_bitvector = high_resolution_clock::now();
    auto querying_bitvector = duration_cast<std::chrono::seconds>(end_bitvector- start_bitvector);
    cout << "ok bitvector : " << querying_bitvector.count() << " seconds."  << endl;

    // VECTOR RANK/NB_APPARITION

    auto start_struct = high_resolution_clock::now();
    vector<uint32_t> kmer_app(num_kmers);
    ifstream fichier_struct(read_file, ios::in);
    if(fichier_struct){
        string ligne_struct, kmer_struct;
        while(!fichier_struct.eof()){
            bool is_first = true;
            getline(fichier_struct,ligne_struct);
            if (ligne_struct[0] != '>' && !ligne_struct.empty()){
                uint16_t i = k;
                while(i<ligne_struct.length()){
                    if(is_first){
                        kmer_struct = ligne_struct.substr(0, k);
                        kmer_int = str2numstrand(kmer_struct);
                        revComp_int = str2numstrand(revcomp(kmer_struct));
                        kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                        is_first = false;
                    }
                    else{
                        kmer_int = shiftnucl(kmer_int, ligne_struct[i], k);
                        revComp_int = shiftnucl_revcomp(revComp_int, ligne_struct[i], k);
                        kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                        i++;
                    }
                    auto rank_hash = bitvector.rank(kmer_hash, *rs_idx);
                    kmer_app[rank_hash-1] ++;
                    } 
                }
            }
        }

    fichier_struct.close();
    auto end_struct = high_resolution_clock::now();
    auto querying_struct = duration_cast<std::chrono::seconds>(end_struct - start_struct);
    cout << "ok struct : " << querying_struct.count() << " seconds."  << endl;

    // VECTOR WITH START POSITIONS FOR READS

    uint32_t total_count = 0;
    for(vector<uint32_t>::iterator i = kmer_app.begin() ; i != kmer_app.end(); i++){
        total_count += *i;
        vect_pos.push_back(total_count);
    }
    cout << "ok vectpos" << endl;
    
    index.resize(total_count);
    uint cpt = 0;
        
    auto start_index = high_resolution_clock::now();
    uint32_t num_read_index = 0;
    ifstream fichier_index(read_file, ios::in);
    if(fichier_index){
        string ligne_index, kmer_index;
        while(!fichier_index.eof()){
            bool is_first = true;
            getline(fichier_index,ligne_index);
            if (ligne_index[0] != '>' && !ligne_index.empty()){
                uint16_t i = k;
                num_read_index ++;

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
                    auto rank_hash = bitvector.rank(kmer_hash, *rs_idx);
                    uint32_t start_pos = vect_pos[rank_hash - 1];
                    uint32_t length_max = vect_pos[rank_hash] - vect_pos[rank_hash-1];

                    bool present = false;
                    uint tmp = start_pos;
                    while(tmp<start_pos+length_max && !present){
                        if(index[tmp] == false){
                            index[tmp] = num_read_index;
                            present = true; 
                        }
                        tmp++;
                    }
                } 
            }
        }
    }
    fichier_index.close();

    auto end_index = high_resolution_clock::now();
    auto querying_index = duration_cast<std::chrono::seconds>(end_index - start_index);
    cout << "ok index : " << querying_index.count() << " seconds."  << endl;
    return index;
}


vector<uint8_t> Index::index_fasta_rankselect_compressed(const string& read_file, uint16_t k, uint16_t bitvector_size){
    bitvector.resize(pow(2,bitvector_size));
    uint32_t bvsize = bitvector.size();
    vector<uint32_t>* index=new vector<uint32_t>;
    uint32_t num_kmers = 1;
    uint32_t kmer_int, revComp_int, kmer_hash, revComp_hash;
    auto start_bitvector = high_resolution_clock::now();
    using value_type = bm::sparse_vector<int, bm::bvector<> >::value_type;
    sparse_vector_u32* occ_of_kmer = new sparse_vector_u32();
    occ_of_kmer -> resize(pow(2,bitvector_size));
    ifstream fichier(read_file, ios::in);
    if(fichier){

        // BLOOM FILTER
        
        while(!fichier.eof()){
            string ligne, kmer;
            getline(fichier,ligne);
            bool is_first = true;
            if (ligne[0] != '>' && !ligne.empty()){
                uint16_t i = k;
                while(i<ligne.length()){
                    if(is_first){
                        kmer = ligne.substr(0, k);
                        kmer_int = str2numstrand(kmer);
                        revComp_int = str2numstrand(revcomp(kmer));
                        kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                        is_first = false;
                        i++;
                    }
                    else{
                        kmer_int = shiftnucl(kmer_int, ligne[i], k);
                        revComp_int = shiftnucl_revcomp(revComp_int, ligne[i], k);
                        kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                        i++;
                    }
                    if(bitvector[kmer_hash] != 1){
                        bitvector.set(kmer_hash);
                        num_kmers++;
                    }
                    occ_of_kmer -> inc(kmer_hash);
                }
            }
        }
        set_nb_kmers(num_kmers);
        bitvector.build_rs_index(rs_idx.get());
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

    for(auto it = occ_of_kmer->begin(); it != occ_of_kmer -> end(); it++){
		if(*it != 0){
            total_count += *it;
            vect_pos.push_back(total_count);
        }
	}//ON PEUX SUREMENT OPTIMISER CA AVEC bm::sparse_vector_find_first_mismatch

    cout << "ok vectpos" << endl;
    occ_of_kmer -> clear();
    delete(occ_of_kmer);
    index->resize(total_count);

    // INDEX
        
    auto start_index = high_resolution_clock::now();
    uint32_t num_read_index = 0;
    ifstream fichier_index(read_file, ios::in);
    if(fichier_index){
        string ligne_index, kmer_index;
            while(!fichier_index.eof()){
                bool is_first = true;
                getline(fichier_index,ligne_index);
                if (ligne_index[0] != '>' && !ligne_index.empty()){
                    uint16_t i = k;
                    num_read_index ++;

                    while(i<ligne_index.length()){
                        if(is_first){
                            kmer_index = ligne_index.substr(0, k);
                            kmer_int = str2numstrand(kmer_index);
                            revComp_int = str2numstrand(revcomp(kmer_index));
                            kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                            is_first = false;
                            i++;
                        }
                        else{
                            kmer_int = shiftnucl(kmer_int, ligne_index[i], k);
                            revComp_int = shiftnucl_revcomp(revComp_int, ligne_index[i], k);
                            kmer_hash = find_canonical(kmer_int, revComp_int, bvsize);
                            i++;
                        }
                        auto rank_hash = bitvector.rank(kmer_hash, *rs_idx);
                        uint32_t start_pos = vect_pos[rank_hash - 1];
                        uint32_t nb_apparition = vect_pos[rank_hash] - vect_pos[rank_hash-1];

                        (*index)[start_pos]++;
                        uint32_t indice_to_insert = start_pos + (*index)[start_pos];
                        if(indice_to_insert < (start_pos + nb_apparition)){
                            (*index)[indice_to_insert] = num_read_index;
                        }
                        else{
                            (*index)[start_pos] = num_read_index;
                        }
                    } 
                }
        }
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
    vector<uint32_t>::const_iterator begin = index->begin();
    for(uint32_t i = 1; i<vect_pos.size(); i++){
        first = begin + vect_pos[i-1];
        last = begin + vect_pos[i];
        
        vector<uint32_t> sub_index(first, last);
        
        // compression du bloc
        compressed_vector = vector<unsigned char>(sub_index.size()*16 + 100); 
        scomp=p4ndenc32(sub_index.data(), sub_index.size() , compressed_vector.data());
        compressed_vector.resize(scomp);
        // ajout de chaque élément compressé dans l'index
        for(uint32_t j=0; j<scomp; j++){ 
            compressed_id.push_back(compressed_vector[j]);
            //(*index)[tmp_indice_vect_pos+j] = compressed_vector[j];
        }
        // modif du vecteur avec les positions d'insertion des numéros de read pour chaque bloc
        tmp_indice_vect_pos += scomp;
        vect_pos[i-1] = tmp_indice_vect_pos;
        
    }
    delete index;
    cout << "ok index compressé" << endl;

    auto end_index = high_resolution_clock::now();
    auto querying_index = duration_cast<std::chrono::seconds>(end_index - start_index);
    cout << "ok index : " << querying_index.count() << " seconds."  << endl;
    bm::bvector<>::statistics st;
    bitvector.calc_stat(&st);
    cout << "Taille du Bloom filter : " << st.memory_used << endl;
    cout << "Taille du vecteur positions : " << vect_pos.size() * 4 <<endl;
    cout << "Taille du vecteur identifiants : " << compressed_id.size() << endl;
    return compressed_id;
}





vector<uint32_t> Index::query_kmer_rankselect(const vector<uint32_t>* index, const string& kmer){
    uint32_t kmer_int = str2numstrand(kmer);
    uint32_t revComp_int = str2numstrand(revcomp(kmer));
    uint64_t kmer_hash = find_canonical(kmer_int, revComp_int, bitvector.size());
    uint length = 0;
    uint start = 0;
    if (bitvector[kmer_hash] == 1){ // if the kmer exists
        auto rank_hash = bitvector.rank(kmer_hash, *rs_idx);

        vector<uint32_t>::const_iterator first, last;
        vector<uint32_t>::const_iterator begin = index -> begin();
        first = begin + vect_pos[rank_hash-1];
        last = begin + vect_pos[rank_hash];
        vector<unsigned char> to_decompress(first, last);

        vector<uint32_t> uncompressed_vector;
        uint scomp;

        uncompressed_vector = vector<uint32_t>(to_decompress.size()*16 + 100);
        scomp=p4nddec32(to_decompress.data(), to_decompress.size() , uncompressed_vector.data());
        uncompressed_vector.resize(scomp);
        return uncompressed_vector;
    }
    else{
        cout << "K-mer not found." << '\n';
        vector<uint32_t> v;
        return v;
    }
}



vect_occ_read Index::query_sequence_rankselect(const vector<uint32_t>* index, const string& sequence){

    vect_occ_read res;
    unordered_map<uint32_t, uint32_t> tmp_reads;
    bool is_first = true;
    for (long unsigned int i=0; i<=sequence.length()-this->k; i++){
        string kmer = sequence.substr(i, this->k);
        vector<uint32_t> reads_of_kmer = query_kmer_rankselect(index, kmer);
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




vect_occ_read Index::query_fasta_rankselect(const vector<uint32_t>* index, const string& filename){
    
    vect_occ_read res;
    unordered_map<uint32_t, uint32_t> tmp_reads;

    ifstream fichier(filename, ios::in);

    if(fichier){   
        string ligne;
        while( !fichier.eof()){
            getline(fichier,ligne);
            if (ligne[0] != '>' && !ligne.empty()){
                vect_occ_read occ_read = query_sequence_rankselect(index, ligne);
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
        cerr << "Error opening the file." << endl;
        return res;
    }
}


void Index::write_index_rankselect(const vector<uint32_t>& index, const string& output_file){

    ofstream outp(output_file);
    cout << "Writing the file..." << endl;
    if (outp){
        uint rank, pos_begin, length;
        uint kmer_hash = bitvector.get_first();
        bool first = true;

        while(kmer_hash != 0 || first){
            rank = bitvector.rank(kmer_hash, *rs_idx);
            pos_begin = vect_pos[rank-1];
            length = vect_pos[rank] - vect_pos[rank-1];
            outp << kmer_hash << " : ";
            for(uint i=pos_begin; i<pos_begin+length; i++){
                outp << index[i] << ",";
            }
            outp << "\n";
            first = false;
            kmer_hash = bitvector.get_next(kmer_hash);
        }
    }
    outp.close();
}



uint32_t Index::uniqueKmers(){
    uint32_t nb_unique = 0;
    for (uint i=0; i<vect_pos.size()-1; i++) {
        if((vect_pos[i+1] - vect_pos[i]) == 1){
            nb_unique ++;
        }
    }
    return nb_unique;
}