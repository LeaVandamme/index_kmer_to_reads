#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include "../BitMagic/src/bm.h"
#include "../BitMagic/src/bmundef.h"
#include "../BitMagic/src/bmsparsevec.h"
#include "../include/fastDelta.h"
#include "../headers/index.h"
#include "../headers/utils.h"
#include "../headers/bloomfilter.h"
#include "../TurboPFor-Integer-Compression/vp4.h"
#include "../SIMDCompressionAndIntersection/include/codecfactory.h"
#include "../SIMDCompressionAndIntersection/include/intersection.h"

using namespace std;
using namespace chrono;
using namespace SIMDCompressionLib;

int main(int argc, char *argv[])
{
    if (argc > 9 || argc < 2){
		cout << "[Fasta file (reads)] [Kmer length] [Bitvector Size] -k [Fasta file (1 kmer)] -s [Fasta file (1 sequence)] -f [Fasta file (several sequences)]" <<endl;
		exit(0);
    }

    else{

        string filename(argv[1]);
        uint16_t k = stoi(argv[2]);
        uint16_t bitvector_size = stoi(argv[3]);
        
        Index index = Index(filename, k);

        auto start_indexing = high_resolution_clock::now();
        
        vector<uint8_t> read_vector = index.index_fasta_rankselect_compressed(filename, k, bitvector_size);
        

        auto end_indexing = high_resolution_clock::now();
        auto indexing = duration_cast<seconds>(end_indexing - start_indexing);
        double pct_unique_kmers = (double)index.uniqueKmers()/index.get_nb_kmers()*100;

        cout << "\n";
        cout << "Indexing takes " << indexing.count() << " seconds." << endl;
        cout << "K-mer indexed : " << intToString(index.get_nb_kmers()) << " | " << "Unique k-mers : " << pct_unique_kmers << "%" << endl;
        uint64_t memory_used_index = getMemorySelfMaxUsed();
        cout << "Resource usage (indexing) : " << intToString(memory_used_index) << " Ko" << endl;
        cout << "Resource usage (for 1 k-mer) : " << koToBit((double)memory_used_index/index.get_nb_kmers()) << " Bits" << endl;

        int opt;

        while ((opt = getopt(argc, argv, "k:s:f:")) != EOF)
        {
            switch (opt)
            {

                case 'k':{

                        // QUERYING KMER
                    string kmer_to_search = seq_from_fasta(optarg);

                    auto start_querying_kmer = high_resolution_clock::now();
                    vector<uint32_t> res_query_kmer = index.query_kmer_rankselect(read_vector, kmer_to_search);
                    if(res_query_kmer.empty()){
                        cout << "The index does not contain this k-mer.";
                    }
                    auto end_querying_kmer = high_resolution_clock::now();
                    auto querying_kmer = duration_cast<seconds>(end_querying_kmer - start_querying_kmer);
                    cout << "Querying the k-mer takes " << querying_kmer.count() << " seconds." << endl;
                    break;
                }

                case 's':{
                    
                    // QUERYING SEQUENCE
                    string seq = seq_from_fasta(optarg);
                    auto start_querying_seq = high_resolution_clock::now();
                    index.query_sequence_rankselect(read_vector, seq);
                    auto end_querying_seq = high_resolution_clock::now();
                    auto querying_seq = duration_cast<seconds>(end_querying_seq - start_querying_seq);
                    cout << "Querying the sequence takes " << querying_seq.count() << " seconds." << endl;
                    break;
                }


                case 'f':{

                    // QUERYING FASTA
                    auto start_querying_fasta = high_resolution_clock::now();
                    index.query_fasta_rankselect(read_vector, optarg);
                    auto end_querying_fasta = high_resolution_clock::now();
                    auto querying_fasta = duration_cast<seconds>(end_querying_fasta - start_querying_fasta);
                    cout << "Querying the Fasta file takes " << querying_fasta.count() << " seconds." << endl;
                    break;
                }

                default:{
                    cout<<"[Fasta file (reads)] [Kmer length] -k [Fasta file (1 kmer)] -s [Fasta file (1 sequence)] -f [Fasta file (several sequences)]" <<endl;
                    exit(0);
                }
            }
        }
    }
    return 0;
}