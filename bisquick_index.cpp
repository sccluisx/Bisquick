//
// Created by lramirez on 18.11.18.
//

#include "bisquick_index.h"
#include "utils.h"

#include <seqan/seq_io.h>

#include <boost/range/irange.hpp>

using namespace seqan;

//CG_Index cg_index;
std::shared_ptr<CG_Index> cg_index = std::make_shared<CG_Index>();
extern int c_id_counter;
CharString current_file;
std::map<Dna5String, map_value> compresskmap;
std::unordered_map<seqan::DnaString, std::vector<map_value>, DnaStringHasher> compKMap;
std::unordered_map<unsigned, unsigned> cgpos2index;

std::vector<int> find_cpgs(seqan::CharString haystack) {
    CharString needle = "CG";
    std::vector<int> cpg_pos;
    Finder<CharString> finder(haystack);
    Pattern<CharString, Horspool> pattern(needle);
    while (find(finder, pattern)) {
        uint cpos = beginPosition(finder);
        uint gpos = endPosition(finder);
        // save the position and the chromosome where "CG" is found
        //cg_index->seq_id.emplace_back(id);
        cpg_pos.emplace_back(cpos);
    }
    return cpg_pos;
}


int create_index(std::vector<std::string> fastafile) {
    int ksize = bisquickOptions->kmersize;
    current_file = fastafile[0];
    bisquickOptions->currentFile = current_file;
    CharString id;
    Dna5String seq;
    createFaiFile(toCString(current_file));
    SeqFileIn seqFileIn(toCString(current_file));
    readRecord(id, seq, seqFileIn);

    CharString haystack = seq;
    CharString haystack_reversec = seq;
    reverseComplement(haystack_reversec);
    CharString needle = "CG";
    Finder<CharString> finder(haystack);
    Pattern<CharString, Horspool> pattern(needle);
    int ncpgs = 0; // count cgs
    cg_index->seq_id = id;
    while (find(finder, pattern)) {
        uint cpos = beginPosition(finder);
        uint gpos = endPosition(finder);
        // save the position and the chromosome where "CG" is found
        cg_index->position_seq.emplace_back(cpos);
        cg_index->c_id.emplace_back(ncpgs);
        cgpos2index.insert({cpos, ncpgs});
        ncpgs++;
    }

    std::cout << "cpgs:" << ncpgs++ << std::endl;
    // std::cout << "first cpg:" << cg_index->position_seq[0] << std::endl;
    int beg1 = cg_index->position_seq[0] - bisquickOptions->kmersize + 2;
    int end1 = cg_index->position_seq[0] + 1;
    cg_index->methylated.resize(cg_index->position_seq.size());
    cg_index->unmethylated.resize(cg_index->position_seq.size());
    cg_index->meth_ratio.resize(cg_index->position_seq.size());
    seqan::CharString ckmer1 = infixWithLength(seq, beg1, bisquickOptions->kmersize);
    // std::cout << "first kmer:"<< ckmer1 << std::endl;
    // std::cout << "compress kmer:"<< compressKmer(ckmer1)<< std::endl;
    // std::cout << "compress Dnakmer:"<< DnaString(compressKmer(ckmer1))<< std::endl;

//    compKMap.insert({seqan::CharString("CTTT"),10});
//    for (auto itr = compKMap.begin(); itr != compKMap.end(); ++itr) {
//        std::cout << itr->first
//             << '\t' << itr->second << '\n';
//    }
    for (auto &&cpos: cg_index->position_seq) {
        int beg = cpos - bisquickOptions->kmersize + 2;
        int end = cpos + 1;
        if (beg < 0)
            beg = 0;
//        std::cout << "[" << beg << " - " << cpos << " - " << end << "]"
//                  << infixWithLength(seq, beg, bisquickOptions->kmersize) << std::endl;
        for (auto i : boost::irange(beg, end)) {
            seqan::CharString ckmer = infixWithLength(seq, i, bisquickOptions->kmersize);
            std::vector<int> cpg_list = find_cpgs(ckmer);
//            std::cout << "[" << i << "]";
//            std::cout << ckmer << ": ";
//            std::for_each(cpg_list.begin(), cpg_list.end(), [](const unsigned int i) { std::cout << i << " "; });
            map_value mapval = {
                    .beg = i,
                    .cpg_pos = cpg_list,
                    .sense = true,
            };

            if (compKMap.find(DnaString(compressKmer(ckmer))) == compKMap.end()) {
                //std::cout<<"Nuevo: "<<DnaString(compressKmer(ckmer))<<std::endl;
                std::vector<map_value> new_val;
                new_val.emplace_back(mapval);
                compKMap.insert({DnaString(compressKmer(ckmer)), new_val});

            } else {
                compKMap[DnaString(compressKmer(ckmer))].emplace_back(mapval);
                //std::cout<<"Viejo: "<<DnaString(compressKmer(ckmer))<<std::endl;
            }
            //compKMap.insert({DnaString(compressKmer(ckmer)), mapval});
//            std::cout << std::endl;
        }
//        std::cout << std::endl;
    }


/*        std::for_each(cg_index->position_seq.begin(), cg_index->position_seq.end(), [seq](int cpos){
        //get all kmers for every cg position calculate the borders
        unsigned beg = cpos - bisquickOptions->kmersize+2;
        unsigned end = cpos + 1;
        if(beg<0)
            beg = 0;
        for (auto i : boost::irange(beg, end)) {
            //infixWithLength(seq,i,bisquickOptions->kmersize);
            //auto cpg_list = find_cpgs(substr);
            std::cout<<"str: "<<" :"<<infixWithLength(seq,i,bisquickOptions->kmersize);
            //std::cout<<" "<<compressKmer(substr)<<" :";
            //std::for_each(cpg_list.begin(), cpg_list.end(), [] (const unsigned int i) {std::cout << i << " ";});
            std::cout<<std::endl;
            //compKMap.insert({compressKmer(substr),3});
        }
       // std::cout<<infixWithLength(seq,cg_index->position_seq[o]-bisquickOptions->kmersize+2,bisquickOptions->kmersize)<<std::endl;
    });*/
//    std::cout<<std::endl;
    int numclash, maxclash = 0;
    seqan::DnaString maxclashstr;
    for (auto itr = compKMap.begin(); itr != compKMap.end(); ++itr) {
        if (itr->second.size() > 1)
            numclash++;
        if (itr->second.size() > maxclash) {
            maxclashstr = itr->first;
            maxclash = itr->second.size();
        }

    }
    std::cout << std::endl << "***/*/*/*/*/*/*/*/*/**/*/*/*" << std::endl;
    std::cout << "numbe of non unique: " << numclash << std::endl;
    std::cout << "max clash_str: " << maxclashstr << " con: " << maxclash << std::endl;
    std::cout << std::endl << "***/*/*/*/*/*/*/*/*/**/*/*/*" << std::endl;
//        std::cout << itr->first << " value:"<<std::endl;
//        std::cout
//             << "tam:\t" << itr->second.size() << '\n';
////        std::for_each(itr->second.cpg_pos.begin(), itr->second.cpg_pos.end(), [] (const unsigned int i) {std::cout << i << " ";});
////        std::cout<< "\nsense:\t" << itr->second.sense << '\n';
//    }
    return 0;
}

/*
    // For each position where a cg exists we will generate all the kmers
    std::for_each(cg_index->position_seq.begin(), cg_index->position_seq.end(), [seq](int cpos) {
        //std::cout<<infixWithLength(seq,cg_index->position_seq[o]-bisquickOptions->kmersize+2,bisquickOptions->kmersize)<<std::endl;
        std::cout<<"("<<cpos<<")"<<std::endl;
        auto substr = infixWithLength(seq,cpos,bisquickOptions->kmersize);
        std::cout<<substr<<std::endl;
//        for(int i = 0; i< bisquickOptions->kmersize; i++)
//            std::cout<<"#"<<substr[i]<<std::endl;
        //std::cout<<infixWithLength(seq,o,bisquickOptions->kmersize)<<std::endl;
        //std::cout<<substr<<std::endl;
        //auto lambda =  []() { std::cout << "Code within a lambda expression" << std::endl; };
        //lambda();
        Dna5String t1= substr;
        //std::cout << "substr is: " << typeid(t1).name() << '\n';
        //std::unordered_map<Dna5String, int> umap;
        // umap[substr]=2;
        //std::unordered_map<Dna5String, int> umap;
//        for(int i=0;i<bisquickOptions->kmersize-1; i++){
//            for(int j = 0; j<i; j++)
//                std::cout<<" ";
//            std::cout<<infixWithLength(seq,o-bisquickOptions->kmersize+2+i,bisquickOptions->kmersize)<<std::endl;
//        }
    });*/



/*    for_each (myvector.begin(), myvector.end(), myfunction);

    std::cout<<"Testing: "<<cg_index->position_seq[3]<<std::endl;
    for(int i=0;i<ksize-1; i++){
        for(int j = 0; j<i; j++)
            std::cout<<" ";
        std::cout<<infixWithLength(seq,cg_index->position_seq[3]-ksize+2+i,ksize)<<std::endl;
    }*/
//int ls;
//std::cout<<cg_index->position_seq[3];
//printvec(cg_index->c_id,VAR_NAME(c_id));
//bisquickOptions->currfilesize = seqan::length(seq);
//    return 0;
//}

int createIntervalIndex() {

    for (int i = 0; i < cg_index->position_seq.size(); i++) {
        int ini = (int) cg_index->position_seq[i] - (bisquickOptions->kmersize - 1);
        std::cout << "<<" << ini
                  << "--" << cg_index->position_seq[i] << "--" <<
                  cg_index->position_seq[i] + (bisquickOptions->kmersize - 1) << std::endl;

    }


    return 0;
}

