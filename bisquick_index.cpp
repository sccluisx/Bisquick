//
// Created by lramirez on 11.11.2019.
//

#include "bisquick_index.h"
#include "utils.h"

#include <seqan/seq_io.h>

#include <boost/range/irange.hpp>

using namespace seqan;

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
    int beg1 = cg_index->position_seq[0] - bisquickOptions->kmersize + 2;
    int end1 = cg_index->position_seq[0] + 1;
    cg_index->methylated.resize(cg_index->position_seq.size());
    cg_index->unmethylated.resize(cg_index->position_seq.size());
    cg_index->meth_ratio.resize(cg_index->position_seq.size());
    seqan::CharString ckmer1 = infixWithLength(seq, beg1, bisquickOptions->kmersize);
    for (auto &&cpos: cg_index->position_seq) {
        int beg = cpos - bisquickOptions->kmersize + 2;
        int end = cpos + 1;
        if (beg < 0)
            beg = 0;
        for (auto i : boost::irange(beg, end)) {
            seqan::CharString ckmer = infixWithLength(seq, i, bisquickOptions->kmersize);
            std::vector<int> cpg_list = find_cpgs(ckmer); // get all cpg positions with in the kmer
            map_value mapval = {
                    .beg = i, //position
                    .cpg_pos = cpg_list,
                    .sense = true,
            };

            if (compKMap.find(DnaString(compressKmer(ckmer))) == compKMap.end()) {
                std::vector<map_value> new_val;
                new_val.emplace_back(mapval);
                compKMap.insert({DnaString(compressKmer(ckmer)), new_val});

            } else {
                compKMap[DnaString(compressKmer(ckmer))].emplace_back(mapval);
            }
        }
    }
    long numclash, maxclash = 0;
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
    return 0;
}



