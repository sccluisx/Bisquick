//
// Created by lramirez on 06.11.19.
//

#include <boost/range/irange.hpp>
#include "reads_processing.h"

#include "utils.h"

extern std::shared_ptr<CG_Index> cg_index;
extern std::shared_ptr<BisquickOptions> bisquickOptions;
extern std::unordered_map<seqan::DnaString, std::vector<map_value>, DnaStringHasher> compKMap;
extern std::unordered_map<unsigned, unsigned> cgpos2index;

int process_reads() {
    // std::cout<<"Reads file:"<<bisquickOptions->readsfile<<std::endl;
    seqan::CharString seqFileName = bisquickOptions->readsfile;
    seqan::SeqFileIn seqFileIn;
    if (!seqan::open(seqFileIn, seqan::toCString(seqFileName))) {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    try { seqan::readRecords(ids, seqs, seqFileIn); }
    catch (seqan::Exception const &e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    for (unsigned i = 0; i < seqan::length(ids); ++i) {
    //for (unsigned i = 0; i < 10; ++i) {
        if (i % (seqan::length(ids)/100) == 0) {
            std::cout << "read [" << i << "/" << seqan::length(ids) << "]:\t" << ids[i] << '\n';
        }
        process_read(seqs[i]);
    }
    return 0;
}

void process_read(seqan::Dna5String seq) {
    int seqlen = seqan::length(seq);
    for (auto i : boost::irange(0, seqlen - bisquickOptions->kmersize + 1)) {
        std::vector<kmer_match> possible_calls;
        seqan::Dna5String kmer = infixWithLength(seq, i, bisquickOptions->kmersize);
        possible_calls = process_kmer(kmer);
        possible_calls.erase(unique(possible_calls.begin(), possible_calls.end()), possible_calls.end());
        if (possible_calls.size() > 1) {
            for (auto &possible_call : possible_calls) {
                int c_read_pos = possible_call.k_beg - i;
                seqan::Dna5String seqref = recoverKmerFai(c_read_pos, seqlen - bisquickOptions->kmersize);
                possible_call.tie_break_missmatch = countdiff_chars(seqan::DnaString(seqref), seqan::DnaString(seq));
            }
        }
        auto best_call = std::min_element(possible_calls.begin(),possible_calls.end());
        if (!possible_calls.empty()) {
            for (auto cpg: best_call->cpgs) {
                int caddress = cgpos2index[cpg.cpos + best_call->k_beg];
                cg_index->methylated[caddress] = cg_index->methylated[caddress] + cpg.methylated;
                cg_index->unmethylated[caddress] = cg_index->unmethylated[caddress] + cpg.unmethylated;
            }
        }
    }

}

std::vector<kmer_match> process_kmer(seqan::Dna5String kmer) {
    seqan::DnaString base4kmerStr(kmer);
    seqan::DnaString compresskmer = seqan::DnaString(compressKmer(base4kmerStr));
    seqan::Dna5 t_nucleotide('T');
    std::vector<kmer_match> possible_calls;
    if (compKMap.find(compresskmer) == compKMap.end()) {
    } else {
        // When the reduced kmer is in the map
        // check all possible matches
        for (int i = 0; i < compKMap[compresskmer].size(); ++i) {
            kmer_match possible_match;
            possible_match.k_beg = compKMap[compresskmer][i].beg;
            int miss_matchs = 0;
            int max_miss_matchs = compKMap[compresskmer][i].cpg_pos.size() + bisquickOptions->missmatch_tol;
            seqan::Dna5String refkmer = recoverKmerFai(compKMap[compresskmer][i].beg);
            int diffs = countdiff_chars(seqan::DnaString(refkmer), seqan::DnaString(kmer));
            if (diffs <= max_miss_matchs) {
                std::vector<meth_caller_struct> curr_cgs(compKMap[compresskmer][i].cpg_pos.size(), {0, 0, 0});
                // check  if there is a conversion in al cpg positions of the kmer
                for (int j = 0; j < compKMap[compresskmer][i].cpg_pos.size(); ++j) {
                    int cpos = compKMap[compresskmer][i].cpg_pos[j];
                    curr_cgs[j].cpos = cpos;
                    if (kmer[cpos] == t_nucleotide) {
                        curr_cgs[j].unmethylated++;
                        miss_matchs++;
                    } else {
                        curr_cgs[j].methylated++;
                    }
                }
                if (diffs == miss_matchs) {
                    possible_match.cpgs.resize(curr_cgs.size());
                    possible_match.cpgs = curr_cgs;
                    possible_calls.emplace_back(possible_match);
                }
            }
        }
    }
    return possible_calls;
}

int countdiff_chars(seqan::DnaString kmerref, seqan::DnaString kmerread) {
    int sum = 0;
    for (int i = 0; i < length(kmerref); i++) {
        sum = sum + (std::bitset<2>(kmerref[i].value) ^ std::bitset<2>(kmerread[i].value)).count();
    }
    return sum;
}

seqan::Dna5String recoverKmerFai(unsigned beg, unsigned morechars) {
    auto current_file = seqan::toCString(bisquickOptions->currentFile);
    auto seq_id = cg_index->seq_id;
    unsigned end = beg + bisquickOptions->kmersize + morechars;
    // Try to load index and create on the fly if necessary.
    seqan::FaiIndex faiIndex;
    if (!seqan::open(faiIndex, current_file)) {
        if (!seqan::build(faiIndex, current_file)) {
            std::cerr << "ERROR: Index could not be loaded or built.\n";
            return 0;
        }
        if (!seqan::save(faiIndex))    // Name is stored from when reading.
        {
            std::cerr << "ERROR: Index could not be written do disk.\n";
            return 0;
        }
    }
    // Translate sequence name to index.
    unsigned idx = 0;
    if (!seqan::getIdByName(idx, faiIndex, seq_id)) {
        std::cerr << "ERROR: Index does not know about sequence " << seq_id << "\n";
        return 0;
    }
    // Convert positions into integers.
    unsigned beginPos = beg;
    unsigned endPos = end;
    // Make sure begin and end pos are on the sequence and begin <= end.
    if (beginPos > seqan::sequenceLength(faiIndex, idx))
        beginPos = seqan::sequenceLength(faiIndex, idx);
    if (endPos > seqan::sequenceLength(faiIndex, idx))
        endPos = seqan::sequenceLength(faiIndex, idx);
    if (beginPos > endPos)
        endPos = beginPos;
    // Finally, get infix of sequence.
    seqan::Dna5String sequenceInfix;
    seqan::readRegion(sequenceInfix, faiIndex, idx, beginPos, endPos);
    return sequenceInfix;
}
