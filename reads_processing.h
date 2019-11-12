//
// Created by lramirez on 11.11.19.
//

#ifndef BISQCK_READS_PROCESSING_H
#define BISQCK_READS_PROCESSING_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/modifier.h>

#include "utils.h"
typedef struct {
    int cpos;
    int methylated;
    int unmethylated;
} meth_caller_struct;


class kmer_match{
public:
    int k_beg;
    int tie_break_missmatch;
    std::vector<meth_caller_struct> cpgs;
    bool operator < (const kmer_match& other) const
    {
        return (tie_break_missmatch < other.tie_break_missmatch);
    }
};



inline bool operator== (const kmer_match &obj1, const kmer_match &obj2)
{
    return obj1.k_beg == obj2.k_beg &&
           obj1.cpgs.size() == obj2.cpgs.size() &&
           obj1.cpgs[0].cpos == obj2.cpgs[0].cpos;
}

int process_reads();
seqan::Dna5String recoverKmerFai(unsigned beg,  unsigned morechars = 0);
void process_read(seqan::Dna5String seq);
std::vector<kmer_match>process_kmer(seqan::Dna5String kmer);
int countdiff_chars(seqan::DnaString kmerref, seqan::DnaString kmerread);
#endif //BISQCK_READS_PROCESSING_H

