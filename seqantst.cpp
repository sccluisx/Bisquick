#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <unordered_map>
#include <functional>
#include <map>
#include "utils.h"
#include <bitset>
#include <seqan/bed_io.h>
#include <boost/lexical_cast.hpp>


using namespace seqan;

struct Dna5StringHasher {
    size_t operator()(const seqan::Dna5String  & obj) const
    {
        std::string Dna5StdString(obj.data_begin,obj.data_end);
        std::hash<std::string> hash_fn;
        size_t hash = hash_fn(Dna5StdString);
        return hash;
    }
};

void write_bed(){
    CharString bedFileOutName = "/home/lramirez/thesis/bisqck/salida.bed";

    BedFileOut bedFileOut(toCString(bedFileOutName),'w');

    std::string str = boost::lexical_cast<std::string>(5.22);
    BedRecord<Bed4> record;

    // Fill and write out the first record.
    record.ref = "chr7";
    record.beginPos = 127471195;
    record.endPos = 127472363;
    record.data = CharString(str);
    writeRecord(bedFileOut, record);

    // Fill and write out the second record.

}

int main(int argc, char const ** argv) {
    OutputAndConsole file("output.txt");
    file << "test" ;
    return  0;
}

int main_s(int argc, char const ** argv) {

    auto current_file = argv[1];
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
    find(finder, pattern);
    std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

    if (argc != 5)
    {
        std::cerr << "USAGE: query_fai FILE.fa SEQ BEGIN END\n";
        return 0;
    }

    // Try to load index and create on the fly if necessary.
    FaiIndex faiIndex;
    if (!open(faiIndex, argv[1]))
    {
        if (!build(faiIndex, argv[1]))
        {
            std::cerr << "ERROR: Index could not be loaded or built.\n";
            return 0;
        }
        if (!save(faiIndex))    // Name is stored from when reading.
        {
            std::cerr << "ERROR: Index could not be written do disk.\n";
            return 0;
        }
    }

    // Translate sequence name to index.
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, argv[2]))
    {
        std::cerr << "ERROR: Index does not know about sequence " << argv[2] << "\n";
        return 0;
    }

    // Convert positions into integers.
    unsigned beginPos = 0, endPos = 0;
    if (!lexicalCast(beginPos, argv[3]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[3] << " into an unsigned.\n";
        return 0;
    }
    if (!lexicalCast(endPos, argv[4]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[4] << " into an unsigned.\n";
        return 0;
    }

    // Make sure begin and end pos are on the sequence and begin <= end.
    if (beginPos > sequenceLength(faiIndex, idx))
        beginPos = sequenceLength(faiIndex, idx);
    if (endPos > sequenceLength(faiIndex, idx))
        endPos = sequenceLength(faiIndex, idx);
    if (beginPos > endPos)
        endPos = beginPos;

    // Finally, get infix of sequence.
    Dna5String sequenceInfix;
    readRegion(sequenceInfix, faiIndex, idx, beginPos, endPos);
    std::cout << sequenceInfix << "\n";

    return 0;
}

int main_ccompare(int argc, char const **argv){
    Dna5String seq1("TTAA");
    Dna5 t_nucleotide('T');
    if(seq1[0]==t_nucleotide)
    {
        std::cout<<"Igual a C";
    }
    return 0;
}

int mainallreads(int argc, char const **argv){
    CharString seqFileName = "/home/lramirez/thesis/indata/chrm21/reads/simulated.fastq";

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(seqFileName)))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    try
    {
        readRecords(ids, seqs, seqFileIn);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    for (unsigned i = 0; i < length(ids); ++i)
        std::cout << ids[i] << '\t' << seqs[i] << '\n';

    return 0;
}

int mainxor(int argc, char const **argv){
    double start = seqan::cpuTime();
    DnaString seqFileName1 = "ACGT";
    DnaString seqFileName2 = "CAGT";
    DnaString seqFileName3 = "GCAT";
    DnaString seqFileName4 = "TCGA";
    //std::cout<<(seqFileName1.data_begin & seqFileName1.data_begin);
    std::string s = "ACGTTT"; // Some string.
    Dna5String prb  = "ACGTTT";
    Dna5String prb2 = "GCGTAN";
    DnaString as = "A";
    DnaString cs = "C";
    DnaString gs = "G";
    DnaString ts = "T";
    DnaString ns = "N";

    CharString tst = "ACGT";
    std::cout << (std::bitset<2>(as[0].value))<<std::endl;
    std::cout << (std::bitset<2>(cs[0].value))<<std::endl;
    std::cout << (std::bitset<2>(gs[0].value))<<std::endl;
    std::cout << (std::bitset<2>(ts[0].value))<<std::endl;
    std::cout << (std::bitset<2>(ns[0].value))<<std::endl;

    std::cout << (std::bitset<2>(as[0].value) ^ std::bitset<2>(ns[0].value) ).count()<< " ";
    std::cout << (std::bitset<2>(cs[0].value) ^ std::bitset<2>(ns[0].value) ).count() << " ";
    std::cout << (std::bitset<2>(gs[0].value) ^ std::bitset<2>(ns[0].value) ).count() << " ";
    std::cout << (std::bitset<2>(ts[0].value) ^ std::bitset<2>(ns[0].value) ).count() << " ";


    DnaString tst2 = tst;
    std::cout<<tst2<<std::endl;
    int sum = 0;
    for (int i = 0; i < length(prb); i++) {
        sum = sum + (std::bitset<3>(prb[i].value) ^ std::bitset<3>(prb2[i].value)).count();
        std::cout << (std::bitset<3>(prb[i].value) ^ std::bitset<3>(prb2[i].value) ) << " ";
    }

    std::cout<<sum<<std::endl;
    double end = seqan::cpuTime();
    std::cout.precision(17);
    std::cout<<"Execution time: "<<end-start;
    return 0;
}


int main_countcpg(int argc, char const **argv){
    double start = seqan::cpuTime();
    CharString seqFileName = "/home/lramirez/thesis/indata/chrm21/ref/chr21.fa";
    auto current_file = seqFileName;
    CharString id;
    Dna5String seq;
    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seq, seqFileIn);
    CharString haystack = seq;
    CharString haystack_reversec = seq;
    reverseComplement(haystack_reversec);
    CharString needle = "CG";
    Finder<CharString> finder(haystack);
    Pattern<CharString, Horspool> pattern(needle);
    int ncpgs = 0; // count cgs
    while (find(finder, pattern)) {
        uint cpos = beginPosition(finder);
        uint gpos = endPosition(finder);
        // save the position and the chromosome where "CG" is found
        //cg_index->seq_id.emplace_back(id);
        cg_index->position_seq.emplace_back(cpos);
        cg_index->c_id.emplace_back(ncpgs);
        ncpgs++;
    }
    std::cout<<"Total cpg:"<<ncpgs<<std::endl;
    double end = seqan::cpuTime();
    std::cout.precision(17);
    std::cout<<"Execution time: "<<end-start;

    return 0;
}

int main3(int argc, char const ** argv){
    CharString seq= "CGTC";
    CharString seq2= "CGTC";
    //unsigned long long tst1 = (unsigned long long)seq;

    //std::cout<<hamming_distance(seq,seq2);

//    std::unordered_map<Dna5String,int, Dna5StringHasher > u1;
//    u1.insert({\
//                {seq,2},
//                {seq,3},
//            {seq2,4}
//            });
//    Dna5String seq3= "acgttaaa";
//
//    std::cout<<u1[seq];
    return 0;
}


int hamming_distance(unsigned long long x, unsigned long long y)
{
    return __builtin_popcountll(x ^ y);
}

int main2(int argc, char const ** argv)
{
    if (argc != 5)
    {
        std::cerr << "USAGE: query_fai FILE.fa SEQ BEGIN END\n";
        return 0;
    }

    // Try to load index and create on the fly if necessary.
    FaiIndex faiIndex;
    if (!open(faiIndex, argv[1]))
    {
        if (!build(faiIndex, argv[1]))
        {
            std::cerr << "ERROR: Index could not be loaded or built.\n";
            return 0;
        }
        if (!save(faiIndex))    // Name is stored from when reading.
        {
            std::cerr << "ERROR: Index could not be written do disk.\n";
            return 0;
        }
    }

    // Translate sequence name to index.
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, argv[2]))
    {
        std::cerr << "ERROR: Index does not know about sequence " << argv[2] << "\n";
        return 0;
    }

    // Convert positions into integers.
    unsigned beginPos = 0, endPos = 0;
    if (!lexicalCast(beginPos, argv[3]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[3] << " into an unsigned.\n";
        return 0;
    }
    if (!lexicalCast(endPos, argv[4]))
    {
        std::cerr << "ERROR: Cannot cast " << argv[4] << " into an unsigned.\n";
        return 0;
    }

    // Make sure begin and end pos are on the sequence and begin <= end.
    if (beginPos > sequenceLength(faiIndex, idx))
        beginPos = sequenceLength(faiIndex, idx);
    if (endPos > sequenceLength(faiIndex, idx))
        endPos = sequenceLength(faiIndex, idx);
    if (beginPos > endPos)
        endPos = beginPos;

    // Finally, get infix of sequence.
    Dna5String sequenceInfix;
    readRegion(sequenceInfix, faiIndex, idx, beginPos, endPos);
    std::cout << sequenceInfix << "\n";


    return 0;
}

/*
////

#include <iostream>
#include <algorithm>
#include <string>

#include <seqan/find.h>
#include <seqan/modifier.h>
using namespace seqan;

int main()
{
    CharString needle2 = "GGGGaaaaaaaatttatataNt";
    CharString needle = needle2;

    reverseComplement(needle);
    reverseComplement(needle);

    for(auto it = needle.data_begin;  it!= needle.data_end; it++){
        std::cout << ' ' << *it;
    }
    std::cout<<std::endl;
    std::replace(needle.data_begin, needle.data_end, 'G', 'U');

    for(auto it = needle.data_begin;  it!= needle.data_end; it++){
        std::cout << ' ' << *it;
    }


    std::cout<<std::endl;
    std::cout<<needle;
    std::cout<<std::endl;
    std::cout<<needle2;
    return 0;
}*/
