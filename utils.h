/**
 * @file utils.h
 * @version 0.1
 * @date 20/11/2018
 * @author Luis Enrique Ramirez Chavez
 * @title Utility functions and structures
 */

#ifndef BISQCK_UTILS_H
#define BISQCK_UTILS_H

#include <iostream>
#include <algorithm>
#include <memory>
#include <seqan/arg_parse.h>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "bisquick_index.h"

/**
 * @brief  Structure to contain the options from command line
 */
struct BisquickOptions
{
    int kmersize;
    int missmatch_tol;
    seqan::CharString genomePath;
    seqan::CharString currentFile;
    seqan::CharString currentFileFai;
    seqan::CharString readsfile;
    seqan::CharString output;

    int currfilesize;
    inline BisquickOptions() :
            kmersize(30),
            missmatch_tol(0),
            currentFile(""),
            currentFileFai(""),
            //genomePath("/home/lramirez/thesis/indata/toy"),
            genomePath("/home/lramirez/thesis/bisqck/Release/indata/chrm21/ref/"),
            readsfile("/home/lramirez/thesis/bisqck/Release/indata/chrm21/reads/simulated.fastq"),
            output("./")

    {}
} ;

extern int c_id_counter;



extern std::shared_ptr<BisquickOptions> bisquickOptions;
extern std::shared_ptr<CG_Index> cg_index;


seqan::ArgumentParser::ParseResult
parseCommandLine(std::shared_ptr<BisquickOptions> options, int argc, char const ** argv);

/**
 * @brief The function hasEnding check if the string ends with a
 *          particular substring
 * @param fullString full string that we want check
 * @param ending the ending part
 * @return True if the string ends with the substring
 */

bool hasEnding (std::string const &fullString, std::string const &ending);


int createFaiFile(char const * fastafile);

int kmerRecover(char const * fastafile, seqan::CharString seqid, int ini, int fin);

std::vector<std::string> open_directory(char * curdir);

/**
 * @brief This functions help to convert all the kmers
 * @param kmer
 * @return The reduced kmer convertings all C's in T's
 */
seqan::CharString compressKmer(seqan::CharString kmer);
// --------------------------------------------------------------------------
// Function  print_mainstructs()
// --------------------------------------------------------------------------
// for debug
#define VAR_NAME(n) #n
template <class T>
void printvec(std::vector<T> const &myvec, const std::string &vecname ) {
    typename std::vector<T>::const_iterator it;
    std::cout <<vecname<< " vector contains:";
    if (typeid(unsigned char) == typeid(T))
        for( it = myvec.begin(); it < myvec.end(); it++) {
            std::cout << " " <<(seqan::CharString) *it;
        } else {
        for( it = myvec.begin(); it < myvec.end(); it++) {
            std::cout << " " << *it;
        }
    }
    std::cout << std::endl;
}

struct CharStringHasher {
    size_t operator()(const seqan::CharString  & obj) const
    {
        std::string CharStdString(obj.data_begin,obj.data_end);
        std::hash<std::string> hash_fn;
        size_t hash = hash_fn(CharStdString);
        return hash;
    }
};

struct DnaStringHasher {
    size_t operator()(const seqan::DnaString  & obj) const
    {
        std::string DnaStdString(obj.data_begin,obj.data_end);
        std::hash<std::string> hash_fn;
        size_t hash = hash_fn(DnaStdString);
        return hash;
    }
};




void output_write_bed_file();


struct OutputAndConsole : std::ofstream
{
    OutputAndConsole(const std::string& fileName)
            : std::ofstream(fileName)
            , fileName(fileName)
    {};

    const std::string fileName;
};

template <typename T>
OutputAndConsole& operator<<(OutputAndConsole& strm, const T& var)
{
    std::cout << var;
    static_cast<std::ofstream&>(strm) << var;
    return strm;
};



#endif //BISQCK_UTILS_H
