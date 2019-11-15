/** @file main.cpp
 *  @title Main bisquick file
 *  @brief Main code of bisquick
 *  This function contains the main() function of bisquick
 * @version 0.1
 * @date 11/11/2019
 * @author Luis Enrique Ramirez Chavez
 * @bug No known bugs.
 */
#include <iostream>
#include <iomanip>
#include <memory>
#include <fstream>

#include <dirent.h>
#include <seqan/seq_io.h>
#include <boost/range/irange.hpp>
#include "bisquick_index.h"
#include "utils.h"
#include "reads_processing.h"

extern std::shared_ptr<CG_Index> cg_index;
extern seqan::CharString current_file;
extern std::shared_ptr<BisquickOptions> bisquickOptions;
int c_id_counter;
extern std::unordered_map<seqan::DnaString, std::vector<map_value>, DnaStringHasher> compKMap;


int main(int argc, char const **argv) {
    std::string output_sum_file(seqan::toCString(bisquickOptions->output));
    output_sum_file += "_summary_file_k_";
    output_sum_file += std::to_string(bisquickOptions->kmersize);
    output_sum_file += ".txt";
    OutputAndConsole t1(output_sum_file);

    c_id_counter = 0;
    // Parse the command line.
    auto parsecmdln = parseCommandLine(bisquickOptions, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (parsecmdln != seqan::ArgumentParser::PARSE_OK)
        return parsecmdln == seqan::ArgumentParser::PARSE_ERROR;
    std::cout << "|||Bisquick - fast methylation estimates|||" << std::endl;
    std::cout << "Bisquick options:" << std::endl;
    std::cout << "Kmer Size: " << bisquickOptions->kmersize << std::endl;
    std::cout << "Genome directory: " << seqan::toCString(bisquickOptions->genomePath) << std::endl;
    std::cout << "Reads directory: " << seqan::toCString(bisquickOptions->readsfile) << std::endl;


    auto fastafile = open_directory(seqan::toCString(bisquickOptions->genomePath));
    double start = seqan::cpuTime();
    create_index(fastafile);
    std::cout << "Index built" << std::endl;

    std::cout << "Processing reads" << std::endl;

    std::cout << "size of compressed kmer hash table: " << compKMap.size() << std::endl;

    //**************************************End of Building Index*******************************///
    double end_index = seqan::cpuTime();
    std::cout.precision(17);
    std::cout << "Execution building index time: " << end_index - start << std::endl;

    int k = process_reads();
    double end_reads = seqan::cpuTime();
    std::cout.precision(17);
    std::cout << "Execution processing reads time: " << end_reads - end_index << std::endl;
    double cpgs_covered = 0;
    double total_meth_ratio = 0;
    for (ulong i = 0; i < cg_index->methylated.size(); ++i) {
        if (cg_index->methylated[i] + cg_index->unmethylated[i] > 0) {
            cg_index->meth_ratio[i] =
                    (double) cg_index->methylated[i] / (double) (cg_index->methylated[i] + cg_index->unmethylated[i]);
            total_meth_ratio += cg_index->meth_ratio[i];
            cpgs_covered = cpgs_covered + 1;
        }
    }

    std::cout << "total meth ratio: " << total_meth_ratio / cpgs_covered << std::endl;
    std::cout << "total cpgs covered: " << cpgs_covered / cg_index->methylated.size() << std::endl;
    output_write_bed_file();
    double end = seqan::cpuTime();
    std::cout.precision(17);
    std::cout << "Execution time: " << end - start << std::endl;
    return 0;
}
