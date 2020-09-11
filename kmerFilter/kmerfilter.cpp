#include <iostream>
#include <set>
#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <getopt.h>

void complement(std::string &seq)
{
    auto lambda = [](const char c) {
        switch (c)
        {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
	case 'N':
	    return 'N';
        default:
            throw std::domain_error("Invalide Nucleotide");
        }
    };
    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);
}

void MakeFilter(std::set<std::string> &filter, std::string &&filter_list_path, const bool stranded, const unsigned int k_length)
{
    std::cerr << "Making filter..." << std::endl;

    assert(filter.empty());

    std::ifstream filter_list(filter_list_path);
    std::string line;

    if (!filter_list.is_open())
    {
        std::cerr << "ERROR: Filter list file " << filter_list_path << " does not exist..." << std::endl;
        exit(EXIT_FAILURE);
    }
    while(std::getline(filter_list, line))
//    for (std::getline(filter_list, line); !filter_list.eof(); std::getline(filter_list, line))
    {
        if (line[0] == '>' || line.size() < k_length)
            continue;

        for (size_t seq_len = line.size(), start_pos = 0; start_pos < seq_len - k_length + 1; start_pos++)
        {
            std::string word = line.substr(start_pos, k_length);

            if (word.size() != k_length)
	    {
		std::cerr << word << "\t" << word.size() << "\t" << line << std::endl;
		exit(1);
	    }

            filter.insert(word);
            if (!stranded)
            {
                reverse(word.begin(), word.end());
                complement(word);
                filter.insert(word);
            }
        }
    }
    filter_list.close();

    std::cerr << "\t-> Filter size: " << filter.size() << std::endl;
}

void Filter(std::string &&kmer_count_path, const std::set<std::string> &filter_set, const bool keep_solid, const bool with_header)
{
    std::ifstream kmer_count_file(kmer_count_path);
    std::string line;

    if (!kmer_count_file.is_open())
    {
        std::cerr << "ERROR: Filter list file " << kmer_count_path << " does not exist..." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (with_header)
    {
        std::getline(kmer_count_file, line);
        std::cout << line << std::endl;
    }

    for (std::getline(kmer_count_file, line); !kmer_count_file.eof(); std::getline(kmer_count_file, line))
    {
        std::istringstream conv(line);
        std::string word;
        conv >> word;
        bool is_in_set = (filter_set.find(word) != filter_set.cend());

        if (is_in_set && keep_solid)
        {
            std::cout << line << std::endl;
        }
        else if (!is_in_set && !keep_solid)
        {
            std::cout << line << std::endl;
        }
    }
    kmer_count_file.close();
}

int main(int argc, char **argv)
{
    std::string filter_path;
    std::string raw_count_path;
    std::set<std::string> filter_set;
    bool stranded = true, keep_solid = true, with_header = true, has_filter_path = false;
    unsigned int k_length = 31;

    int opt;

    while ((opt = getopt(argc, argv, "hnk:lf:N")) != -1)
    {
        switch (opt)
        {
        case 'n':
            stranded = false;
            break;
        case 'k':
            k_length = atoi(optarg);
            break;
        case 'l':
            keep_solid = false;
            break;
        case 'f':
            filter_path = optarg;
            has_filter_path = true;
            break;
        case 'N':
            with_header = false;
            break;
        case 'h':
            std::cerr << "========= kmerfilter helper =========" << std::endl;
            std::cerr << "Usage: kmerfilter [-n] [-k k_length] [-l] -f filter_path kmer_count_table_path > output_sub_table_path" << std::endl << std::endl;
            std::cerr << "Parameter:    -n         if the k-mers are generated from unstranded RNA-seq data" << std::endl;
            std::cerr << "              -k INT     the length of k-mers [31]" << std::endl;
            std::cerr << "              -l         if to REMOVE the k-mers in filter list (to keep \"liquid\" in real filtering experiment)" << std::endl;
            std::cerr << "                         do not put this parameter if to SELECT the k-mers in filter list " 
                      << "(to keep \"solid\" in real filtering experiment)" << std::endl;
            std::cerr << "              -f STRING  the filter list path, could be either a fasta file or a list WITHOUT header" << std::endl;
            exit(EXIT_SUCCESS);
        }
    }
    if (!has_filter_path)
    {
        std::cerr << "ERROR: Filter list path is mandatory !" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (argc == optind)
    {
        std::cerr << "ERROR: The quantification table of k-mers is mandatory !" << std::endl;
        exit(EXIT_FAILURE);
    }
    raw_count_path = argv[optind++];

    std::cerr << "Length of k-mer: " << k_length << std::endl;
    std::cerr << "Stranded mode: " << (stranded ? "true" : "false") << std::endl;
    std::cerr << "Filter path: " << filter_path << std::endl;
    std::cerr << "K-mer count path: " << raw_count_path << std::endl;

    if (keep_solid)
    {
        std::cerr << "Filtering mode: keeping solid" << std::endl;
    }
    else
    {
        std::cerr << "Filtering mode: keeping liquid" << std::endl;
    }

    //exit(EXIT_FAILURE);

    MakeFilter(filter_set, std::move(filter_path), stranded, k_length);
    Filter(std::move(raw_count_path), filter_set, keep_solid, with_header);

    return EXIT_SUCCESS;
}
