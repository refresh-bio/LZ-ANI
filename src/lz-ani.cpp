// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

// all2all --in-file-names fl/ds50_0.txt -out aaa_new -t 32 -bs -cd 128 -reg 48 -mml 6 -mdl 11 -aw 16 -am 6 -ar 2  -filter kmer-db.db/ds50_0.a2a 3 --verbose 2

#include "defs.h"
#include "worker.h"
#include "s_worker.h"

#include "parallel-queues.h"

#include "app.h"
#include "lz_matcher.h"

#include <iostream>
#include <iomanip>
#include <istream>
#include <fstream>
#include <chrono>
#include <queue>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <thread>
#include <mutex>
#include <cstdio>
#include <future>
#include <atomic>
#include <fstream>
#include <algorithm>
#include <barrier>
#include <condition_variable>
#include <filesystem>
#include <random>
#include <iostream>


using namespace std::chrono;
using namespace refresh;
using namespace std;

struct pair_hash
{
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

enum class working_mode_t {none, all2all};

working_mode_t working_mode;
vector<string> input_file_names;
string input_one_name;
vector<pair<string,string>> input_pair_names;
string output_file_name;

CApplication app;

queue<pair<string, string>> q_files_pairs;
vector<pair<string, uint64_t>> v_files_all2all;
vector<pair<int, int>> v_files_all2all_order;
map<pair<string, string>, CFatResults> m_results;

mutex mtx_queue;
mutex mtx_res;

vector<pair<seq_t, int>> v_buffer_seqs;
vector<thread> v_threads;
vector<future<void>> v_fut;
//unordered_set<pair<int, int>, pair_hash> filter;
set<pair<int, int>> filter;
vector<int> filter_id_mapping;
string filter_name;
double filter_thr;
bool buffer_data = true;

CParams2 params;
//CParams2 params2;

bool parse_params(int argc, char** argv);
//bool parse_params2(int argc, char** argv);
vector<string> load_input_names(const string& fn);

void load_filter();
bool load_tasks_all2all();
void usage();

void run_all2all_threads_mode();
void run_all2all_sparse();

void split(const std::string& str, std::vector<std::string>& parts, char sep);

// ****************************************************************************
void usage()
{
	cerr << "lz-ani <mode> [options]\n";
	cerr << "Modes:\n";
	cerr << "   all2all                - all to all\n";
	cerr << "Options:\n";
	cerr << "   --in-fasta <file_name> - FASTA file (for multisample-fasta mode)\n";
	cerr << "   --in-txt <file_name>   - text file with FASTA file names\n";

	cerr << "   -out <file_name>       - output file name\n";
	cerr << "   -t <val>               - no of threads (default: " << params.no_threads << ")\n";
	cerr << "   -mml <val>             - min. match length (default: " << params.min_match_len << ")\n";
	cerr << "   -mdl <val>             - min. distant length (default: " << params.min_distant_match_len << ")\n";
	cerr << "   -cd <val>              - max. dist. between close matches (default: " << params.close_dist << ")\n";
	cerr << "   -mlrin <val>           - max. literal run len. in match (default: " << params.max_lit_run_in_match << ")\n";
	cerr << "   -cov <val>             - min. coverage threshold (default: " << params.min_coverage << ")\n";
	cerr << "   -reg <val>             - min. considered region length (default: " << params.min_region_len << ")\n";
	cerr <<	"   -aw <val>              - approx. window length (default: " << params.approx_window << ")\n";
	cerr << "   -am <val>              - max. no. of mismatches in approx. window (default: " << params.approx_mismatches << ")\n";
	cerr << "   -ar <val>              - min. length of run ending approx. extension (default: " << params.approx_run_len << ")\n";
	cerr << "   -filter <file_name> <min_val> - filtering file (kmer-db output) and threshold\n";
	cerr << "   --verbose <int>        - verbosity level (default: " << params.verbosity_level << ")\n";
	cerr << "   --multisample-fasta    - multi sample FASTA input (default: " << params.multisample_fasta << ")\n";

	cerr << "   --output-type <type>   - one of: 'single-file', 'split-files' (default: " << "" << ")\n";
	cerr << "   --store-total-ani      - store total ANI (default: " << "" << ")\n";
	cerr << "   --store-ani            - store ANI in both directions for (default: " << "" << ")\n";
	cerr << "   --store-cov            - store coverage in both directions for (default: " << "" << ")\n";
	cerr << "   --store-shorter-ani    - store ANI for shorter sequence (default: " << "" << ")\n";
	cerr << "   --store-shorter-cov    - store coverage for shorter sequence (default: " << "" << ")\n";
	cerr << "   --store-full-seq-ids   - store full sequence ids in main file (default: " << "" << ")\n";
}

// ****************************************************************************
vector<string> load_input_names(const string& fn)
{
	ifstream ifs(fn);
	vector<string> vec;

	if (!ifs.is_open())
	{
		cerr << "Cannot open file: " << fn << endl;
		return vec;
	}

	vec.assign(istream_iterator<string>(ifs), istream_iterator<string>());

	return vec;		
}

// ****************************************************************************
bool parse_params(int argc, char** argv)
{
	if (argc < 3)
	{
		usage();
		return false;
	}

	working_mode = working_mode_t::none;

	if (argv[1] == "all2all"s)
		working_mode = working_mode_t::all2all;
	else
	{
		cerr << "Unknown mode: " << argv[1] << endl;
		usage();
		return false;
	}

	for (int i = 2; i < argc;)
	{
		string par = string(argv[i]);

		if (par == "--in-txt"s && i + 1 < argc)
		{
			input_file_names = load_input_names(argv[i+1]);
			if (input_file_names.empty())
				return false;

			input_one_name.clear();

			i += 2;
		}
		else if (par == "--in-fasta"s && i + 1 < argc)
		{
			input_one_name = argv[i + 1];
			input_file_names.clear();
			i += 2;
		}
		else if (par == "-out"s)
		{
			if (i + 1 >= argc)
			{
				cerr << "Unknown out name\n";
				return false;
			}

			output_file_name = argv[i + 1];
			i += 2;
		}
		else if (par == "-t")
		{
			params.no_threads = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-mml")
		{
			params.min_match_len = atoi(argv[i + 1]);
			params.min_close_match_len = params.min_match_len;
			i += 2;
		}
		else if (par == "-mdl")
		{
			params.min_distant_match_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-cd")
		{
			params.close_dist = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-mlrim")
		{
			params.max_lit_run_in_match = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-cov")
		{
			params.min_coverage = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-reg")
		{
			params.min_region_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-aw")
		{
			params.approx_window = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-am")
		{
			params.approx_mismatches = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-ar")
		{
			params.approx_run_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-filter" && i + 2 < argc)
		{
			filter_name = argv[i + 1];
			filter_thr = atof(argv[i + 2]);
			i += 3;
		}
		else if (par == "--verbose"s && i + 1 < argc)
		{
			params.verbosity_level = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "--output-type"s && i + 1 < argc)
		{
			string par_type = argv[i + 1];
			if (par_type == "single-file"s)
				params.output_type = output_type_t::single_file;
			else if (par_type == "split-files"s)
				params.output_type = output_type_t::split_files;
			else
			{
				cerr << "Unknown output-type: " << par_type << endl;
				usage();
				exit(0);
			}
			i += 2;
		}
		else if (par == "--multisample-fasta"s)
		{
			params.multisample_fasta = true;
			++i;
		}
		else if (par == "--store-total-ani"s)
		{
			params.store_total_ani = true;
			++i;
		}
		else if (par == "--store-ani"s)
		{
			params.store_ani = true;
			++i;
		}
		else if (par == "--store-cov"s)
		{
			params.store_cov = true;
			++i;
		}
		else if (par == "--store-shorter-ani"s)
		{
			params.store_shorter_ani = true;
			++i;
		}
		else if (par == "--store-shorter-cov"s)
		{
			params.store_shorter_cov = true;
			++i;
		}
		else if (par == "--store-full-seq-ids"s)
		{
			params.store_full_seq_ids = true;
			++i;
		}
		else
		{
			cerr << "Unknown parameter: " << string(argv[i]) << endl;
			usage();
			exit(0);
		}
	}

	if (working_mode == working_mode_t::all2all && (input_file_names.empty() && input_one_name.empty()))
	{
		cerr << "Input file names not provided\n";
		return false;
	}

	return true;
}

// ****************************************************************************
void split(const std::string& str, std::vector<std::string>& parts, char sep)
{
	parts.clear();

	std::string s;

	for (auto c : str)
	{
		if (c == sep)
		{
			parts.emplace_back(s);
			s.clear();
		}
		else
			s.push_back(c);
	}

	if (!s.empty())
		parts.emplace_back(s);
}

// ****************************************************************************
void load_filter()
{
	ifstream ifs(filter_name);

	if (!ifs.is_open())
	{
		cerr << "Cannot open filter file: " << filter_name << endl;
		exit(0);
	}

	string line;
	vector<string> parts;
	vector<string> elem;

	getline(ifs, line);
	getline(ifs, line);

	for (int i = 0; !ifs.eof(); ++i)
	{
		getline(ifs, line);
		split(line, parts, ',');

		if (parts.size() <= 2)
			continue;

		for (const auto& p : parts)
		{
			split(p, elem, ':');
			if (elem.size() == 2)
			{
				int id = stoi(elem[0]) - 1;			// In kmer-db output indices are 1-based
				int val = stoi(elem[1]);

				if (val >= filter_thr)
					filter.insert(minmax(i, id));
			}
		}
	}

	cerr << "Filter size: " << filter.size() << endl;
}

// ****************************************************************************
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
		return 0;

	params.adjust_threads();
	params.multisample_fasta = true;

	CLZMatcher lzm(params);

	if (!lzm.load_sequences())
		return 0;
	
	if (!lzm.load_filter())
		return 0;

	if (!lzm.compare_sequences())
		return 0;

	lzm.reorder_sequences();

	lzm.run_all2all();

	return 0;
}
