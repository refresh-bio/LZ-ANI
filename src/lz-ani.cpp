// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

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

enum class working_mode_t {none, all2all, pairs};

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

CParams params;

bool parse_params(int argc, char** argv);
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
	cerr << "   all2all              - all to all\n";
	cerr << "   pairs                - pairs\n";
	cerr << "Options:\n";
//	cerr << "   -bs                  - turn on all sequences buffering  (default: " << buffer_data << ")\n";		// hidden option
	cerr << "   --in-file-names <file_name> - file with list of input file names (all2all and one2all mode)\n";
	cerr << "   --one-file-name <file_name> - file name matched to the all files (one2all mode)\n";
	cerr << "   --one-file-name <file_name> - multi FASTA file name (all2all mode)\n";
	cerr << "   --in-pair-names <file_name> - file with list of pairs (space-separated) names\n";

	cerr << "   -out <file_name>     - output file name\n";
	cerr << "   -t <val>             - no of threads (default: " << params.no_threads << ")\n";
	cerr << "   -mml <val>           - min. match length (default: " << params.min_match_len << ")\n";
	cerr << "   -mdl <val>           - min. distant length (default: " << params.min_distant_match_len << ")\n";
	cerr << "   -cd <val>            - max. dist. between close matches (default: " << params.close_dist << ")\n";
	cerr << "   -mlrin <val>         - max. literal run len. in match (dafault: " << params.max_lit_run_in_match << ")\n";
	cerr << "   -cov <val>           - min. coverage threshold (default: " << params.min_coverage << ")\n";
	cerr << "   -reg <val>           - min. considered region length (default: " << params.min_region_len << ")\n";
	cerr <<	"   -aw <val>            - approx. window length (default: " << params.approx_window << ")\n";
	cerr << "   -am <val>            - max. no. of mismatches in approx. window (default: " << params.approx_mismatches << ")\n";
	cerr << "   -ar <val>            - min. length of run ending approx. extension (default: " << params.approax_run_len << ")\n";
	cerr << "   -filter <file_name> <min_val> - filtering in all-to-all mode\n";
//	cerr << "   --dense-output       - store results (in all to all mode) in a dense form\n";
	cerr << "   --verbose <int>      - verbosity level (default: " << params.verbosity_level << ")\n";
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

		if (par == "--in-file-names"s && i + 1 < argc)
		{
			input_file_names = load_input_names(argv[i+1]);
			if (input_file_names.empty())
				return false;

			i += 2;
		}
		else if (par == "--one-file-name"s && i + 1 < argc)
		{
			input_one_name = argv[i + 1];
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
		else if (par == "-bs")
		{
			buffer_data = true;
			i += 1;
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
			params.approax_run_len = atoi(argv[i + 1]);
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
/*		else if (par == "--dense-output"s)
		{
			params.output_dense_matrix = true;
		}*/
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
void prepare_worker_base(CSharedWorker* wb, int id)
{
	wb->clear_ref();

	if (!wb->load_reference(v_files_all2all[id].first, &(v_buffer_seqs[id])))
	{
		cerr << "Cannot read: " << v_files_all2all[id].first << endl;
		exit(0);
	}

	atomic<bool> long_ready = false;

	std::future<void> fut = std::async(std::launch::async, [&] {
		wb->prepare_kmers_ref_long();
		long_ready = true;
		long_ready.notify_one();
		wb->prepare_ht_long();
		});

	wb->prepare_kmers_ref_short();
	atomic_wait(&long_ready, false);
	wb->prepare_ht_short();

	fut.get();
}

// ****************************************************************************
void run_all2all_threads_mode()
{
	CLZMatcher lz_matcher(params);

	if (!input_file_names.empty())
		lz_matcher.init_data_storage(input_file_names);
	else
		lz_matcher.init_data_storage(input_one_name);

	lz_matcher.set_filter(filter_name, filter_thr);

	lz_matcher.run_all2all(output_file_name);

	return;
}

// ****************************************************************************
void run_all2all_sparse()
{
	CLZMatcher lz_matcher(params);

	if (!input_file_names.empty())
		lz_matcher.init_data_storage(input_file_names);
	else
		lz_matcher.init_data_storage(input_one_name);

	lz_matcher.reorder_input_files_sparse();
	lz_matcher.set_filter_map(filter_name, filter_thr);

	lz_matcher.run_all2all_sparse(output_file_name);

	return;
}

// ****************************************************************************
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
		return 0;

	if (params.no_threads == 0)
	{
		params.no_threads = thread::hardware_concurrency();
		if (!params.no_threads)
			params.no_threads = 1;
	}

	if (working_mode == working_mode_t::all2all)
	{
/*		if (!load_tasks_all2all())
			return 0;*/
//		run_all2all_threads_mode();
		run_all2all_sparse();
	}

	return 0;
}
