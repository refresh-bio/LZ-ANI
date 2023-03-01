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

enum class working_mode_t {none, all2all, one2all, pairs};

working_mode_t working_mode;
vector<string> input_file_names;
string input_one_name;
vector<pair<string,string>> input_pair_names;
string output_file_name;

CApplication app;

queue<pair<string, string>> q_files_pairs;
vector<pair<string, uint64_t>> v_files_all2all;
vector<pair<int, int>> v_files_all2all_order;
//vector<string> v_files_one2all;
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
int filter_thr;
bool buffer_data = true;

CParams params;

bool parse_params(int argc, char** argv);
vector<string> load_input_names(const string& fn);



void load_filter();
void load_tasks_pairs();
bool load_tasks_all2all();
//void load_tasks_one2all();
void usage();
void run_pairs_mode();
void run_one2all_mode();
void split(const std::string& str, std::vector<std::string>& parts, char sep);

// ****************************************************************************
void usage()
{
	cerr << "lz-ani <mode> [options]\n";
	cerr << "Modes:\n";
	cerr << "   all2all              - all to all\n";
	cerr << "   one2all              - one to all\n";
	cerr << "   pairs                - pairs\n";
	cerr << "Options:\n";
//	cerr << "   -bs                  - turn on all sequences buffering  (default: " << buffer_data << ")\n";		// hidden option
	cerr << "   --in-file-names <file_name> - file with list of input file names (all2all and one2all mode)\n";
	cerr << "   --one-file-name <file_name> - file name matched to the all files (one2all mode)\n";
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
//	cerr << "   -rng-f <val>         - start of range of sequences to compare with all other (default: " << RANGE_FROM << ")\n";	// hidden option
//	cerr << "   -rng-t <val>         - end of range of sequences to compare with all other (default: " << RANGE_TO << ")\n";		// hidden option
	cerr << "   -flt <file_name> <min_val> - filtering in all-to-all mode\n";
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
	else if(argv[1] == "one2all"s)
		working_mode = working_mode_t::one2all;
	else if(argv[1] == "pairs"s)
		working_mode = working_mode_t::pairs;
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
		else if (par == "--in-pair-names"s && i + 1 < argc)
		{
			auto vec = load_input_names(argv[i + 1]);
			if (input_file_names.empty())
				return false;

			input_pair_names.clear();
			vector<string> tmp;
			for (auto& x : vec)
			{
				split(x, tmp, ' ');
				if (tmp.size() < 2)
				{
					cerr << "Wrong pair of file names: " << x << endl;
					return false;
				}

				input_pair_names.emplace_back(tmp[0], tmp[1]);
			}

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
/*		else if (par == "-rng-f")
		{
			RANGE_FROM = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-rng-t")
		{
			RANGE_TO = atoi(argv[i + 1]);
			i += 2;
		}*/
		else if (par == "-filter" && i + 2 < argc)
		{
			filter_name = argv[i + 1];
			filter_thr = atoi(argv[i + 2]);
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

/*	if (!filter_name.empty())
		load_filter();*/


	if (working_mode == working_mode_t::all2all && input_file_names.empty())
	{
		cerr << "Input file names not provided\n";
		return false;
	}
	else if (working_mode == working_mode_t::one2all && input_file_names.empty())
	{
		cerr << "Input file names not provided\n";
		return false;
	}
	else if (working_mode == working_mode_t::one2all && input_one_name.empty())
	{
		cerr << "One input file name not provided\n";
		return false;
	}
	else if (working_mode == working_mode_t::pairs && input_pair_names.empty())
	{
		cerr << "Input pair names not provided\n";
		return false;
	}

	return true;
}

// ****************************************************************************
/*void load_tasks_pairs()
{
	FILE *f = fopen(input_name.c_str(), "rb");
	if (!f)
	{
		cerr << "Error: Cannot load " + input_name + "\n";
		exit(0);
	}
	
	while (!feof(f))
	{
		char s[256], t[256];
		if (fscanf(f, "%s%s", s, t) == EOF)
			break;
		if (feof(f))
			break;

		q_files_pairs.push(make_pair(string(s), string(t)));
	}

	fclose(f);
}*/

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
bool load_tasks_all2all()
{
/*	FILE* f = fopen(input_name.c_str(), "rb");
	if (!f)
	{
		cerr << "Error: Cannot load " + input_name + "\n";
		exit(0);
	}*/

	unordered_map<string, int> map_file_list;
	int id = 0;
	std::error_code ec;

	for(const auto &fn : input_file_names)
	{
		auto fs = std::filesystem::file_size(std::filesystem::path(fn), ec);

		if (ec)
		{
			cerr << fn << " : " << ec.message() << endl;
			return false;
		}

		v_files_all2all.emplace_back(string(fn), fs);
		map_file_list[fn] = id++;
	}
		
	stable_sort(v_files_all2all.begin(), v_files_all2all.end(), [](const auto& x, const auto& y) {
		return x.second > y.second; });

	filter_id_mapping.resize(v_files_all2all.size());

	for (int i = 0; i < v_files_all2all.size(); ++i)
		filter_id_mapping[i] = map_file_list[v_files_all2all[i].first];

	if (buffer_data)
		v_buffer_seqs.resize(v_files_all2all.size(), make_pair(seq_t(), 0));
}

// ****************************************************************************
/*void load_tasks_one2all()
{
	load_tasks_all2all();
	v_files_one2all.swap(v_files_all2all);
}*/

// ****************************************************************************
void run_pairs_mode()
{
	for (int i = 0; i < params.no_threads; ++i)
	{
		v_threads.push_back(thread([&] {
			CWorker worker(params);

			while (true)
			{
				pair<string, string> task;
				{
					lock_guard<mutex> lck(mtx_queue);
					if (q_files_pairs.empty())
						return;
					task = q_files_pairs.front();
					q_files_pairs.pop();
				}

				CFatResults res;

				worker.clear();

				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				if (!worker.load_data(task.first, task.second))
				{
					lock_guard<mutex> lck(mtx_res);
					cout << task.first << " " << task.second << " - Error!" << endl;
					continue;
				}

				// 1 -> 2
				worker.prepare_kmers();
				worker.prepare_ht_short();
				worker.prepare_ht_long();

				worker.parse();
				//				worker.export_parsing();

				worker.calc_ani(res, 1);

				// 2 -> 1
				swap(task.first, task.second);
				worker.swap_data();

				worker.prepare_kmers();
				worker.prepare_ht_short();
				worker.prepare_ht_long();

				worker.parse();

				worker.calc_ani(res, 2);

				swap(task.first, task.second);

				res.ani[0] = (res.ani[1] + res.ani[2]) / 2;
				res.coverage[0] = (res.coverage[1] + res.coverage[2]) / 2;

				high_resolution_clock::time_point t2 = high_resolution_clock::now();

				res.time = duration_cast<duration<double>>(t2 - t1).count();

				{
					lock_guard<mutex> lck(mtx_res);
					m_results[task] = res;
					cout << task.first << " " << task.second <<
						" - ANI: " << 100 * res.ani[0] << " (" << 100 * res.ani[1] << " : " << 100 * res.ani[2] << ") " <<
						"   cov: " << res.coverage[0] << " (" << res.coverage[1] << " : " << res.coverage[2] << ") " <<
						"    time: " << res.time << endl;
				}
			}
			}));
	}

	for (auto& x : v_threads)
		x.join();

	FILE * f = fopen(output_file_name.c_str(), "w");
	FILE * g = fopen((output_file_name + ".csv").c_str(), "w");
	if (!f || !g)
	{
		cerr << "Cannot open " << output_file_name << endl;
		exit(0);
	}

	setvbuf(f, nullptr, _IOFBF, 32 << 20);
	setvbuf(g, nullptr, _IOFBF, 32 << 20);

	fprintf(g, "ref_name,query_name,ref_size,query_size,sym_in_matches1,sym_in_literals1,sym_in_matches2,sym_in_literals2,coverage,coverage1,coverage2,ani,ani1,ani2,time\n");

	for (auto& x : m_results)
	{
		fprintf(f, "%s %s : cov:%8.3f  ani:%8.3f\n", x.first.first.c_str(), x.first.second.c_str(), 100 * x.second.coverage[0], 100 * x.second.ani[0]);
		fprintf(g, "%s,%s,%d,%d,%d,%d,%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", x.first.first.c_str(), x.first.second.c_str(),
			x.second.ref_size, x.second.query_size,
			x.second.sym_in_matches[1], x.second.sym_in_literals[1], x.second.sym_in_matches[2], x.second.sym_in_literals[2],
			100 * x.second.coverage[0], 100 * x.second.coverage[1], 100 * x.second.coverage[2],
			100 * x.second.ani[0], 100 * x.second.ani[1], 100 * x.second.ani[2],
			x.second.time);
	}

	fclose(f);
	fclose(g);
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
	if (v_files_all2all.size() < 2)
	{
		cerr << "Too few input sequences\n";
		return;
	}

	CLZMatcher lz_matcher(params);

	lz_matcher.set_filter(filter_name, (uint32_t) filter_thr);
	lz_matcher.run_all2all(input_file_names, output_file_name);

	return;

	auto t_start = high_resolution_clock::now();

	vector<pair<int, string>> q_fn_data;
	map<pair<int, int>, CFatResults> p_results;
	vector<vector<pair<pair<int, int>, CFatResults>>> loc_results;
	atomic<int> a_fn_data;

	size_t q_size = min((size_t)16ull, v_files_all2all.size());

	parallel_queue<pair<int, CSharedWorker*>> q_to_prepare(q_size);
	parallel_queue<pair<int, CSharedWorker*>> q_ready(q_size);

	for (size_t i = 0; i < q_size; ++i)
		q_to_prepare.push(make_pair((int)i, new CSharedWorker(params)));

	cout << "All-2-All fast mode\n";	fflush(stdout);

	a_fn_data = 0;

	q_fn_data.clear();
	q_fn_data.reserve(v_files_all2all.size());
	for (int i = 0; i < (int)v_files_all2all.size(); ++i)
		q_fn_data.push_back(make_pair(i, v_files_all2all[i].first));

	loc_results.resize(params.no_threads);

	cout << "Prefetch input files\n";		fflush(stdout);
	atomic<int> fid = 0;

	vector<future<void>> v_fut;
	v_fut.reserve(params.no_threads);

	for (int i = 0; i < params.no_threads; ++i)
		v_fut.push_back(async(std::launch::async, [&] {
		CSharedWorker s_worker(params);

		while (true)
		{
			int cid = fid.fetch_add(1);
			if (fid >= v_files_all2all.size())
				break;

			if (!s_worker.load_data(v_files_all2all[cid].first, &(v_buffer_seqs[cid])))
			{
				lock_guard<mutex> lck(mtx_res);
				cout << "Cannot load: " << v_files_all2all[cid].first << endl;
				exit(0);
			}
		}
		}));

	for (auto& f : v_fut)
		f.get();


	v_files_all2all_order.clear();

	for (int i = 0; i < v_files_all2all.size(); ++i)
		v_files_all2all_order.emplace_back((int)v_buffer_seqs[i].first.size(), i);

	sort(v_files_all2all_order.begin(), v_files_all2all_order.end(), greater<pair<int, int>>());

	cerr << "ANI calculation\n";		fflush(stdout);

	barrier bar(params.no_threads + 1);

	thread thr_worker_base([&] {
		pair<int, CSharedWorker*> task;
		while (q_to_prepare.pop(task))
		{
			prepare_worker_base(task.second, task.first);
			q_ready.push(move(task));
		}

		q_ready.mark_completed();
		});

	pair<int, CSharedWorker*> ref_task;
	atomic<int> global_task_no = -1;

	thread thr_base_switcher([&] {
		for (int task_no = 0; task_no < v_files_all2all.size(); ++task_no)
		{
			q_ready.pop(ref_task);

			a_fn_data = 0;
			global_task_no = task_no;

			cerr << "Task " << task_no << "\r";	fflush(stdout);

			bar.arrive_and_wait();
			bar.arrive_and_wait();

			if (task_no + q_size < v_files_all2all.size())
				q_to_prepare.push(make_pair(task_no + q_size, ref_task.second));
			else
			{
				if (task_no + q_size == v_files_all2all.size())
					q_to_prepare.mark_completed();

				delete ref_task.second;
			}
		}

		global_task_no = v_files_all2all.size();

		bar.arrive_and_wait();
		});


	vector<thread> thr_workers;
	thr_workers.reserve(params.no_threads);

	for (int i = 0; i < params.no_threads; ++i)
	{
		thr_workers.emplace_back([&, i] {

		CSharedWorker s_worker(params);
		int thread_id = i;

		vector<pair<pair<int, int>, CFatResults>> res_loc;

		while (true)
		{
			bar.arrive_and_wait();

			if (global_task_no >= v_files_all2all.size())
				break;

			int task_no = global_task_no;

			s_worker.share_from(ref_task.second);

			int flt_id_task_no = filter_id_mapping[task_no];

			while (true)
			{
				pair<int, string> task;
				int cur_id = a_fn_data.fetch_add(1);

				if (cur_id >= (int)q_fn_data.size())
					break;

				cur_id = v_files_all2all_order[cur_id].second;

				task = q_fn_data[cur_id];

				int flt_id_task = filter_id_mapping[task.first];

				CFatResults res;

				high_resolution_clock::time_point t1 = high_resolution_clock::now();

				if (task_no == task.first)
				{
					res.ani[1] = 1;
					res.coverage[1] = 1;
					res.sym_in_literals[1] = 1;
					res.sym_in_matches[1] = 1;
					res.query_size = 1;
				}
				else if (!filter.empty() && filter.count(minmax(flt_id_task_no, flt_id_task)) == 0)
				{
					res.ani[1] = 0;
					res.coverage[1] = 0;
					res.sym_in_literals[1] = 0;
					res.sym_in_matches[1] = 0;
					res.query_size = 1;
				}
				else
				{
					s_worker.clear_data();

					if (!s_worker.load_data_fast(task.second, &(v_buffer_seqs[task.first])))
					{
						lock_guard<mutex> lck(mtx_res);
						cout << task.second << " - Error!" << endl;
						continue;
					}

					s_worker.prepare_kmers_data();
					s_worker.parse();

					s_worker.calc_ani(res, 1);

					high_resolution_clock::time_point t2 = high_resolution_clock::now();

					res.time = duration_cast<duration<double>>(t2 - t1).count();

					res_loc.emplace_back(make_pair(task_no, task.first), res);

					if (params.verbosity_level > 1)
					{
						cout << to_string(task_no) + " "s + to_string(task.first) +
							" - ANI: " + to_string(100 * res.ani[1]) +
							"   cov: " + to_string(res.coverage[1]) +
							"    time: " + to_string(res.time) + "\n";
					}
				}
			}

			bar.arrive_and_wait();
		}

		loc_results[thread_id] = move(res_loc);

		s_worker.share_from(nullptr);
			});
	}
		
	thr_worker_base.join();
	thr_base_switcher.join();
	for (auto& thr : thr_workers)
		thr.join();

	cerr << endl;

	cerr << "Merging local results\n";		fflush(stdout);

	size_t tot_loc_size = 0;
	for (auto& lr : loc_results)
		tot_loc_size += lr.size();

	//	p_results.reserve(tot_loc_size);

	for (auto& lr : loc_results)
	{
		p_results.insert(lr.begin(), lr.end());
		lr.clear();
		lr.shrink_to_fit();
	}

	auto t_stop = high_resolution_clock::now();

	cerr << "Processing time: " << duration<double>(t_stop - t_start).count() << " s" << endl;

	cerr << "Saving results\n";		fflush(stdout);

	FILE* fr1 = fopen((output_file_name + ".ani.csv").c_str(), "wb");
	FILE* fr2 = fopen((output_file_name + ".cov.csv").c_str(), "wb");
	FILE* fr3 = fopen((output_file_name + ".tani.csv").c_str(), "wb");

	if (!fr1 || !fr2 || !fr3)
	{
		cerr << "Cannot create res files\n";
		exit(0);
	}

	setvbuf(fr1, nullptr, _IOFBF, 32 << 20);
	setvbuf(fr2, nullptr, _IOFBF, 32 << 20);
	setvbuf(fr3, nullptr, _IOFBF, 32 << 20);

	fprintf(fr1, ",");
	fprintf(fr2, ",");
	fprintf(fr3, ",");

	for (int i = 0; i < (int)v_files_all2all.size(); ++i)
	{
		fprintf(fr1, "%s,", v_files_all2all[i].first.c_str());
		fprintf(fr2, "%s,", v_files_all2all[i].first.c_str());
		fprintf(fr3, "%s,", v_files_all2all[i].first.c_str());
	}

	fprintf(fr1, "\n");
	fprintf(fr2, "\n");
	fprintf(fr3, "\n");

	for (int task_no = 0; task_no < v_files_all2all.size(); ++task_no)
	{
		// Obliczanie total ANI z obu stron
		for (int i = 0; i < task_no; ++i)
		{
			auto p1 = make_pair(task_no, i);
			auto p2 = make_pair(i, task_no);

			auto res1 = p_results.find(p1);
			if (res1 != p_results.end())
			{
				auto res2 = p_results.find(p2);

				double len1 = res1->second.sym_in_matches[1];
				double q1 = res1->second.query_size;
				double len2 = res2->second.sym_in_matches[1];
				double q2 = res2->second.query_size;

				res1->second.total_ani = (len1 + len2) / (q1 + q2);

/*				if (q1 < q2)
					res1->second.total_ani = len1 / q1;
				else
					res1->second.total_ani = len2 / q2;*/
				if (res1->second.total_ani > 1)
				{
					cerr << v_files_all2all[task_no].first << " : " << v_files_all2all[i].first << res1->second.total_ani << endl;
					res1->second.total_ani = 1;
				}
			}
		}

		fprintf(fr1, "%s,", v_files_all2all[task_no].first.c_str());
		fprintf(fr2, "%s,", v_files_all2all[task_no].first.c_str());
		fprintf(fr3, "%s,", v_files_all2all[task_no].first.c_str());

		for (int i = 0; i < (int)v_files_all2all.size(); ++i)
		{
			auto p = p_results.find(make_pair(task_no, i));
			if (p != p_results.end())
			{
				fprintf(fr1, "%.5f,", p->second.ani[1]);
				fprintf(fr2, "%.5f,", p->second.coverage[1]);
				//			fprintf(fr3, "%.5f,", p_results[make_pair(task_no, i)].ani[1] * p_results[make_pair(task_no, i)].coverage[1]);
				fprintf(fr3, "%.5f,", p->second.total_ani);
			}
			else
			{
				//				putc('0', fr1);				putc('0', fr2);				putc('0', fr3);
				//				putc(',', fr1);				putc(',', fr2);				putc(',', fr3);
				fputs("0.00000,", fr1);
				fputs("0.00000,", fr2);
				fputs("0.00000,", fr3);
			}
		}

		fprintf(fr1, "\n");
		fprintf(fr2, "\n");
		fprintf(fr3, "\n");
	}

	fclose(fr1);
	fclose(fr2);
	fclose(fr3);


	FILE* f = fopen(output_file_name.c_str(), "w");
	FILE* g = fopen((output_file_name + ".csv").c_str(), "w");
	if (!f || !g)
	{
		cerr << "Cannot open " << output_file_name << endl;
		exit(0);
	}

	setvbuf(f, nullptr, _IOFBF, 32 << 20);
	setvbuf(g, nullptr, _IOFBF, 32 << 20);

	fprintf(g, "ref_name,query_name,ref_size,query_size,sym_in_matches1,sym_in_literals1,coverage,ani,time\n");

	for (auto& x : p_results)
	{
		fprintf(f, "%s %s : cov:%8.3f  ani:%8.3f\n", v_files_all2all[x.first.first].first.c_str(), v_files_all2all[x.first.second].first.c_str(), 100 * x.second.coverage[1], 100 * x.second.ani[1]);
		fprintf(g, "%s,%s,%d,%d,%d,%d,%.5f,%.5f,%.5f\n", v_files_all2all[x.first.first].first.c_str(), v_files_all2all[x.first.second].first.c_str(),
			x.second.ref_size, x.second.query_size,
			x.second.sym_in_matches[1], x.second.sym_in_literals[1],
			100 * x.second.coverage[1],
			100 * x.second.ani[1],
			x.second.time);
	}

	fclose(f);
	fclose(g);
}

// ****************************************************************************
/*void run_one2all_mode()
{
	// Prepare one-2-all pairs
	for (auto i = 0; i < (int) v_files_one2all.size(); ++i)
		q_files_pairs.push(make_pair(v_files_one2all[0], v_files_one2all[i]));

	for (int i = 0; i < no_threads; ++i)
	{
		v_threads.push_back(thread([&] {
			CWorker worker;

			while (true)
			{
				pair<string, string> task;
				{
					lock_guard<mutex> lck(mtx_queue);
					if (q_files_pairs.empty())
						return;
					task = q_files_pairs.front();
					q_files_pairs.pop();
				}

				CFatResults res;

				worker.clear();

				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				if (!worker.load_data(task.first, task.second))
				{
					lock_guard<mutex> lck(mtx_res);
					cout << task.first << " " << task.second << " - Error!" << endl;
					continue;
				}

				// 1 -> 2
				worker.prepare_kmers();
				worker.prepare_ht_short();
				worker.prepare_ht_long();

				worker.parse();
				//				worker.export_parsing();

				worker.calc_ani(res, 1);

				// 2 -> 1
				swap(task.first, task.second);
				worker.swap_data();

				worker.prepare_kmers();
				worker.prepare_ht_short();
				worker.prepare_ht_long();

				worker.parse();

				worker.calc_ani(res, 2);

				swap(task.first, task.second);

				res.ani[0] = (res.ani[1] + res.ani[2]) / 2;
				res.coverage[0] = (res.coverage[1] + res.coverage[2]) / 2;

				high_resolution_clock::time_point t2 = high_resolution_clock::now();

				res.time = duration_cast<duration<double>>(t2 - t1).count();

				{
					lock_guard<mutex> lck(mtx_res);
					m_results[task] = res;

					cout << task.first << " " << task.second <<
						" - ANI: " << 100 * res.ani[0] << " (" << 100 * res.ani[1] << " : " << 100 * res.ani[2] << ") " <<
						"   cov: " << res.coverage[0] << " (" << res.coverage[1] << " : " << res.coverage[2] << ") " <<
						"    time: " << res.time << endl;
				}
			}
			}));
	}

	for (auto& x : v_threads)
		x.join();


	FILE* fr1 = fopen((output_name + "1.ani.csv").c_str(), "wb");
	FILE* fr2 = fopen((output_name + "1.cov.csv").c_str(), "wb");
	FILE* fr3 = fopen((output_name + "2.ani.csv").c_str(), "wb");
	FILE* fr4 = fopen((output_name + "2.cov.csv").c_str(), "wb");

	if (!fr1 || !fr2 || !fr3 || !fr4)
	{
		cerr << "Cannot create res files\n";
		exit(0);
	}

	// 0th -> all
	fprintf(fr1, ",");
	fprintf(fr2, ",");

	for (int i = 0; i < (int)v_files_one2all.size(); ++i)
	{
		fprintf(fr1, "%s,", v_files_one2all[i].c_str());
		fprintf(fr2, "%s,", v_files_one2all[i].c_str());
	}

	fprintf(fr1, "\n");
	fprintf(fr2, "\n");

	fprintf(fr1, "%s,", v_files_one2all[0].c_str());
	fprintf(fr2, "%s,", v_files_one2all[0].c_str());
	
	for (int i = 0; i < (int)v_files_one2all.size(); ++i)
	{
		fprintf(fr1, "%.5f,", m_results[make_pair(v_files_one2all[0], v_files_one2all[i])].ani[1]);
		fprintf(fr2, "%.5f,", m_results[make_pair(v_files_one2all[0], v_files_one2all[i])].coverage[1]);
	}

	fprintf(fr1, "\n");
	fprintf(fr2, "\n");

	// all -> 0th
	fprintf(fr3, ",");
	fprintf(fr4, ",");
	fprintf(fr3, "%s,", v_files_one2all[0].c_str());
	fprintf(fr4, "%s,", v_files_one2all[0].c_str());
	fprintf(fr3, "\n");
	fprintf(fr4, "\n");

	for (int i = 0; i < (int)v_files_one2all.size(); ++i)
	{
		fprintf(fr3, "%s,", v_files_one2all[i].c_str());
		fprintf(fr4, "%s,", v_files_one2all[i].c_str());
		fprintf(fr3, "%.5f,", m_results[make_pair(v_files_one2all[0], v_files_one2all[i])].ani[2]);
		fprintf(fr4, "%.5f,", m_results[make_pair(v_files_one2all[0], v_files_one2all[i])].coverage[2]);
		fprintf(fr3, "\n");
		fprintf(fr4, "\n");
	}

	fclose(fr1);
	fclose(fr2);
	fclose(fr3);
	fclose(fr4);
}
*/
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
		if (!load_tasks_all2all())
			return 0;
		run_all2all_threads_mode();
	}
/*	else if (is_one2all)
	{
		load_tasks_one2all();
		run_one2all_mode();
	}*/
	else
	{
/*		load_tasks_pairs();
		run_pairs_mode();*/
	}

	return 0;
}
