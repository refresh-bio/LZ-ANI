// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "defs.h"
#include "worker.h"
#include "s_worker.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <cstdio>
#include <future>

using namespace std::chrono;

queue<pair<string, string>> q_files_pairs;
vector<string> v_files_all2all;
vector<string> v_files_one2all;
map<pair<string, string>, CResults> m_results;

mutex mtx_queue;
mutex mtx_res;

vector<pair<seq_t, int>> v_buffer_seqs;
vector<thread> v_threads;
int no_threads = 0;
string input_name;
string output_name;
bool is_all2all = false;
bool is_one2all = false;
bool buffer_data = false;

int verbosity_level = 1;

int MIN_MATCH_LEN = DEF_MIN_MATCH_LEN;
int MIN_CLOSE_MATCH_LEN = DEF_MIN_CLOSE_MATCH_LEN;
int MIN_DISTANT_MATCH_LEN = DEF_MIN_DISTANT_MATCH_LEN;
int CLOSE_DIST = DEF_CLOSE_DIST;
//int LONG_LITERAL_RUN_LEN = DEF_LONG_LITERAL_RUN_LEN;
int MAX_LIT_RUN_IN_MATCH = DEF_MAX_LIT_RUN_IN_MATCH;
double MIN_COVERAGE = DEF_MIN_COVERAGE;
int MIN_REGION_LEN = DEF_MIN_REGION_LEN;
int APPROX_WINDOW = DEF_APPROX_WINDOW;
int APPROX_MISMATCHES = DEF_APPROX_MISMATCHES;
int APPROX_RUNLEN = DEF_APPROX_RUNLEN;
int RANGE_FROM = DEF_RANGE_FROM;
int RANGE_TO = DEF_RANGE_TO;

void load_params(int argc, char** argv);
void load_tasks_pairs();
void load_tasks_all2all();
void load_tasks_one2all();
void usage();
void run_pairs_mode();
void run_one2all_mode();

// ****************************************************************************
void usage()
{
	cerr << "lz-ani [options] <in_list> <output_file>\n";
	cerr << "Options:\n";
	cerr << "   -a2a         - turn on all to all mode (default: " << is_all2all << ")\n";
	cerr << "   -o2a         - turn on one to all mode (default: " << is_one2all << ")\n";
	cerr << "   -bs          - turn on all sequences buffering  (default: " << buffer_data << ")\n";
	cerr << "   -t <val>     - no of threads (default: " << no_threads << ")\n";
	cerr << "   -mml <val>   - min. match length (default: " << MIN_MATCH_LEN << ")\n";
	cerr << "   -mdl <val>   - min. distant length (default: " << MIN_DISTANT_MATCH_LEN << ")\n";
	cerr << "   -cd <val>    - max. dist. between close matches (default: " << CLOSE_DIST << ")\n";
	cerr << "   -mlrin <val> - max. literal run len. in match (dafault: " << MAX_LIT_RUN_IN_MATCH << ")\n";
	cerr << "   -cov <val>   - min. coverage threshold (default: " << MIN_COVERAGE << ")\n";
	cerr << "   -reg <val>   - min. considered region length (default: " << MIN_REGION_LEN << ")\n";
	cerr <<	"   -aw <val>    - approx. window length (default: " << APPROX_WINDOW << ")\n";
	cerr << "   -am <val>    - max. no. of mismatches in approx. window (default: " << APPROX_MISMATCHES << ")\n";
	cerr << "   -ar <val>    - min. length of run ending approx. extension (default: " << APPROX_RUNLEN << ")\n";
	cerr << "   -rng-f <val> - start of range of sequences to compare with all other (default: " << RANGE_FROM << ")\n";
	cerr << "   -rng-t <val> - end of range of sequences to compare with all other (default: " << RANGE_TO << ")\n";
}

void load_params(int argc, char** argv)
{
	if (argc < 3)
	{
		usage();
		exit(0);
	}

	input_name = string(argv[argc - 2]);
	output_name = string(argv[argc - 1]);

	for (int i = 1; i < argc - 2;)
	{
		string par = string(argv[i]);

		if (par == "-a2a")
		{
			is_all2all = true;
			is_one2all = false;
			i += 1;
		}
		else if (par == "-o2a")
		{
			is_one2all = true;
			is_all2all = false;
			i += 1;
		}
		else if (par == "-bs")
		{
			buffer_data = true;
			i += 1;
		}
		else if (par == "-t")
		{
			no_threads = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-mml")
		{
			MIN_MATCH_LEN = atoi(argv[i + 1]);
			MIN_CLOSE_MATCH_LEN = MIN_MATCH_LEN;
			i += 2;
		}
		else if (par == "-mdl")
		{
			MIN_DISTANT_MATCH_LEN = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-cd")
		{
			CLOSE_DIST = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-mlrim")
		{
			MAX_LIT_RUN_IN_MATCH = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-cov")
		{
			MIN_COVERAGE = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-reg")
		{
			MIN_REGION_LEN = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-aw")
		{
			APPROX_WINDOW = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-am")
		{
			APPROX_MISMATCHES = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-ar")
		{
			APPROX_RUNLEN = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-rng-f")
		{
			RANGE_FROM = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-rng-t")
		{
			RANGE_TO = atoi(argv[i + 1]);
			i += 2;
		}
		else
		{
			cerr << "Unknown parameter: " << string(argv[i]) << endl;
			usage();
			exit(0);
		}
	}
}

// ****************************************************************************
void load_tasks_pairs()
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
}

// ****************************************************************************
void load_tasks_all2all()
{
	FILE* f = fopen(input_name.c_str(), "rb");
	if (!f)
	{
		cerr << "Error: Cannot load " + input_name + "\n";
		exit(0);
	}

	while (!feof(f))
	{
		char s[256];
		if (fscanf(f, "%s", s) == EOF)
			break;
		if (feof(f))
			break;

		v_files_all2all.push_back(string(s));
	}

	if (buffer_data)
		v_buffer_seqs.resize(v_files_all2all.size(), make_pair(seq_t(), 0));

	fclose(f);
}

// ****************************************************************************
void load_tasks_one2all()
{
	load_tasks_all2all();
	v_files_one2all.swap(v_files_all2all);
}
// ****************************************************************************
void run_pairs_mode()
{
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

				CResults res;

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

	FILE * f = fopen(output_name.c_str(), "w");
	FILE * g = fopen((output_name + ".csv").c_str(), "w");
	if (!f || !g)
	{
		cerr << "Cannot open " << output_name << endl;
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
void run_all2all_mode()
{
	CSharedWorker *s_worker_base = new CSharedWorker;
	queue<pair<int, string>> q_fn_data;
	map<pair<int, int>, CResults> p_results;

	FILE* fr1 = fopen((output_name + ".ani.csv").c_str(), "wb");
	FILE* fr2 = fopen((output_name + ".cov.csv").c_str(), "wb");
	FILE* fr3 = fopen((output_name + ".tani.csv").c_str(), "wb");

	if (!fr1 || !fr2 || !fr3)
	{
		cerr << "Cannot create res files\n";
		exit(0);
	}

	fprintf(fr1, ",");
	fprintf(fr2, ",");
	fprintf(fr3, ",");

	for (int i = 0; i < (int) v_files_all2all.size(); ++i)
	{
		fprintf(fr1, "%s,", v_files_all2all[i].c_str());
		fprintf(fr2, "%s,", v_files_all2all[i].c_str());
		fprintf(fr3, "%s,", v_files_all2all[i].c_str());
	}

	fprintf(fr1, "\n");
	fprintf(fr2, "\n");
	fprintf(fr3, "\n");

	fclose(fr1);
	fclose(fr2);
	fclose(fr3);

	cout << "All-2-All mode\n";	fflush(stdout);

//	for (int task_no = 0; task_no < v_files_all2all.size(); ++task_no)
	for (int task_no = RANGE_FROM; task_no < min(RANGE_TO+1, (int) v_files_all2all.size()); ++task_no)
	{
		cout << "Task " << task_no << endl;	fflush(stdout);
		s_worker_base->clear_ref();

		if (!s_worker_base->load_reference(v_files_all2all[task_no], buffer_data ? &(v_buffer_seqs[task_no]) : nullptr))
		{
			cerr << "Cannot read: " << v_files_all2all[task_no] << endl;
			continue;
		}

		for (int i = 0; i < (int) v_files_all2all.size(); ++i)
			q_fn_data.push(make_pair(i, v_files_all2all[i]));
				
		std::future<void> fut = std::async(std::launch::async, [&] {
			s_worker_base->prepare_kmers_ref_short();
			s_worker_base->prepare_ht_short(); 
			});

		s_worker_base->prepare_kmers_ref_long();
		s_worker_base->prepare_ht_long();
		fut.get();

		v_threads.clear();
		for (int i = 0; i < no_threads; ++i)
		{
			v_threads.push_back(thread([&] {
				CSharedWorker s_worker;

				s_worker.share_from(s_worker_base);

				while (true)
				{
					pair<int, string> task;
					{
						lock_guard<mutex> lck(mtx_queue);
						if (q_fn_data.empty())
						{
							s_worker.share_from(nullptr);
							return;
						}
						task = q_fn_data.front();
						q_fn_data.pop();
					}

					CResults res;

					s_worker.clear_data();

					high_resolution_clock::time_point t1 = high_resolution_clock::now();
					if (!s_worker.load_data(task.second, buffer_data ? &(v_buffer_seqs[task.first]) : nullptr))
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
					
					{
						lock_guard<mutex> lck(mtx_res);
						
						p_results[make_pair(task_no, task.first)] = res;

						if(verbosity_level > 1)
							cout << to_string(task_no) + " "s + to_string(task.first) +
								" - ANI: " + to_string(100 * res.ani[1]) +
								"   cov: "  + to_string(res.coverage[1]) + 
								"    time: " + to_string(res.time) + "\n";
					}
				}

				s_worker.share_from(nullptr);

				}));
		}

		for (auto& x : v_threads)
			x.join();
			
		FILE* fr1 = fopen((output_name + ".ani.csv").c_str(), "ab");
		FILE* fr2 = fopen((output_name + ".cov.csv").c_str(), "ab");
		FILE* fr3 = fopen((output_name + ".tani.csv").c_str(), "ab");

		setvbuf(fr1, nullptr, _IOFBF, 1 << 20);
		setvbuf(fr2, nullptr, _IOFBF, 1 << 20);
		setvbuf(fr3, nullptr, _IOFBF, 1 << 20);

		fprintf(fr1, "%s,", v_files_all2all[task_no].c_str());
		fprintf(fr2, "%s,", v_files_all2all[task_no].c_str());
		fprintf(fr3, "%s,", v_files_all2all[task_no].c_str());

		for (int i = 0; i < (int)v_files_all2all.size(); ++i)
		{
			fprintf(fr1, "%.5f,", p_results[make_pair(task_no, i)].ani[1]);
			fprintf(fr2, "%.5f,", p_results[make_pair(task_no, i)].coverage[1]);
			fprintf(fr3, "%.5f,", p_results[make_pair(task_no, i)].ani[1] * p_results[make_pair(task_no, i)].coverage[1]);
		}

		fprintf(fr1, "\n");
		fprintf(fr2, "\n");
		fprintf(fr3, "\n");

		fclose(fr1);
		fclose(fr2);
		fclose(fr3);
	}

	delete s_worker_base;

	FILE* f = fopen(output_name.c_str(), "w");
	FILE* g = fopen((output_name + ".csv").c_str(), "w");
	if (!f || !g)
	{
		cerr << "Cannot open " << output_name << endl;
		exit(0);
	}

	setvbuf(f, nullptr, _IOFBF, 32 << 20);
	setvbuf(g, nullptr, _IOFBF, 32 << 20);

	fprintf(g, "ref_name,query_name,ref_size,query_size,sym_in_matches1,sym_in_literals1,coverage,ani,time\n");

	for (auto& x : p_results)
	{
		fprintf(f, "%s %s : cov:%8.3f  ani:%8.3f\n", v_files_all2all[x.first.first].c_str(), v_files_all2all[x.first.second].c_str(), 100 * x.second.coverage[1], 100 * x.second.ani[1]);
		fprintf(g, "%s,%s,%d,%d,%d,%d,%.5f,%.5f,%.5f\n", v_files_all2all[x.first.first].c_str(), v_files_all2all[x.first.second].c_str(),
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
void run_one2all_mode()
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

				CResults res;

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

// ****************************************************************************
int main(int argc, char **argv)
{
	load_params(argc, argv);

	if (no_threads == 0)
	{
		no_threads = thread::hardware_concurrency();
		if (!no_threads)
			no_threads = 1;
	}

	if (is_all2all)
	{
		load_tasks_all2all();
		run_all2all_mode();
	}
	else if (is_one2all)
	{
		load_tasks_one2all();
		run_one2all_mode();
	}
	else
	{
		load_tasks_pairs();
		run_pairs_mode();
	}

	return 0;
}
