// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

#include "defs.h"
#include "worker.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <queue>
#include <map>
#include <thread>
#include <mutex>

using namespace std::chrono;

queue<pair<string, string>> q_files;
map<pair<string, string>, CResults> m_results;

mutex mtx_queue;
mutex mtx_res;

vector<thread> v_threads;
int no_threads = 0;
string output_name;

void load_tasks(int argc, char **argv);

// ****************************************************************************
void load_tasks(int argc, char **argv)
{
	if (argc < 3)
	{
		cerr << "Usage: ani-entropy <file_list> <output_name> [no_threads]\n";
		exit(0);
	}
	
	FILE *f = fopen(argv[1], "rb");
	if (!f)
	{
		cerr << "Error: Cannot load " + string(argv[1]) + "\n";
		exit(0);
	}
	
	while (!feof(f))
	{
		char s[256], t[256];
		fscanf(f, "%s%s", s, t);
		if (feof(f))
			break;

		q_files.push(make_pair(string(s), string(t)));
	}

	fclose(f);

	output_name = string(argv[2]);

	if(argc == 4)
		no_threads = atoi(argv[3]);
}

// ****************************************************************************
int main(int argc, char **argv)
{
	load_tasks(argc, argv);

	if (no_threads == 0)
	{
		no_threads = thread::hardware_concurrency();
		if (!no_threads)
			no_threads = 1;
	}

	for (int i = 0; i < no_threads; ++i)
	{
		v_threads.push_back(thread([&] {
			CWorker worker;

			while (true)
			{
				pair<string, string> task;
				{
					lock_guard<mutex> lck(mtx_queue);
					if (q_files.empty())
						return;
					task = q_files.front();
					q_files.pop();
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
				worker.prepare_ht_short();
				worker.prepare_ht_long();

				worker.prepare_pf();
				worker.parse_new();
//				worker.parsing_postprocess();
				worker.export_parsing();

				worker.calc_ani(res, 1);

				// 2 -> 1
				swap(task.first, task.second);
				worker.swap_data();
				worker.prepare_ht_short();
				worker.prepare_ht_long();

				worker.prepare_pf();
				worker.parse_new();
//				worker.parsing_postprocess();

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

	for (auto &x : v_threads)
		x.join();

	FILE *f = fopen(output_name.c_str(), "w");
	FILE *g = fopen((output_name+".csv").c_str(), "w");
	if (!f || !g)
	{
		cerr << "Cannot open " << output_name << endl;
		exit(0);
	}

	setvbuf(f, nullptr, _IOFBF, 32 << 20);
	setvbuf(g, nullptr, _IOFBF, 32 << 20);

	fprintf(g, "ref_name,query_name,ref_size,query_size,sym_in_matches1,sym_in_literals1,sym_in_matches2,sym_in_literals2,coverage,coverage1,coverage2,ani,ani1,ani2,time\n");

	for (auto &x : m_results)
	{
		fprintf(f, "%s %s : cov:%8.3f  ani:%8.3f\n", x.first.first.c_str(), x.first.second.c_str(), 100 * x.second.coverage[0], 100 * x.second.ani[0]);
		fprintf(g, "%s,%s,%d,%d,%d,%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", x.first.first.c_str(), x.first.second.c_str(),
			x.second.ref_size, x.second.query_size, 
			x.second.sym_in_matches[1], x.second.sym_in_literals[1], x.second.sym_in_matches[2], x.second.sym_in_literals[2],
			100 * x.second.coverage[0], 100 * x.second.coverage[1], 100 * x.second.coverage[2],
			100 * x.second.ani[0], 100 * x.second.ani[1], 100 * x.second.ani[2],
			x.second.time);
	}

	fclose(f);
	fclose(g);

	return 0;
}
