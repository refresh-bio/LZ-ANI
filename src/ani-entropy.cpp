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

				worker.clear();

				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				if (!worker.load_data(task.first, task.second))
				{
					lock_guard<mutex> lck(mtx_res);
					cout << task.first << " " << task.second << " - Error!" << endl;
					continue;
				}

				worker.prepare_ht_short();
				worker.prepare_ht_long();

				worker.prepare_pf();
				worker.parse();
				worker.parsing_postprocess();
//				worker.export_parsing();

				CResults res;
				worker.calc_ani(res);

				high_resolution_clock::time_point t2 = high_resolution_clock::now();

				res.time = duration_cast<duration<double>>(t2 - t1).count();

				{
					lock_guard<mutex> lck(mtx_res);
					m_results[task] = res;
					cout << task.first << " " << task.second << " - ANI: " << 100 * res.ani << "   cov: " << res.coverage << "    time: " << res.time << endl;
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

	fprintf(g, "ref_name,query_name,ref_size,query_size,sym_in_matches,sym_in_literals,coverage,ani,time\n");

	for (auto &x : m_results)
	{
		fprintf(f, "%s %s : cov:%8.3f  ani:%8.3f\n", x.first.first.c_str(), x.first.second.c_str(), 100 * x.second.coverage, 100 * x.second.ani);
		fprintf(g, "%s,%s,%d,%d,%d,%d,%.3f,%.3f,%.3f\n", x.first.first.c_str(), x.first.second.c_str(),
			x.second.ref_size, x.second.query_size, x.second.sym_in_matches, x.second.sym_in_literals, 100 * x.second.coverage, 100 * x.second.ani, x.second.time);
	}

	fclose(f);
	fclose(g);

	return 0;
}
