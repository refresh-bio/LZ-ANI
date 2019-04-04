// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

#include "defs.h"
#include "files.h"
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

				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				worker.load_data(task.first, task.second);
				worker.prepare_ht();
				worker.parse();
				worker.parsing_postprocess();
//				worker.export_parsing();

				CResults res;
				worker.calc_ani(res);

				high_resolution_clock::time_point t2 = high_resolution_clock::now();

				duration<double> time_span = duration_cast<duration<double>>(t2 - t1);


			}
		}));
	}

	for (auto &x : v_threads)
		x.join();

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

/*	cout << "*** Loading...";
	if (!load_data(argc, argv))
	{
		return 0;
	}
	cout << " ok\n";

	cout << "*** HT preparation...";
	prepare_ht();
	cout << " ok\n";

	cout << "*** Parsing...";
	parse();
	parsing_postprocess();
	export_parsing();

	cout << " ok\n";

	calc_ani();*/

	cout << "*** End of work\n";

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Running time: " << time_span.count() << " seconds.\n";

//	getchar();

	return 0;
}
