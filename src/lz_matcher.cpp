#include <iostream>
#include <istream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <utility>
#include <atomic>
#include <mutex>
#include <barrier>

#include "parallel-queues.h"
#include "lz_matcher.h"
#include "utils.h"
#include "s_worker.h"

using namespace refresh;


// ****************************************************************************
bool CLZMatcher::run_all2all(const string &output_file_name)
{
	times.clear();
	times.emplace_back(high_resolution_clock::now(), "");

//	input_file_names = _input_file_names;

	if (!reorder_input_files())
		return false;

	times.emplace_back(high_resolution_clock::now(), "Reordering input files");

	prefetch_input_files();
	times.emplace_back(high_resolution_clock::now(), "Prefetching input files");

	results.clear();


	barrier bar(params.no_threads + 1);
	vector<vector<pair<pair_id_t, CResults>>> loc_results(params.no_threads);
	atomic<size_t> a_fn_data;

	// Initailize queues for maintaining reference data
	size_t q_size = min((size_t)16ull, input_file_desc.size());

	parallel_queue<pair<size_t, CSharedWorker*>> q_to_prepare(q_size);
	parallel_queue<pair<size_t, CSharedWorker*>> q_ready(q_size);

	for (size_t i = 0; i < q_size; ++i)
		q_to_prepare.push(make_pair(i, new CSharedWorker(params, data_storage)));


	// Start thread preparing reference sequences
	thread thr_worker_base([&] {
		pair<size_t, CSharedWorker*> task;
		while (q_to_prepare.pop(task))
		{
			prepare_worker_base(task.second, task.first);
			q_ready.push(move(task));
		}

		q_ready.mark_completed();
		});

	pair<size_t, CSharedWorker*> ref_task;
	atomic<size_t> global_task_no = 0;

	// Start thread for reference switching
	thread thr_base_switcher([&] {
		for (size_t task_no = 0; task_no < input_file_desc.size(); ++task_no)
		{
			q_ready.pop(ref_task);

			a_fn_data = 0;
			global_task_no = task_no;
			
			if(params.verbosity_level > 1)
			{
				cerr << "Task " << task_no << "\r";	
				fflush(stdout);
			}

			bar.arrive_and_wait();
			bar.arrive_and_wait();

			if (task_no + q_size < input_file_desc.size())
				q_to_prepare.push(make_pair(task_no + q_size, ref_task.second));
			else
			{
				if (task_no + q_size == input_file_desc.size())
					q_to_prepare.mark_completed();

				delete ref_task.second;
			}
		}

		global_task_no = input_file_desc.size();

		bar.arrive_and_wait();
		});


	// Start worker threads
	vector<thread> thr_workers;
	thr_workers.reserve(params.no_threads);

	for (int i = 0; i < params.no_threads; ++i)
	{
		thr_workers.emplace_back([&, i] {
			CSharedWorker s_worker(params, data_storage);
			int thread_id = i;

			vector<pair<pair_id_t, CResults>> res_loc;

			while (true)
			{
				bar.arrive_and_wait();

				if (global_task_no >= input_file_desc.size())
					break;

				size_t ref_id = global_task_no;

				s_worker.share_from(ref_task.second);

				while (true)
				{
					pair<size_t, string> task;
					size_t cur_id = a_fn_data.fetch_add(1);

					if (cur_id >= input_file_desc.size())
						break;

					if (!filter_set.empty())
					{
						if (filter_set.count(encode_pair_id_mm(ref_id, cur_id)) == 0)
							continue;
					}

					if (ref_id == cur_id)
						res_loc.emplace_back(make_pair(encode_pair_id(ref_id, cur_id), CResults(1, 0, 1)));
					else if (!filter_set.empty() && filter_set.count(encode_pair_id_mm(ref_id, cur_id)) == 0)
					{
						res_loc.emplace_back(make_pair(encode_pair_id(ref_id, cur_id), CResults(0, 0, 0)));
					}
					else
					{
						s_worker.clear_data();

						if (!s_worker.load_data_fast(input_file_desc[cur_id]))
						{
							// !!! TODO: Set stop token
							continue;
						}

						s_worker.prepare_kmers_data();
						s_worker.parse();

						auto results = s_worker.calc_stats();

						res_loc.emplace_back(make_pair(encode_pair_id(ref_id, cur_id), results));
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

	for (auto& lr : loc_results)
	{
		results.insert(lr.begin(), lr.end());
		lr.clear();
		lr.shrink_to_fit();
	}

	times.emplace_back(high_resolution_clock::now(), "LZ matching");

	store_results(output_file_name);

	show_timinigs_info();

	return true;
}

// ****************************************************************************
bool CLZMatcher::store_results(const string &output_file_name)
{
	ofstream ofs(output_file_name);

	if (!ofs.is_open())
	{
		cerr << "Cannot open output file: " << output_file_name << endl;
		return false;
	}

	// Store params
	ofs << params.str();

	// Store file desc
	ofs << "[no_input_sequences]" << endl;
	ofs << input_file_desc.size() << endl;

	ofs << "[input_sequences]" << endl;
	for (const auto& x : input_file_desc)
		ofs << x.file_name << " " << x.seq_name << " " << x.seq_size - x.n_parts * params.close_dist << " " << x.n_parts << endl;

	// Store mapping results in sparse form
	uint32_t q1, q2;

	ofs << "[lz_similarities]" << endl;
	if (params.output_mode == output_mode_t::sym_in_matches_literals)
	{
		for (const auto& x : results)
		{
			if (!is_pair_id_mm(x.first))
				continue;

			auto r = results.find(swap_pair_id(x.first));
			tie(q1, q2) = decode_pair_id(x.first);

			if (q1 != q2)
				ofs << to_string(q1) + " " + to_string(q2) + " " +
				to_string(r->second.sym_in_matches) + " " + to_string(r->second.sym_in_literals) + " " + to_string(r->second.no_components) +
				" " +
				to_string(x.second.sym_in_matches) + " " + to_string(x.second.sym_in_literals) + " " + to_string(x.second.no_components) +
				"\n";
		}
	}
	else if (params.output_mode == output_mode_t::total_ani)
	{
		for (const auto& x : results)
		{
			if (!is_pair_id_mm(x.first))
				continue;

			auto r = results.find(swap_pair_id(x.first));
			tie(q1, q2) = decode_pair_id(x.first);

			double total_ani = (x.second.sym_in_matches + r->second.sym_in_matches) /
				(input_file_desc[q1].seq_size - input_file_desc[q1].n_parts * params.close_dist + input_file_desc[q2].seq_size - input_file_desc[q2].n_parts * params.close_dist);

			if (q1 != q2)
				ofs << to_string(q1) + " " + to_string(q2) + " " + to_string(total_ani) + "\n";
		}
	}

	return true;
}

// ****************************************************************************
bool CLZMatcher::reorder_input_files()
{
	if (params.verbosity_level > 1)
	{
		cerr << "Reordering input files according to size\n";
	}

	if (!filter_set.empty() && filter_set.size() != data_storage.size())
	{
		cerr << "Sizes of filter and input file names are different!\n";
		return false;
	}

	input_file_desc.clear();

	// Check files sizes
	for (const auto& ds_item : data_storage)
		input_file_desc.emplace_back(ds_item.name, remove_path_from_file(ds_item.name), ds_item.size, 0, 0);

	// Sort files from the largest one - just for better parallelization of calculations
	stable_sort(input_file_desc.begin(), input_file_desc.end(), [](const auto& x, const auto& y)
		{
			if (x.file_size != y.file_size)
				return x.file_size > y.file_size;
			return x.seq_name < y.seq_name;
		});

	if (filter_vec.empty())
		return true;

	// Check presence of sequences in filter
	unordered_map<string, int> m_filter_names;
	for (size_t i = 0; i < filter_genome_names.size(); ++i)
		m_filter_names[filter_genome_names[i].first] = i;

	for (size_t i = 0; i < input_file_desc.size(); ++i)
	{
		auto p = m_filter_names.find(input_file_desc[i].seq_name);

		if (p == m_filter_names.end())
		{
			cerr << "File name not present in filter: " << input_file_desc[i].seq_name << endl;
			return false;
		}

		filter_genome_names[p->second].second = i;
	}

	// Renumerate filter set
	filter_set.clear();
	uint32_t g1, g2;

	for (auto x : filter_vec)
	{
		tie(g1, g2) = decode_pair_id(x);

		// Assign new ids (according to sorted ordering)
		g1 = filter_genome_names[g1].second;
		g2 = filter_genome_names[g2].second;

		filter_set.insert(encode_pair_id_mm(g1, g2));
	}

	filter_vec.clear();
	filter_vec.shrink_to_fit();

	return true;
}

// ****************************************************************************
bool CLZMatcher::set_filter(const string& _filter_name, const double _filter_thr)
{
	filter_name = _filter_name;
	filter_thr = _filter_thr;

	return load_filter();
}

// ****************************************************************************
bool CLZMatcher::load_filter()
{
	ifstream ifs(filter_name);

	if (!ifs.is_open())
	{
		cerr << "Cannot open filter file: " << filter_name << endl;
		return false;
	}

	string line;
	vector<string> parts;
	vector<string> elem;

	getline(ifs, line);		// genome names
	auto vec = split(line, ',');
	if (vec.size() <= 2)
	{
		cerr << "Incorrect kmer-db filter file\n";
		return false;
	}

	// Load filter genome names with mappings to the lz_matcher input order
/*	filter_genome_names.resize(vec.size() - 2);
	for (size_t i = 0; i < vec.size() - 2; ++i)
		filter_genome_names[i] = make_pair(strip_at_space(vec[i + 2]), -1);		// currently the matching is unknown
		*/

	filter_genome_names.resize(vec.size() - 1);
	for (size_t i = 0; i < vec.size() - 1; ++i)
		filter_genome_names[i] = make_pair(strip_at_space(vec[i + 1]), -1);		// currently the matching is unknown

	getline(ifs, line);		// no. k-mers

	for (int i = 0; !ifs.eof(); ++i)
	{
		getline(ifs, line);
//		cout << i << ":" << line << endl;
//		fflush(stdout);
		parts = split(line, ',');

		if (parts.size() <= 2)
			continue;

//		for (const auto& p : parts)
		for(size_t j = 2; j < parts.size(); ++j)
		{
			const auto p = parts[j];
			elem = split(p, ':');
			if (elem.size() == 2)
			{
				int id = stoi(elem[0]) - 1;			// In kmer-db output indices are 1-based
				double val = stod(elem[1]);

				if (val >= filter_thr)
					filter_vec.emplace_back(encode_pair_id_mm(i, id));
			}
		}
	}

	if(params.verbosity_level > 1)
		cerr << "Filter size: " << filter_vec.size() << endl;

	return true;
}

// ****************************************************************************
bool CLZMatcher::init_data_storage(const vector<string>& input_file_names)
{
	return data_storage.prepare_from_many_files(input_file_names);
}

// ****************************************************************************
bool CLZMatcher::init_data_storage(const string& input_file_name)
{
	return data_storage.prepare_from_multi_fasta(input_file_name);
}

// ****************************************************************************
bool CLZMatcher::prefetch_input_files()
{
	if (!params.buffer_input_data)
		return true;

	if (params.verbosity_level > 1)
	{
		cerr << "Prefetching input files\n";
	}

	mutex mtx;
	atomic<int> fid = 0;

	vector<future<void>> v_fut;
	v_fut.reserve(params.no_threads);

	for (int i = 0; i < params.no_threads; ++i)
		v_fut.push_back(async(std::launch::async, [&] {
		CSharedWorker s_worker(params, data_storage);

		while (true)
		{
			int cid = fid.fetch_add(1);
			if (cid >= input_file_desc.size())
				break;

			if (!s_worker.load_data(input_file_desc[cid]))
			{
				lock_guard<mutex> lck(mtx);
				cout << "Cannot load: " << input_file_desc[cid].file_name << endl;
				exit(0);
			}
		}
		}));

	for (auto& f : v_fut)
		f.get();

	return true;
}

// ****************************************************************************
void CLZMatcher::show_timinigs_info()
{
	if (params.verbosity_level <= 1)
		return;

	cerr << "Timings\n";

	for (size_t i = 1; i < times.size(); ++i)
	{
		cerr << times[i].second << " : " << duration<double>(times[i].first - times[i - 1].first).count() << "s\n";
	}

}

// ****************************************************************************
bool CLZMatcher::prepare_worker_base(CSharedWorker* wb, uint32_t id)
{
	wb->clear_ref();

	if (!wb->load_reference(input_file_desc[id]))
	{
		cerr << "Cannot read: " << input_file_desc[id].file_name << endl;
		return false;
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

	return true;
}


// EOF
