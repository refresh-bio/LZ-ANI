#include <iostream>
#include <fstream>
#include "lz_matcher2.h"
#include "parser.h"
#include "parallel-queues.h"
#include "conversion.h"

// ****************************************************************************
bool CLZMatcher2::load_sequences()
{
	if(params.multisample_fasta)
		return seq_reservoir.load_multifasta(params.input_file_names);
	else
		return seq_reservoir.load_fasta(params.input_file_names);
}

// ****************************************************************************
bool CLZMatcher2::load_filter()
{
	return filter.load_filter(params.filter_file_name, params.filter_thr, params.no_threads);
}


// ****************************************************************************
bool CLZMatcher2::compare_sequences()
{
	if (filter.is_empty())
		return true;

	auto seq_sn = seq_reservoir.get_sequence_names();
	auto flt_sn = filter.get_sequence_names();

	if (seq_sn.size() != flt_sn.size() || seq_sn != flt_sn)
	{
		cerr << "Input sequences and filter sequences are different!" << endl;
		return false;
	}

	filter.clear_sequence_names();

	return true;
}

// ****************************************************************************
void CLZMatcher2::reorder_sequences()
{
	auto reordering_map = seq_reservoir.reorder_items();
	filter.reorder_items(reordering_map, params.no_threads);

	cerr << "Reordered" << endl;
}

// ****************************************************************************
void CLZMatcher2::show_timinigs_info()
{
	if (params.verbosity_level <= 1)
		return;

	cerr << "Timings\n";

	for (size_t i = 1; i < times.size(); ++i)
		cerr << times[i].second << " : " << duration<double>(times[i].first - times[i - 1].first).count() << "s\n";
}

// ****************************************************************************
void CLZMatcher2::run_all2all()
{
	cerr << "All2all sparse" << endl;

	times.clear();
	times.emplace_back(high_resolution_clock::now(), "");

	results.clear();
	results.resize(seq_reservoir.size());

	// Start worker threads
	vector<thread> thr_workers;
	thr_workers.reserve(params.no_threads);

	atomic<uint64_t> global_task_no = 0;
	atomic<uint64_t> global_no_pairs = 0;

	for (int i = 0; i < params.no_threads; ++i)
	{
		thr_workers.emplace_back([&, i] {
			CParser parser(params);
			int thread_id = i;
			uint64_t local_task_no = 0;

			VecIdResults res_row;

			while (true)
			{
				local_task_no = global_task_no.fetch_add(1);

				if (local_task_no >= seq_reservoir.size())
					break;

				res_row.clear();

				// Prepare reference
				auto sr_iter = seq_reservoir.get_sequence(local_task_no);
				parser.prepare_reference(seq_view(sr_iter->data, sr_iter->len), sr_iter->no_parts);

				auto to_print = global_no_pairs.fetch_add((uint64_t)filter.get_row(local_task_no).size());

				cerr << to_string(local_task_no) + " : " + to_string(to_print) + "    \r";

				if (filter.is_empty())
				{
					for (uint64_t id = 0; id < local_task_no; ++id)
					{
						auto sr_iter = seq_reservoir.get_sequence(id);
						parser.prepare_data(seq_view(sr_iter->data, sr_iter->len), sr_iter->no_parts);

						parser.parse();

						auto results = parser.calc_stats();

						res_row.emplace_back(id, results);
					}
				}
				else
				{
					for (auto& id : filter.get_row(local_task_no))
					{
						auto sr_iter = seq_reservoir.get_sequence(id);
						parser.prepare_data(seq_view(sr_iter->data, sr_iter->len), sr_iter->no_parts);

						parser.parse();

						auto results = parser.calc_stats();

						res_row.emplace_back(id, results);
					}

					filter.clear_row(local_task_no);
				}

//				res_row.emplace_back(local_task_no, CResults(1, 0, 1));
				res_row.shrink_to_fit();
				sort(res_row.begin(), res_row.end());

				results[local_task_no] = move(res_row);
			}
			});
	}

	for (auto& thr : thr_workers)
		thr.join();

//	seq_reservoir.release();

	times.emplace_back(high_resolution_clock::now(), "LZ matching");

	cerr << "Storing results" << endl;
	store_results();

	show_timinigs_info();
}

// ****************************************************************************
bool CLZMatcher2::store_results()
{
	ofstream ofs(params.output_file_name);

	if (!ofs.is_open())
	{
		cerr << "Cannot open output file: " << params.output_file_name << endl;
		return false;
	}

	const size_t io_buffer_size = 32 << 20;
	char* io_buffer = new char[io_buffer_size];
	ofs.rdbuf()->pubsetbuf(io_buffer, io_buffer_size);

	// Store params
	ofs << params.str();

	// Store file desc
	ofs << "[no_input_sequences]" << endl;
	ofs << seq_reservoir.size() << endl;

	ofs << "[input_sequences]" << endl;
	for (size_t i = 0; i < seq_reservoir.size(); ++i)
	{
		auto x = seq_reservoir.get_sequence(i);
//		ofs << x->name << " " << x->len - x->no_parts * params.close_dist << " " << x->no_parts << endl;
		ofs << x->name << " " << x->len << " " << x->no_parts << endl;
	}

	// Store mapping results in sparse form
	uint32_t q1, q2;

	ofs << "[lz_similarities]" << endl;

	refresh::parallel_priority_queue<string> par_queue(params.no_threads * 64, params.no_threads);
	atomic<uint64_t> id_global = 0;

	vector<thread> thr_workers;
	thr_workers.reserve(params.no_threads);

	for (int i = 0; i < params.no_threads; ++i)
		thr_workers.emplace_back([&] {
		string str;
		string tmp;

		while (true)
		{
			int my_id = id_global.fetch_add(1);
			if (my_id >= results.size())
			{
				par_queue.mark_completed();
				break;
			}

			str.clear();

			size_t needed = results[my_id].size() * (8 * 18 + 8);
			if (tmp.size() < needed)
				tmp.resize(needed);
			char* ptr = tmp.data();

			CIdResults x(my_id, {});

			for (auto q = lower_bound(results[my_id].begin(), results[my_id].end(), x); q != results[my_id].end(); ++q)
			{
				if (my_id >= q->id)
					continue;

				auto p = lower_bound(results[q->id].begin(), results[q->id].end(), x);
				assert(p != results2[q->id].end() && p->id == my_id);

				if (params.output_mode == output_mode_t::sym_in_matches_literals)
				{
					ptr += num2str(my_id, ptr);							*ptr++ = ' ';
					ptr += num2str(q->id, ptr);							*ptr++ = ' ';
					ptr += num2str(p->results.sym_in_matches, ptr);		*ptr++ = ' ';
					ptr += num2str(p->results.sym_in_literals, ptr);	*ptr++ = ' ';
					ptr += num2str(p->results.no_components, ptr);		*ptr++ = ' ';
					ptr += num2str(q->results.sym_in_matches, ptr);		*ptr++ = ' ';
					ptr += num2str(q->results.sym_in_literals, ptr);	*ptr++ = ' ';
					ptr += num2str(q->results.no_components, ptr);		*ptr++ = '\n';

					/*						str += to_string(my_id) + " " + to_string(q->id) + " " +
												to_string(p->results.sym_in_matches) + " " + to_string(p->results.sym_in_literals) + " " + to_string(p->results.no_components) +
												" " +
												to_string(q->results.sym_in_matches) + " " + to_string(q->results.sym_in_literals) + " " + to_string(q->results.no_components) +
												"\n";*/
				}
/*				else if (params.output_mode == output_mode_t::total_ani)
				{
					double total_ani = (q->results.sym_in_matches + p->results.sym_in_matches) /
						(input_file_desc[my_id].seq_size - input_file_desc[my_id].n_parts * params.close_dist + input_file_desc[q->id].seq_size - input_file_desc[q->id].n_parts * params.close_dist);

					str += to_string(my_id) + " " + to_string(q->id) + " " + to_string(total_ani) + "\n";
				}*/
			}

			if (params.output_mode == output_mode_t::sym_in_matches_literals)
				str.assign(tmp.data(), ptr);

			par_queue.push(my_id, move(str));
		}
			});

	string to_print;

	while (par_queue.pop(to_print))
		ofs.write(to_print.data(), to_print.size());

	for (auto& t : thr_workers)
		t.join();

	ofs.close();
	delete[] io_buffer;

	return true;
}


// EOF
