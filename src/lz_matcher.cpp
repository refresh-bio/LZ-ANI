#include <iostream>
#include <fstream>
#include "lz_matcher.h"
#include "parser.h"
#include "parallel-queues.h"
#include "conversion.h"


char NumericConversions::digits[];
NumericConversions::_si NumericConversions::_init;
uint64_t NumericConversions::powers10[];


// ****************************************************************************
bool CLZMatcher::load_sequences()
{
	if(params.multisample_fasta)
		return seq_reservoir.load_multifasta(params.input_file_names);
	else
		return seq_reservoir.load_fasta(params.input_file_names, params.close_dist);
}

// ****************************************************************************
bool CLZMatcher::load_filter()
{
	if (params.filter_file_name.empty())
		return true;

	return filter.load_filter(params.filter_file_name, params.filter_thr, params.no_threads);
}

// ****************************************************************************
bool CLZMatcher::compare_sequences()
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
void CLZMatcher::reorder_sequences()
{
	auto reordering_map = seq_reservoir.reorder_items();
	filter.reorder_items(reordering_map, params.no_threads);

	cerr << "Reordered" << endl;
}

// ****************************************************************************
void CLZMatcher::show_timinigs_info()
{
//	if (params.verbosity_level <= 1)
//		return;

	cerr << "Timings\n";

	for (size_t i = 1; i < times.size(); ++i)
		cerr << times[i].second << " : " << duration<double>(times[i].first - times[i - 1].first).count() << "s\n";

	cerr << "Total time: " << duration<double>(times.back().first - times.front().first).count() << "s\n";
}

// ****************************************************************************
void CLZMatcher::do_matching()
{
	cerr << "All2all sparse" << endl;

	results.clear();
	results.resize(seq_reservoir.size());

	// Start worker threads
	vector<thread> thr_workers;
	thr_workers.reserve(params.no_threads);

	atomic<uint64_t> global_task_no = 0;
	atomic<uint64_t> global_no_pairs = 0;

	for (uint32_t i = 0; i < params.no_threads; ++i)
	{
		thr_workers.emplace_back([&] {
			CParser parser(params);
			uint64_t local_task_no = 0;

			vec_id_results_t res_row;

			while (true)
			{
				local_task_no = global_task_no.fetch_add(1);

				if (local_task_no >= seq_reservoir.size())
					break;

				res_row.clear();

				// Prepare reference
				auto sr_iter = seq_reservoir.get_sequence(local_task_no);
				parser.prepare_reference(seq_view(sr_iter->data, sr_iter->len), sr_iter->no_parts);

				uint64_t to_add = filter.is_empty() ? local_task_no : (uint64_t)filter.get_row(local_task_no).size();
				auto to_print = global_no_pairs.fetch_add(to_add);

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

//				res_row.emplace_back(local_task_no, results_t(1, 0, 1));
				res_row.shrink_to_fit();
				sort(res_row.begin(), res_row.end());

				results[local_task_no] = move(res_row);
			}
			});
	}

	for (auto& thr : thr_workers)
		thr.join();

//	seq_reservoir.release();

}

// ****************************************************************************
bool CLZMatcher::store_results()
{
	cerr << "Storing results" << endl;

	string fn = params.output_file_name;

	if (params.output_type == output_type_t::split_files)
		fn += ".ids.tsv";

	ofstream ofs(fn);

	if (!ofs.is_open())
	{
		cerr << "Cannot open output file: " << fn << endl;
		return false;
	}

	const size_t io_buffer_size = 32 << 20;
	char* io_buffer = new char[io_buffer_size];
	ofs.rdbuf()->pubsetbuf(io_buffer, io_buffer_size);

	if (params.output_type == output_type_t::single_file)
	{
		// Store params
		ofs << params.str();

		// Store file desc
		ofs << "[no_input_sequences]" << endl;
		ofs << seq_reservoir.size() << endl;

		ofs << "[input_sequences]" << endl;
		for (size_t i = 0; i < seq_reservoir.size(); ++i)
		{
			auto x = seq_reservoir.get_sequence(i);
			ofs << x->name << " " << x->len - (x->no_parts - 1) * params.close_dist << " " << x->no_parts << endl;
		}

		ofs << "[lz_similarities]" << endl;
	}
	else
	{
		ofs << "id\tseq_len\tno_parts\n";
		for (size_t i = 0; i < seq_reservoir.size(); ++i)
		{
			auto x = seq_reservoir.get_sequence(i);
			ofs << x->name << "\t" << x->len - (x->no_parts - 1) * params.close_dist << "\t" << x->no_parts << endl;
		}

		ofs.close();

		fn = params.output_file_name + ".data.tsv";
		ofs.open(fn);

		if (!ofs.is_open())
		{
			cerr << "Cannot open output file: " << fn << endl;
			return false;
		}

		bool first = true;
		for (const auto x : params.output_components)
		{
			if (!first)
				ofs << "\t";
			else
				first = false;
			ofs << params.comp_id_name[x];
		}

		ofs << endl;
	}

	refresh::parallel_priority_queue<string> par_queue(params.no_threads * 64, params.no_threads);
	atomic<uint64_t> id_global = 0;

	auto sequence_names = seq_reservoir.get_sequence_names();

	vector<thread> thr_workers;
	thr_workers.reserve(params.no_threads);

	for (uint32_t i = 0; i < params.no_threads; ++i)
		thr_workers.emplace_back([&] {
		string str;
		string tmp;

		while (true)
		{
			uint64_t my_id = id_global.fetch_add(1);
			if (my_id >= results.size())
			{
				par_queue.mark_completed();
				break;
			}

			str.clear();
			const size_t max_line_len = params.output_components.size() * 100;
			str.resize(2 * max_line_len);

/*			size_t needed = results[my_id].size() * (8 * 18 + 8);
			if (tmp.size() < needed)
				tmp.resize(needed);*/
			char* ptr = tmp.data();

			id_results_t x(my_id, {});

			for (auto q = lower_bound(results[my_id].begin(), results[my_id].end(), x); q != results[my_id].end(); ++q)
			{
				if (my_id >= q->id)
					continue;

				if (tmp.size() - (ptr - tmp.data()) < 2 * max_line_len)
				{
					tmp.resize((size_t) (tmp.size() * 1.3) + 2 * max_line_len);
					ptr = tmp.data();
				}

				auto p = lower_bound(results[q->id].begin(), results[q->id].end(), x);
				assert(p != results[q->id].end() && p->id == my_id);

				if (params.output_type == output_type_t::single_file)
				{
					ptr += num2str(my_id, ptr);							*ptr++ = ' ';
					ptr += num2str(q->id, ptr);							*ptr++ = ' ';
					ptr += num2str(p->results.sym_in_matches, ptr);		*ptr++ = ' ';
					ptr += num2str(p->results.sym_in_literals, ptr);	*ptr++ = ' ';
					ptr += num2str(p->results.no_components, ptr);		*ptr++ = ' ';
					ptr += num2str(q->results.sym_in_matches, ptr);		*ptr++ = ' ';
					ptr += num2str(q->results.sym_in_literals, ptr);	*ptr++ = ' ';
					ptr += num2str(q->results.no_components, ptr);		*ptr++ = '\n';
				}
				else if (params.output_type == output_type_t::split_files)
				{
					vector<CSeqReservoir::item_t>::iterator item[2] = {seq_reservoir.get_sequence(my_id), seq_reservoir.get_sequence(q->id)};

					string names[2] = { sequence_names[my_id] , sequence_names[q->id] };
					uint32_t ids[2] = { (uint32_t) my_id, q->id };
					uint32_t len[2] = { item[0]->len - (item[0]->no_parts - 1) * (uint32_t) params.close_dist, item[1]->len - (item[1]->no_parts - 1) * (uint32_t) params.close_dist };
					int32_t si_mat[2] = { p->results.sym_in_matches, q->results.sym_in_matches };
					int32_t si_lit[2] = { p->results.sym_in_literals, q->results.sym_in_literals };
					int no_reg[2] = { p->results.no_components, q->results.no_components };

					double total_ani = (double) (si_mat[0] + si_mat[1]) / (len[0] + len[1]);
					double global_ani[2] = { (double) si_mat[0] / len[0], (double) si_mat[1] / len[1] };
					double local_ani[2] = {
						si_mat[0] + si_lit[0] != 0 ? (double) si_mat[0] / (si_mat[0] + si_lit[0]) : 0,
						si_mat[1] + si_lit[1] != 0 ? (double) si_mat[1] / (si_mat[1] + si_lit[1]) : 0 };
					double cov[2] = { (double)(si_mat[0] + si_lit[0]) / len[0], (double)(si_mat[1] + si_lit[1]) / len[1] };

					for (int i = 0; i < 2; ++i)
					{
						for(auto oc : params.output_components)
							if (oc == output_component_t::seq_idx1)
							{
								ptr += num2str(ids[i], ptr);						*ptr++ = '\t';
							}
							else if (oc == output_component_t::seq_idx2)
							{
								ptr += num2str(ids[!i], ptr);						*ptr++ = '\t';
							}
							else if (oc == output_component_t::seq_id1)
							{
								strcpy(ptr, names[i].c_str());
								ptr += names[i].length();
								*ptr++ = '\t';
							}
							else if (oc == output_component_t::seq_id2)
							{
								strcpy(ptr, names[!i].c_str());
								ptr += names[!i].length();
								*ptr++ = '\t';
							}
							else if (oc == output_component_t::cov)
							{
								ptr += num2str(cov[i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::global_ani)
							{
								ptr += num2str(global_ani[i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::len1)
							{
								ptr += num2str(len[i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::len2)
							{
								ptr += num2str(len[!i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::len_ratio)
							{
								if (len[0] && len[1])
									ptr += num2str((double)len[i] / len[!i], ptr);
								else
									ptr += num2str(0, ptr);
								*ptr++ = '\t';
							}
							else if (oc == output_component_t::local_ani)
							{
								ptr += num2str(local_ani[i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::no_reg)
							{
								ptr += num2str(no_reg[i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::sim)
							{
								ptr += num2str(local_ani[i] * cov[i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::sym_lit)
							{
								ptr += num2str(si_lit[i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::sym_mat)
							{
								ptr += num2str(si_mat[i], ptr);		*ptr++ = '\t';
							}
							else if (oc == output_component_t::total_ani)
							{
								ptr += num2str(total_ani, ptr);		*ptr++ = '\t';
							}

						if (!params.output_components.empty())
							--ptr;			// overwrite last \t

						*ptr++ = '\n';
					}
				}
			}

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

// ****************************************************************************
bool CLZMatcher::run_all2all()
{
	times.clear();
	times.emplace_back(high_resolution_clock::now(), "");

	if (!load_sequences())
		return false;

	times.emplace_back(high_resolution_clock::now(), "Loading sequences");

	if (!load_filter())
		return false;

	times.emplace_back(high_resolution_clock::now(), "Loading filter");

	if (!compare_sequences())
		return false;

	times.emplace_back(high_resolution_clock::now(), "Comparing sequence and filter compatibility");

	reorder_sequences();

	times.emplace_back(high_resolution_clock::now(), "Reordering sequences");

	do_matching();

	times.emplace_back(high_resolution_clock::now(), "LZ matching");

	store_results();

	times.emplace_back(high_resolution_clock::now(), "Storing results");

	show_timinigs_info();

	return true;
}

// EOF
