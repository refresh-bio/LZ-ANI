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
	if (params.verbosity_level <= 1)
		return;

	cerr << "Timings\n";

	for (size_t i = 1; i < times.size(); ++i)
		cerr << times[i].second << " : " << duration<double>(times[i].first - times[i - 1].first).count() << "s\n";
}

// ****************************************************************************
void CLZMatcher::run_all2all()
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

	for (uint32_t i = 0; i < params.no_threads; ++i)
	{
		thr_workers.emplace_back([&, i] {
			CParser parser(params);
			uint32_t thread_id = i;
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

	times.emplace_back(high_resolution_clock::now(), "LZ matching");

	cerr << "Storing results" << endl;
	store_results();

	show_timinigs_info();
}

// ****************************************************************************
bool CLZMatcher::store_results()
{
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

		ofs << "id1\tid2";

		if (params.store_condensed)
		{
			if (params.store_total_ani)
				ofs << "\ttotal_ani";
			if (params.store_global_ani)
				ofs << "\tglobal_ani1\tglobal_ani2";
			if (params.store_local_ani)
				ofs << "\tlocal_ani1\tlocal_ani2";
			if (params.store_coverage)
				ofs << "\tcov1\tcov2";
			if (params.store_regions)
				ofs << "\tno_reg1\tno_reg2";
		}
		else
		{
			if (params.store_total_ani)
				ofs << "\ttotal_ani";
			if (params.store_global_ani)
				ofs << "\tglobal_ani";
			if (params.store_local_ani)
				ofs << "\tlocal_ani";
			if (params.store_coverage)
				ofs << "\tcov";
			if (params.store_regions)
				ofs << "\tno_reg";
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

			size_t needed = results[my_id].size() * (8 * 18 + 8);
			if (tmp.size() < needed)
				tmp.resize(needed);
			char* ptr = tmp.data();

			id_results_t x(my_id, {});

			for (auto q = lower_bound(results[my_id].begin(), results[my_id].end(), x); q != results[my_id].end(); ++q)
			{
				if (my_id >= q->id)
					continue;

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

					/*						str += to_string(my_id) + " " + to_string(q->id) + " " +
												to_string(p->results.sym_in_matches) + " " + to_string(p->results.sym_in_literals) + " " + to_string(p->results.no_components) +
												" " +
												to_string(q->results.sym_in_matches) + " " + to_string(q->results.sym_in_literals) + " " + to_string(q->results.no_components) +
												"\n";*/
				}
				else if (params.output_type == output_type_t::split_files)
				{
					vector<CSeqReservoir::item_t>::iterator item[2] = {seq_reservoir.get_sequence(my_id), seq_reservoir.get_sequence(my_id)};

					string names[2] = { sequence_names[my_id] , sequence_names[q->id] };
					uint32_t ids[2] = { my_id, q->id };
					double len[2] = { item[0]->len - (item[0]->no_parts - 1) * (uint32_t) params.close_dist, item[1]->len - (item[1]->no_parts - 1) * (uint32_t) params.close_dist };
					double si_mat[2] = { p->results.sym_in_matches, q->results.sym_in_matches };
					double si_lit[2] = { p->results.sym_in_literals, q->results.sym_in_literals };
					int no_reg[2] = { p->results.no_components, q->results.no_components };

					double total_ani = (si_mat[0] + si_mat[1]) / (len[0] + len[1]);
					double global_ani[2] = { si_mat[0] / len[0], si_mat[1] / len[1] };
					double local_ani[2] = {
						si_mat[0] + si_lit[0] != 0 ? si_mat[0] / (si_mat[0] + si_lit[0]) : 0,
						si_mat[1] + si_lit[1] != 0 ? si_mat[1] / (si_mat[1] + si_lit[1]) : 0 };
					double cov[2] = { (si_mat[0] + si_lit[0]) / len[0], (si_mat[1] + si_lit[1]) / len[1] };

					if (params.store_condensed)
					{
						if (params.store_full_seq_ids)
						{
							strcpy(ptr, names[0].c_str());
							ptr += names[0].length();
							*ptr++ = '\t';
							strcpy(ptr, names[1].c_str());
							ptr += names[1].length();
							*ptr++ = '\t';
						}
						else
						{
							ptr += num2str(ids[0], ptr);							*ptr++ = '\t';
							ptr += num2str(ids[1], ptr);							*ptr++ = '\t';
						}

						if (params.store_total_ani)
						{
							ptr += num2str(total_ani, ptr);		*ptr++ = '\t';
						}
						if (params.store_global_ani)
						{
							ptr += num2str(global_ani[0], ptr);		*ptr++ = '\t';
							ptr += num2str(global_ani[1], ptr);		*ptr++ = '\t';
						}
						if (params.store_local_ani)
						{
							ptr += num2str(local_ani[0], ptr);		*ptr++ = '\t';
							ptr += num2str(local_ani[1], ptr);		*ptr++ = '\t';
						}
						if (params.store_coverage)
						{
							ptr += num2str(cov[0], ptr);		*ptr++ = '\t';
							ptr += num2str(cov[1], ptr);		*ptr++ = '\t';
						}
						if (params.store_regions)
						{
							ptr += num2str(no_reg[0], ptr);		*ptr++ = '\t';
							ptr += num2str(no_reg[1], ptr);		*ptr++ = '\t';
						}
						*ptr++ = '\n';
					}
					else
					{
						for (int i = 0; i < 2; ++i)
						{
							if (params.store_full_seq_ids)
							{
								strcpy(ptr, names[i].c_str());
								ptr += names[i].length();
								*ptr++ = '\t';
								strcpy(ptr, names[!i].c_str());
								ptr += names[!i].length();
								*ptr++ = '\t';
							}
							else
							{
								ptr += num2str(ids[i], ptr);							*ptr++ = '\t';
								ptr += num2str(ids[!i], ptr);							*ptr++ = '\t';
							}

							if (params.store_total_ani)
							{
								ptr += num2str(total_ani, ptr);		*ptr++ = '\t';
							}
							if (params.store_global_ani)
							{
								ptr += num2str(global_ani[i], ptr);		*ptr++ = '\t';
							}
							if (params.store_local_ani)
							{
								ptr += num2str(local_ani[i], ptr);		*ptr++ = '\t';
							}
							if (params.store_coverage)
							{
								ptr += num2str(cov[i], ptr);		*ptr++ = '\t';
							}
							if (params.store_regions)
							{
								ptr += num2str(no_reg[i], ptr);		*ptr++ = '\t';
							}
							*ptr++ = '\n';
						}
					}



				}
/*				else if (params.output_mode == output_mode_t::total_ani)
				{
					double total_ani = (q->results.sym_in_matches + p->results.sym_in_matches) /
						(input_file_desc[my_id].seq_size - input_file_desc[my_id].n_parts * params.close_dist + input_file_desc[q->id].seq_size - input_file_desc[q->id].n_parts * params.close_dist);

					str += to_string(my_id) + " " + to_string(q->id) + " " + to_string(total_ani) + "\n";
				}*/
			}

//			if (params.output_type == output_type_t::single_file)
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
