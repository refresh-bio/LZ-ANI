// *******************************************************************************************
// This file is a part of LZ-ANI software distributed under GNU GPL 3 license.
// The homepage of the LZ-ANI project is https://github.com/refresh-bio/LZ-ANI
//
// Copyright(C) 2024-2024, S.Deorowicz, A.Gudys
//
// Version: 1.1.0
// Date   : 2024-09-05
// *******************************************************************************************

#include <iostream>
#include <fstream>
#include "lz_matcher.h"
#include "parser.h"
#include "../libs/refresh/parallel_queues/lib/parallel-queues.h"
#include "../libs/refresh/conversions/lib/numeric_conversions.h"

// ****************************************************************************
bool CLZMatcher::load_sequences()
{
	if (params.verbosity_level >= 1)
		cerr << "Loading sequences\n";

	if(params.multisample_fasta)
		return seq_reservoir.load_multifasta(params.input_file_names, params.verbosity_level);
	else
		return seq_reservoir.load_fasta(params.input_file_names, params.max_dist_in_ref, params.verbosity_level);
}

// ****************************************************************************
bool CLZMatcher::load_filter()
{
	if (params.verbosity_level >= 1)

	if (params.filter_file_name.empty())
		return true;

	return filter.load_filter(params.filter_file_name, params.filter_thr, params.no_threads, params.verbosity_level);
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
	auto reordering_map = seq_reservoir.reorder_items(params.verbosity_level);
	filter.reorder_items(reordering_map, params.no_threads, params.verbosity_level);

	if (params.verbosity_level > 1)
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

	cerr << "Total time: " << duration<double>(times.back().first - times.front().first).count() << "s\n";
}

// ****************************************************************************
void CLZMatcher::store_alignment(uint64_t idx1, uint64_t idx2, const vector<region_t>& parsing)
{
	string s1 = seq_reservoir.get_sequence_name(idx1);
	string s2 = seq_reservoir.get_sequence_name(idx2);

	int s1_len = s1.length();
	int s2_len = s2.length();

	int seq1_len = seq_reservoir.get_sequence_length(idx1);
	int seq2_len = seq_reservoir.get_sequence_length(idx2);

	int rc_correction = 2 * seq1_len + 2 * params.max_dist_in_ref + 1;

	// Do not output alignment that should be filtered out. Here only a part of allowed filters can be applied
	
	int32_t si_mat = 0;
	int32_t si_lit = 0;

	for (const auto& region : parsing)
	{
		si_mat += region.num_matches;
		si_lit += region.num_mismatches;
	}

	double global_ani = (double)si_mat / seq2_len;
	double local_ani = si_mat + si_lit != 0 ? (double)si_mat / (si_mat + si_lit) : 0;
	double qcov = (double)(si_mat + si_lit) / seq2_len;

	if (params.output_filter_mask != 0)
	{
		if (global_ani < params.output_filter_vals[(uint32_t)output_component_t::gani])
			return;
		if (local_ani < params.output_filter_vals[(uint32_t)output_component_t::ani])
			return;
		if (qcov < params.output_filter_vals[(uint32_t)output_component_t::qcov])
			return;
	}

	lock_guard<mutex> lck(mtx_alignment);

	char *p = alignment_buffer;

	for (const auto& region: parsing)
	{
		strcpy(p, s2.c_str());		p += s2_len;	*p++ = '\t';
		strcpy(p, s1.c_str());		p += s1_len;	*p++ = '\t';
		p += refresh::real_to_pchar(100.0 * region.num_matches / region.length(), p, 6, '\t');
		p += refresh::int_to_pchar(region.length(), p, '\t');
		p += refresh::int_to_pchar(1 + region.seq_start, p, '\t');
		p += refresh::int_to_pchar(1 + region.seq_end - 1, p, '\t');

		if (region.ref_start < seq1_len)
		{
			p += refresh::int_to_pchar(1 + region.ref_start, p, '\t');
			p += refresh::int_to_pchar(1 + region.ref_end - 1, p, '\t');
		}
		else
		{
			p += refresh::int_to_pchar(rc_correction - (1 + region.ref_start), p, '\t');
			p += refresh::int_to_pchar(rc_correction - (1 + region.ref_end - 1), p, '\t');
		}
		p += refresh::int_to_pchar(region.num_matches, p, '\t');
		p += refresh::int_to_pchar(region.num_mismatches, p, '\n');

		fwrite(alignment_buffer, 1, p - alignment_buffer, f_alignment);
		p = alignment_buffer;
	}
}

// ****************************************************************************
void CLZMatcher::do_matching()
{
	if (params.verbosity_level >= 1)
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

				if (params.verbosity_level >= 2)
					cerr << to_string(local_task_no) + " : " + to_string(to_print) + "    \r";

				if (filter.is_empty())
				{
//					for (uint64_t id = 0; id < local_task_no; ++id)
					for (uint64_t id = 0; id < seq_reservoir.size(); ++id)
					{
						if (id == local_task_no)
							continue;

						auto sr_iter = seq_reservoir.get_sequence(id);
						parser.prepare_data(seq_view(sr_iter->data, sr_iter->len), sr_iter->no_parts);

						parser.parse();

						if(f_alignment)
							store_alignment(local_task_no, id, parser.get_parsing());

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

						if (f_alignment)
							store_alignment(local_task_no, id, parser.get_parsing());

						auto results = parser.calc_stats();

						res_row.emplace_back(id, results);
					}

					filter.clear_row(local_task_no);
				}

				res_row.shrink_to_fit();
				sort(res_row.begin(), res_row.end());

				results[local_task_no] = move(res_row);
			}
			});
	}

	for (auto& thr : thr_workers)
		thr.join();

	if (params.verbosity_level >= 2)
		cerr << endl;
}

// ****************************************************************************
bool CLZMatcher::store_results()
{
	if (params.verbosity_level >= 1)
		cerr << "Storing results" << endl;

	string fn, fn_anis;
	const string sep = "\t";
	const char sep_char = '\t';
	double val_mult = params.output_in_percent ? 100 : 1;

	if (params.output_type == output_type_t::two_tsv)
	{
		fn_anis = params.output_file_name;
		fn = params.output_ids_file_name;

		if (fn.empty())
		{
			auto p = fn_anis.rfind(".");
			if (p == string::npos)
				fn = fn_anis + ".ids";
			else
				fn = fn_anis.substr(0, p) + ".ids" + fn_anis.substr(p, string::npos);
		}
	}
	else if (params.output_type == output_type_t::single_txt)
	{
		fn = params.output_file_name;
	}

	ofstream ofs(fn);

	if (!ofs.is_open())
	{
		cerr << "Cannot open output file: " << fn << endl;
		return false;
	}

	const size_t io_buffer_size = 32 << 20;
	char* io_buffer = new char[io_buffer_size];
	ofs.rdbuf()->pubsetbuf(io_buffer, io_buffer_size);

	if (params.output_type == output_type_t::single_txt)
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
			ofs << x->name << " " << x->len - (x->no_parts - 1) * params.max_dist_in_ref << " " << x->no_parts << endl;
		}

		ofs << "[lz_similarities]" << endl;
	}
	else
	{
		ofs << "id" + sep + "seq_len" + sep + "no_parts\n";
		for (size_t i = 0; i < seq_reservoir.size(); ++i)
		{
			auto x = seq_reservoir.get_sequence(i);
			ofs << x->name << sep << x->len - (x->no_parts - 1) * params.max_dist_in_ref << sep << x->no_parts << endl;
		}

		ofs.close();

		ofs.open(fn_anis);

		if (!ofs.is_open())
		{
			cerr << "Cannot open output file: " << fn << endl;
			return false;
		}

		bool first = true;
		for (const auto x : params.output_components)
		{
			if (!first)
				ofs << sep;
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
			uint64_t ref_id = id_global.fetch_add(1);
			if (ref_id >= results.size())
			{
				par_queue.mark_completed();
				break;
			}

			tmp.clear();
			const size_t max_line_len = params.output_components.size() * 100;
			tmp.resize(2 * max_line_len);

/*			size_t needed = results[my_id].size() * (8 * 18 + 8);
			if (tmp.size() < needed)
				tmp.resize(needed);*/
			char* ptr = tmp.data();

			id_results_t x((uint32_t)ref_id, {});

			for (auto q = lower_bound(results[ref_id].begin(), results[ref_id].end(), x); q != results[ref_id].end(); ++q)
			{
				if (ref_id >= q->id)
					continue;

				size_t cur_size = ptr - tmp.data();

				if (tmp.size() - cur_size < 2 * max_line_len)
				{
					tmp.resize((size_t) (tmp.size() * 1.3) + 2 * max_line_len);
					ptr = tmp.data() + cur_size;
				}

				auto p = lower_bound(results[q->id].begin(), results[q->id].end(), x);
				assert(p != results[q->id].end() && p->id == ref_id);

				if (params.output_type == output_type_t::single_txt)
				{
					ptr += refresh::int_to_pchar(ref_id, ptr, ' ');
					ptr += refresh::int_to_pchar(q->id, ptr, ' ');
					ptr += refresh::int_to_pchar(p->results.sym_in_matches, ptr, ' ');
					ptr += refresh::int_to_pchar(p->results.sym_in_literals, ptr, ' ');
					ptr += refresh::int_to_pchar(p->results.no_components, ptr, ' ');
					ptr += refresh::int_to_pchar(q->results.sym_in_matches, ptr, ' ');
					ptr += refresh::int_to_pchar(q->results.sym_in_literals, ptr, ' ');
					ptr += refresh::int_to_pchar(q->results.no_components, ptr, '\n');
				}
				else if (params.output_type == output_type_t::two_tsv)
				{
					vector<CSeqReservoir::item_t>::iterator item[2] = {seq_reservoir.get_sequence(ref_id), seq_reservoir.get_sequence(q->id)};

					string names[2] = { sequence_names[ref_id], sequence_names[q->id] };
					uint32_t ids[2] = { (uint32_t)ref_id, q->id };
					uint32_t len[2] = { item[1]->len - (item[1]->no_parts - 1) * (uint32_t)params.max_dist_in_ref, item[0]->len - (item[0]->no_parts - 1) * (uint32_t) params.max_dist_in_ref };
					int32_t si_mat[2] = { q->results.sym_in_matches, p->results.sym_in_matches};
					int32_t si_lit[2] = { q->results.sym_in_literals, p->results.sym_in_literals};
					int no_reg[2] = { q->results.no_components, p->results.no_components};

					double total_ani = (double) (si_mat[0] + si_mat[1]) / (len[0] + len[1]);
					double global_ani[2] = { (double) si_mat[0] / len[0], (double) si_mat[1] / len[1] };
					double local_ani[2] = {
						si_mat[0] + si_lit[0] != 0 ? (double) si_mat[0] / (si_mat[0] + si_lit[0]) : 0,
						si_mat[1] + si_lit[1] != 0 ? (double) si_mat[1] / (si_mat[1] + si_lit[1]) : 0 };
					double cov[2] = { (double)(si_mat[0] + si_lit[0]) / len[0], (double)(si_mat[1] + si_lit[1]) / len[1] };

					for (int i = 0; i < 2; ++i)
					{
						if (params.output_filter_mask != 0)
						{
							if (global_ani[i] < params.output_filter_vals[(uint32_t)output_component_t::gani])
								continue;
							if (local_ani[i] < params.output_filter_vals[(uint32_t)output_component_t::ani])
								continue;
							if (total_ani < params.output_filter_vals[(uint32_t)output_component_t::tani])
								continue;
							if (cov[i] < params.output_filter_vals[(uint32_t)output_component_t::qcov])
								continue;
							if (cov[!i] < params.output_filter_vals[(uint32_t)output_component_t::rcov])
								continue;
						}

						for(auto oc : params.output_components)
							if (oc == output_component_t::ridx)
							{
								ptr += refresh::int_to_pchar(ids[i], ptr, sep_char);
							}
							else if (oc == output_component_t::qidx)
							{
								ptr += refresh::int_to_pchar(ids[!i], ptr, sep_char);
							}
							else if (oc == output_component_t::reference)
							{
								strcpy(ptr, names[i].c_str());
								ptr += names[i].length();
								*ptr++ = sep_char;
							}
							else if (oc == output_component_t::query)
							{
								strcpy(ptr, names[!i].c_str());
								ptr += names[!i].length();
								*ptr++ = sep_char;
							}
							else if (oc == output_component_t::qcov)
							{
								ptr += refresh::real_to_pchar(val_mult * cov[i], ptr, 6, sep_char);
							}
							else if (oc == output_component_t::rcov)
							{
								ptr += refresh::real_to_pchar(val_mult * cov[!i], ptr, 6, sep_char);
							}
							else if (oc == output_component_t::gani)
							{
								ptr += refresh::real_to_pchar(val_mult * global_ani[i], ptr, 6, sep_char);
							}
							else if (oc == output_component_t::rlen)
							{
								ptr += refresh::int_to_pchar(len[!i], ptr, sep_char);
							}
							else if (oc == output_component_t::qlen)
							{
								ptr += refresh::int_to_pchar(len[i], ptr, sep_char);
							}
							else if (oc == output_component_t::len_ratio)
							{
								if (len[0] && len[1])
								{
									double len_ratio = 0;
									if (len[i] < len[!i])
										len_ratio = (double) len[i] / len[!i];
									else
										len_ratio = (double) len[!i] / len[i];

									ptr += refresh::real_to_pchar(len_ratio, ptr, 4, sep_char);
								}
								else
								{
									*ptr++ = '0';	*ptr++ = sep_char;
								}
							}
							else if (oc == output_component_t::ani)
							{
								ptr += refresh::real_to_pchar(val_mult * local_ani[i], ptr, 6, sep_char);
							}
							else if (oc == output_component_t::num_alns)
							{
								ptr += refresh::int_to_pchar(no_reg[i], ptr, sep_char);
							}
							else if (oc == output_component_t::nt_mismatch)
							{
								ptr += refresh::int_to_pchar(si_lit[i], ptr, sep_char);
							}
							else if (oc == output_component_t::nt_match)
							{
								ptr += refresh::int_to_pchar(si_mat[i], ptr, sep_char);
							}
							else if (oc == output_component_t::tani)
							{
								ptr += refresh::real_to_pchar(val_mult * total_ani, ptr, 6, sep_char);
							}

						if (!params.output_components.empty())
							--ptr;			// overwrite last sep

						*ptr++ = '\n';
					}
				}
			}

			str.assign(tmp.data(), ptr);

			par_queue.push(ref_id, move(str));
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
