#include <iostream>

#include "filter.h"
#include "utils.h"
#include "../libs/refresh/parallel-queues.h"
#include "../libs/refresh/conversion.h"


// ****************************************************************************
bool CFilter::load_filter(const string& fn, double thr, uint32_t no_threads)
{
	refresh::stream_in_file fil(fn, 16 << 20, 16 << 20);

	if (!fil.is_open())
	{
		cerr << "Cannot open file: " << fn << endl;
		return false;
	}

	refresh::stream_decompression sfil(&fil, 16 << 20);

	string line;
	vector<string> parts;
	vector<string> elem;

	sfil.getline(line);		// genome names
	sequence_names = split(line, ',');
	if (sequence_names.size() <= 2)
	{
		cerr << "Incorrect kmer-db filter file\n";
		return false;
	}

	sequence_names.erase(sequence_names.begin());

	filter.resize(sequence_names.size());

	line.clear();
	line.shrink_to_fit();

	uint64_t no_items = 0;

	cerr << "Start loading filter data (no_threads: " << no_threads << ")" << endl;

	if (no_threads < 4)
	{
		cerr << "Single-threaded ver." << endl;
		for (int i = 0; !sfil.eof(); ++i)
		{
			sfil.getline(line);
			parts = split(line, ',');

			if (parts.size() <= 1)
				continue;

//			uint32_t reo_i = reordering[i];

			//		for (const auto& p : parts)
			for (size_t j = 1; j < parts.size(); ++j)
			{
				const auto p = parts[j];
				elem = split(p, ':');
				if (elem.size() == 2)
				{
					int id = stoi(elem[0]) - 1;			// In kmer-db output indices are 1-based
					double val = stod(elem[1]);

					if (val >= thr)
					{
//						uint32_t reo_id = reordering[id];

//						filter[reo_i].emplace_back(reo_id);
//						filter[reo_id].emplace_back(reo_i);
						filter[i].emplace_back(id);
						filter[id].emplace_back(i);
					}
				}
			}

//			filter[reo_i].shrink_to_fit();
			filter[i].shrink_to_fit();
//			no_items += filter[reo_i].size();
			no_items += filter[i].size();
		}
	}
	else
	{
		cerr << "Multi-threaded ver." << endl;

		refresh::parallel_queue<string> pq_lines(128, 1);
		refresh::parallel_queue<vector<pair<string, string>>> pq_parts(128, 1);
		refresh::parallel_queue<vector<uint32_t>> pq_ids(128, 1);

		thread thr_reader([&]
			{
				while (!sfil.eof())
				{
					string line;
					sfil.getline(line);
					if (line.length() > 2)
						pq_lines.push(move(line));
				}

				pq_lines.mark_completed();
			});

		thread thr_splitter([&]
			{
				string line;
				vector<string> parts;

				while (pq_lines.pop(line))
				{
					parts = split(line, ',');

					vector<pair<string, string>> splitted_parts;

					if (parts.size() == 1)
					{
						pq_parts.push(move(splitted_parts));
						continue;
					}
					else if (parts.size() == 0)
						continue;

					splitted_parts.reserve(parts.size() - 1);

					for (size_t j = 1; j < parts.size(); ++j)
					{
						const auto p = parts[j];
						elem = split(p, ':');
						if (elem.size() == 2)
							splitted_parts.emplace_back(elem[0], elem[1]);
					}

					pq_parts.push(move(splitted_parts));
				}

				pq_parts.mark_completed();
			});

		thread thr_converter([&]
			{
				vector<pair<string, string>> splitted_parts;

				while (pq_parts.pop(splitted_parts))
				{
					vector<uint32_t> ids;
					char* end;

					ids.reserve(splitted_parts.size());

					for (const auto& elem : splitted_parts)
					{
						double val = stod(elem.second);

						if (val >= thr)
							ids.emplace_back(NumericConversions::strtol(elem.first.c_str(), &end) - 1);
					}

					pq_ids.push(move(ids));
				}

				pq_ids.mark_completed();
			});

		vector<uint32_t> ids;

		int row_id = 0;
		vector<uint32_t> no_items_to_extend(filter.size(), 0);

		while (pq_ids.pop(ids))
		{
//			uint32_t reo_i = reordering[row_id];

//			auto& row = filter[reo_i];
			auto& row = filter[row_id];
			row.reserve(ids.size());

			//		for (const auto& p : parts)
			for (auto& elem : ids)
			{
//				uint32_t reo_id = reordering[elem];
//				row.emplace_back(reo_id);
				row.emplace_back(elem);

//				++no_items_to_extend[reo_id];
				++no_items_to_extend[elem];

				//				filter_map[reo_id].emplace_back(reo_i);
			}

			++row_id;

			if(row_id % 1000 == 0)
				cout << row_id << "\r";
		}
		cout << row_id << "\n";

		thr_reader.join();
		thr_splitter.join();
		thr_converter.join();

		vector<uint32_t> first_pass_sizes(filter.size(), 0);

		cerr << "End of loading data" << endl;

		for (size_t i = 0; i < filter.size(); ++i)
		{
			first_pass_sizes[i] = (uint32_t)filter[i].size();
			filter[i].reserve(filter[i].size() + no_items_to_extend[i]);

			no_items += filter[i].size();
		}

		vector<thread> thr_workers;
		thr_workers.reserve(no_threads);

		cerr << "End of resizing filter vectors" << endl;

		for (uint32_t i = 0; i < no_threads; ++i)
			thr_workers.emplace_back([&, i] {
			uint32_t thread_id = i;

			for (size_t j = 0; j < filter.size(); ++j)
			{
				auto& row = filter[j];
				for (uint32_t k = 0; k < first_pass_sizes[j]; ++k)
				{
					auto id = row[k];
					if (id % no_threads == thread_id)
						filter[id].emplace_back(j);
				}
			}
				});

		for (auto& t : thr_workers)
			t.join();
	}

	//if (params.verbosity_level > 1)
		cerr << "Filter size: " << no_items << endl;

	return true;
}

// ****************************************************************************
void CFilter::reorder_items(const vector<uint32_t>& reordering_map, uint32_t no_threads)
{
	if (filter.empty())
		return;

	cerr << "Reordering filter" << endl;

	assert(filter.size() == reordering_map.size());

	// Reorder rows
	{
		vector<vector<id_t>> reo_filter(filter.size());

		for (size_t i = 0; i < reordering_map.size(); ++i)
			reo_filter[reordering_map[i]] = move(filter[i]);

		filter = move(reo_filter);
	}

	// Reorder in rows
	vector<thread> thr_workers;
	thr_workers.reserve(no_threads);
	atomic<uint32_t> row_id = 0;

	for (uint32_t i = 0; i < no_threads; ++i)
		thr_workers.emplace_back([&] {
		while (true)
		{
			uint32_t my_row_id = row_id.fetch_add(64);
			if (my_row_id >= filter.size())
				break;

			uint32_t id_start = my_row_id;
			uint32_t id_end = min(id_start + 64, (uint32_t)filter.size());

			for (uint32_t id = id_start; id < id_end; ++id)
				for (auto& x : filter[id])
					x = reordering_map[x];
		}
		});

	for (auto& t : thr_workers)
		t.join();


}