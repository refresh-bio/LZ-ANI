#include "seq_reservoir.h"

#include <iostream>

// ****************************************************************************
void CSeqReservoir::append(const string& name, const string& seq)
{
#ifdef USE_PACKED_SEQS
	uint8_t* ptr_seq = (uint8_t*) mma_seq.allocate((seq.length() + 1) / 2);

	for (size_t i = 0; i < seq.length() / 2; ++i)
		ptr_seq[i] = (uint8_t) ((dna_code[seq[2 * i]] << 4) + dna_code[seq[2 * i + 1]]);

	if (seq.length() & 1)
		ptr_seq[seq.length() / 2] = (uint8_t) (dna_code[seq[seq.length() - 1]] << 4);
#else
	uint8_t* ptr_seq = (uint8_t*)mma_seq.allocate(seq.length());

	for (size_t i = 0; i < seq.length(); ++i)
		ptr_seq[i] = dna_code[seq[i]];
#endif

	uint32_t id = (uint32_t)items.size();

	auto p = find(name.begin(), name.end(), ' ');

	char *ptr_name = (char*) mma_name.allocate(p - name.begin() + 1);
	copy(name.begin(), p, ptr_name);
	ptr_name[p - name.begin()] = 0;

	items.emplace_back(ptr_name, ptr_seq, seq.length(), 1);
}

// ****************************************************************************
bool CSeqReservoir::load_fasta(const vector<string>& fasta_files, uint32_t sep_len)
{
	cerr << "Loading sequences\n";

	for (const auto& fn : fasta_files)
	{
		auto buf_size = min<size_t>(filesystem::file_size(filesystem::path(fn)), 16ull << 20);

		refresh::stream_in_file sif(fn, buf_size, buf_size);

		if (!sif.is_open())
		{
			cerr << "Cannot open file: " << fn << endl;
			return false;
		}

		refresh::stream_decompression sdec(&sif);

		uint32_t no_parts = 0;
		string line;
		string seq;

		string separator(sep_len, code_N_seq);

		while (true)
		{
			if (sdec.getline(line) < 0)
				break;

			if (line.empty())
				continue;

			if (line.front() == '>')
			{
				++no_parts;
				if (!seq.empty())
					seq.append(separator);
			}
			else
				seq.append(line);
		}

		if (!seq.size())
			seq.resize(seq.size() - sep_len);

		string name = fn;

		if (name.ends_with(".gz"))
			name.resize(name.size() - 3);
		if(name.ends_with(".fna"))
			name.resize(name.size() - 4);
		else if(name.ends_with(".fasta"))
			name.resize(name.size() - 6);

		append(name, seq);

		if (items.size() % 100 == 0)
		{
			cerr << items.size() << "\r";
			fflush(stdout);
		}
	}

	cerr << items.size() << "\r";
	fflush(stdout);

	return true;
}

// ****************************************************************************
bool CSeqReservoir::load_multifasta(const vector<string>& fasta_files)
{
	cerr << "Loading sequences\n";

	for (const auto& fn : fasta_files)
	{
		refresh::stream_in_file sif(fn, 16 << 20, 16 << 20);

		if (!sif.is_open())
		{
			cerr << "Cannot open file: " << fn << endl;
			return false;
		}

		refresh::stream_decompression sdec(&sif);

		string line;
		string name, seq;

		while (true)
		{
			if (sdec.getline(line) < 0)
				break;

			if (line.empty())
				continue;

			if (line.front() == '>')
			{
				if (!name.empty())
					append(name, seq);
				name.assign(line.begin() + 1, line.end());
				seq.clear();
			}
			else
				seq.append(line);

			if (items.size() % 1000 == 0)
			{
				cerr << items.size() << "\r";
				fflush(stdout);
			}
		}

		if (!name.empty())
			append(name, seq);
	}

	cerr << items.size() << "\r";
	fflush(stdout);

	return true;
}

// ****************************************************************************
vector<uint32_t> CSeqReservoir::reorder_items()
{
	cerr << "Reordering sequences" << endl;

	vector<uint32_t> reordering_map;
	reordering_map.resize(items.size());

	vector<item_for_sorting_t> items_for_sorting;
	items_for_sorting.reserve(items.size());

	for (size_t i = 0; i < items.size(); ++i)
		items_for_sorting.emplace_back(items[i].name, items[i].len, items[i].no_parts, (uint32_t)i);

	stable_sort(items_for_sorting.begin(), items_for_sorting.end(), [](const auto& x, const auto& y)
		{
			auto x_size = x.len - x.no_parts * 2;
			auto y_size = y.len - y.no_parts * 2;
			if (x_size != y_size)
				return x_size > y_size;
			return x.name < y.name;
		});


	for (size_t i = 0; i < items_for_sorting.size(); ++i)
		reordering_map[items_for_sorting[i].id] = (uint32_t)i;

	vector<item_t> reo_items;
	reo_items.reserve(items.size());

	for (size_t i = 0; i < items_for_sorting.size(); ++i)
		reo_items.emplace_back(items[items_for_sorting[i].id]);

	items = move(reo_items);

	return reordering_map;
}
