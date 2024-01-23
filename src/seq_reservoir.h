#pragma once
#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <cinttypes>
#include <algorithm>
#include <string_view>

#include "../libs/refresh/memory_monotonic.h"
#include "../libs/refresh/file_wrapper.h"

#include "defs.h"

using namespace std;

//#define USE_PACKED_SEQS

class seq_view
{
	const uint8_t* data;
	uint32_t len;

public:
	seq_view(const uint8_t *data = nullptr, const uint32_t len = 0) :
		data(data), len(len)
	{}

	uint8_t operator[](const uint32_t pos) const
	{
		return data[pos];
	}

	void assign(uint8_t* _data, uint32_t _len)
	{
		data = _data;
		len = len;
	}

	uint32_t size() const
	{
		return len;
	}
};

class packed_seq_view
{
	const uint8_t* data;
	uint32_t len;

public:
	packed_seq_view(const uint8_t *data = 0, const uint32_t len = 0) :
		data(data), len(len)
	{}

	void assign(uint8_t* _data, uint32_t _len)
	{
		data = _data;
		len = len;
	}

	uint8_t operator[](const uint32_t pos) const
	{
		if (pos & 1)
			return data[pos / 2] & 0xf;
		else
			return data[pos / 2] >> 4;
	}

	uint32_t size()	const
	{
		return len;
	}

	static void pack(uint8_t* dest, uint8_t* src, uint32_t len)
	{
		for (uint32_t i = 0; i < len / 2; ++i)
			dest[i] = (src[2 * i] << 4) + src[2 * i + 1];

		if (len & 1)
			dest[len / 2] = src[len - 1] << 4;
	}

	static void unpack(uint8_t * dest, uint8_t * src, uint32_t len)
	{
		for (uint32_t i = 0; i < len / 2; ++i)
		{
			dest[2 * i] = src[i] >> 4;
			dest[2 * i + 1] = src[i] & 0xf;
		}

		if (len & 1)
			dest[len - 1] = src[len / 2] >> 4;
	}
};

class CSeqReservoir
{
public:
	struct item_t
	{
		string_view name;
		const uint8_t* data;
		uint32_t len;
		uint32_t no_parts;

		item_t(const string_view name, const uint8_t* data, const uint32_t len, const uint32_t no_parts) :
			name(name),
			data(data),
			len(len),
			no_parts(no_parts)
		{}
	};

	struct item_for_sorting_t
	{
		string_view name;
		uint32_t len;
		uint32_t no_parts;
		uint32_t id;

		item_for_sorting_t(const string_view name, const uint32_t len, const uint32_t no_parts, uint32_t id) :
			name(name),
			len(len),
			no_parts(no_parts),
			id(id)
		{}
	};

private:
	vector<item_t> items;
	unordered_map<string, size_t> seq_id_map;
	refresh::memory_monotonic_unsafe mma_seq;
	refresh::memory_monotonic_unsafe mma_name;

	array<uint8_t, 256> dna_code;

	void append(const string& name, const string& seq);

public:
	CSeqReservoir() :
		mma_seq(16 << 20, 16),
		mma_name(1 << 20, 16),
		dna_code{}
	{
		fill(dna_code.begin(), dna_code.end(), code_N_seq);
		dna_code['A'] = code_A;
		dna_code['C'] = code_C;
		dna_code['G'] = code_G;
		dna_code['T'] = code_T;
	}

	bool load_fasta(const vector<string>& fasta_files, uint32_t sep_len);
	bool load_multifasta(const vector<string>& fasta_files);

	void release()
	{
		items.clear();
		items.shrink_to_fit();
		seq_id_map.clear();
		mma_seq.release();
		mma_name.release();
	}

	size_t size() const { return items.size(); } 
	bool empty() const { return items.empty(); }

	vector<item_t>::iterator get_sequence(const size_t id)
	{
		if (id >= items.size())
			return items.end();
		else
			return items.begin() + id;
	}

	vector<item_t>::iterator get_sequence(const string &name)
	{
		auto p = seq_id_map.find(name);

		if (p != seq_id_map.end())
			return items.begin() + p->second;
		else
			return items.end();
	}

	bool is_valid_iter(const vector<item_t>::iterator iter)
	{
		return iter != items.end();
	}

	vector<string> get_sequence_names()
	{
		vector<string> vs;

		vs.reserve(items.size());

		for (const auto& x : items)
			vs.emplace_back(x.name);

		return vs;
	}

	vector<uint32_t> reorder_items();

};