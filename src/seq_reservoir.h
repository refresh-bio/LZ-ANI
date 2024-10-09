// *******************************************************************************************
// This file is a part of LZ-ANI software distributed under GNU GPL 3 license.
// The homepage of the LZ-ANI project is https://github.com/refresh-bio/LZ-ANI
//
// Copyright(C) 2024-2024, S.Deorowicz, A.Gudys
//
// Version: 1.1.0
// Date   : 2024-09-05
// *******************************************************************************************

#pragma once
#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <cinttypes>
#include <algorithm>
#include <string_view>

#include "../libs/refresh/allocators/lib/memory_monotonic.h"
#include "../libs/refresh/compression/lib/file_wrapper.h"

#include "defs.h"

using namespace std;

// #define USE_PACKED_SEQS

class seq_view
{
	const uint8_t* data;
	const uint32_t len;
	const internal_packing_t internal_packing;

public:
	seq_view(const uint8_t *data = 0, const uint32_t len = 0, internal_packing_t internal_packing = internal_packing_t::none) :
		data(data), 
		len(len), 
		internal_packing(internal_packing)
	{}

/*	void assign(uint8_t* _data, uint32_t _len, internal_packing_t _internal_packing = internal_packing_t::none)
	{
		data = _data;
		len = _len;
		internal_packing = _internal_packing;
	}*/

	uint8_t operator[](const uint32_t pos) const
	{
		switch (internal_packing)
		{
		case internal_packing_t::none:
			return data[pos];
			break;
		case internal_packing_t::two_in_byte:
			if (pos & 1)
				return data[pos / 2] & 0xf;
			else
				return data[pos / 2] >> 4;
			break;
		case internal_packing_t::three_in_byte:
			uint32_t pos_div_3 = pos / 3;

			switch (pos - 3 * pos_div_3)
			{
			case 0:
				return data[pos_div_3] / 36;
			case 1:
				return (data[pos_div_3] / 6) % 6;
			case 2:
				return data[pos_div_3] % 6;
			}
			break;
		}

		return 0;			// Never should be here
	}

	uint32_t size()	const
	{
		return len;
	}

/*	static void pack2(uint8_t* dest, uint8_t* src, uint32_t len)
	{
		for (uint32_t i = 0; i < len / 2; ++i)
			dest[i] = (src[2 * i] << 4) + src[2 * i + 1];

		if (len & 1)
			dest[len / 2] = src[len - 1] << 4;
	}

	static void unpack2(uint8_t* dest, uint8_t* src, uint32_t len)
	{
		for (uint32_t i = 0; i < len / 2; ++i)
		{
			dest[2 * i] = src[i] >> 4;
			dest[2 * i + 1] = src[i] & 0xf;
		}

		if (len & 1)
			dest[len - 1] = src[len / 2] >> 4;
	}

	static void pack3(uint8_t* dest, uint8_t* src, uint32_t len)
	{
		for (uint32_t i = 0; i < len / 3; ++i)
			dest[i] = 36 * src[3 * i] + 6 * src[3 * i + 1] + src[3 * i + 2];

		switch (len % 3)
		{
		case 2:
			dest[len / 3] = 6 * src[len - 2] + src[len - 1];
			break;
		case 1:
			dest[len / 3] = src[len - 1];
			break;
		case 0:
			// Nothing
			break;
		}
	}

	static void unpack3(uint8_t * dest, uint8_t * src, uint32_t len)
	{
		uint32_t len_div_3 = len / 3;

		for (uint32_t i = 0; i < len_div_3; ++i)
		{
			dest[3 * i] = src[i] / 36;
			dest[3 * i + 1] = src[i] / 6 - 6 * dest[3 * i];
			dest[3 * i + 2] = src[i] - 36 * dest[3 * i] - 6 * dest[3 * i + 1];
		}

		switch (len % 3)
		{
		case 2:
			dest[len - 2] = src[len_div_3] / 6;
			dest[len - 1] = src[len_div_3] - 6 * dest[len - 2];
			break;
		case 1:
			dest[len - 1] = src[len_div_3];
		case 0:
			// Nothing
			break;
		}
	}*/
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
	refresh::memory_monotonic_unsafe mma_seq;
	refresh::memory_monotonic_unsafe mma_name;

	internal_packing_t internal_packing = internal_packing_t::none;
	alphabet_t alphabet = alphabet_t::dna;

	array<uint8_t, 256> dna_code;

	void append(const string& name, const string& seq);

public:
	CSeqReservoir(internal_packing_t internal_packing = internal_packing_t::none /*, alphabet_t alphabet = alphabet_t::dna*/) :
		mma_seq(32 << 20, 16),
		mma_name(1 << 20, 16),
		internal_packing(internal_packing),
//		alphabet(alphabet),
		dna_code{}
	{
		if (alphabet == alphabet_t::dna)
		{
			fill(dna_code.begin(), dna_code.end(), code_N_seq);
			dna_code['A'] = dna_code['a'] = code_A;
			dna_code['C'] = dna_code['c'] = code_C;
			dna_code['G'] = dna_code['g'] = code_G;
			dna_code['T'] = dna_code['t'] = code_T;
		}
		else if (alphabet == alphabet_t::aminoacid)
		{
			assert(0);
		}
	}

	bool load_fasta(const vector<string>& fasta_files, uint32_t sep_len, uint32_t verbosity_level);
	bool load_multifasta(const vector<string>& fasta_files, uint32_t verbosity_level);

	void release()
	{
		items.clear();
		items.shrink_to_fit();
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

	string get_sequence_name(const size_t id)
	{
		return string(items[id].name);
	}

	uint32_t get_sequence_length(const size_t id)
	{
		return items[id].len;
	}

	vector<uint32_t> reorder_items(uint32_t verbosity_level);
};

// EOF