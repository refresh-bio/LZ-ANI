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

#include <cstdint>
#include <vector>
#include <string>
#include "params.h"

const std::string LZ_ANI_VER = "lz-ani 1.1.0";
const std::string LZ_ANI_DATE = "2024-09-05";
const std::string LZ_ANI_AUTHORS = "Sebastian Deorowicz, Adam Gudys";
const std::string LZ_ANI_INFO = LZ_ANI_VER + " (" + LZ_ANI_DATE + ") by " + LZ_ANI_AUTHORS;

const uint8_t code_A = 0;
const uint8_t code_C = 1;
const uint8_t code_G = 2;
const uint8_t code_T = 3;
const uint8_t code_N = 4;
const uint8_t code_N_ref = 4;
const uint8_t code_N_seq = 5;

enum class flag_t {match, match_close, match_distant, literal, run_literals, match_literal};

using id_t = uint32_t;
using seq_t = vector<uint8_t>;

struct factor_t {
	int data_pos;
	flag_t flag;
	int offset;
	int len;

	factor_t(int _data_pos, flag_t _flag, uint32_t _offset, uint32_t _len) :
		data_pos(_data_pos), flag(_flag), offset(_offset), len(_len)
	{}
};

struct results_t
{
	int sym_in_matches;
	int sym_in_literals;
	int no_components;

	results_t() :
		sym_in_matches(0),
		sym_in_literals(0),
		no_components(0)
	{}

	results_t(int sym_in_matches, int sym_in_literals, int no_components) :
		sym_in_matches(sym_in_matches),
		sym_in_literals(sym_in_literals),
		no_components(no_components)
	{}
};

struct region_t
{
	int ref_start;
	int ref_end;
	int seq_start;
	int seq_end;
	int num_matches;
	int num_mismatches;

	region_t() :
		ref_start(-1),
		ref_end(-1),
		seq_start(-1),
		seq_end(-1),
		num_matches(0),
		num_mismatches(0)
	{}

	region_t(int ref_start, int ref_end, int seq_start, int seq_end, int num_matches, int num_mismatches) :
		ref_start(ref_start),
		ref_end(ref_end),
		seq_start(seq_start),
		seq_end(seq_end),
		num_matches(num_matches),
		num_mismatches(num_mismatches)
	{}

	void clear()
	{
		ref_start = -1;
		ref_end = -1;
		seq_start = -1;
		seq_end = -1;
		num_matches = 0;
		num_mismatches = 0;
	}

	bool empty()	const
	{
		return seq_start == seq_end;
	}

	int length()	const
	{
		return seq_end - seq_start;
	}

	void update_ref_start(int cand_ref_start)
	{
		if (ref_start < 0 || cand_ref_start < ref_start)
			ref_start = cand_ref_start;
	}

	void update_ref_end(int cand_ref_end)
	{
		if (ref_end < 0 || cand_ref_end > ref_end)
			ref_end = cand_ref_end;
	}

	void update_seq_start(int cand_seq_start)
	{
		if (seq_start < 0 || cand_seq_start < seq_start)
			seq_start = cand_seq_start;
	}

	void update_seq_end(int cand_seq_end)
	{
		if (seq_end < 0 || cand_seq_end > seq_end)
			seq_end = cand_seq_end;
	}

	void extend_region(int len)
	{
		ref_end += len;
		seq_end += len;
	}

	void update_matches(int num)
	{
		num_matches += num;
	}

	void update_mismatches(int num)
	{
		num_mismatches += num;
	}
};

struct id_results_t
{
	uint32_t id;
	results_t results;

	id_results_t() :
		id(0),
		results()
	{}
	
	id_results_t(uint32_t id, const results_t& results) :
		id(id),
		results(results)
	{}

	bool operator<(const id_results_t& rhs)
	{
		return id < rhs.id;
	}
};

using vec_id_results_t = vector<id_results_t>;

// EOF
