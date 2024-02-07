#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include "params.h"

const std::string LZ_ANI_VER = "lz_ani 0.4";
const std::string LZ_ANI_DATE = "2024-02-07";
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
