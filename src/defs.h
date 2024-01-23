#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include "params.h"

using namespace std;

const uint8_t code_A = 0;
const uint8_t code_C = 1;
const uint8_t code_G = 2;
const uint8_t code_T = 3;
const uint8_t code_N = 4;
const uint8_t code_N_ref = 4;
const uint8_t code_N_seq = 5;

enum class flag_t {match, match_close, match_distant, literal, run_literals, match_literal};

const int HT_EMPTY = -1;
const int HT_FAIL = -1;

using id_t = uint32_t;
using seq_t = vector<uint8_t>;

struct factor_t {
	int data_pos;
	flag_t flag;
	int offset;
	int len;
	uint8_t symbol;

	factor_t(int _data_pos, flag_t _flag, uint32_t _offset, uint32_t _len, uint8_t _symbol) :
		data_pos(_data_pos), flag(_flag), offset(_offset), len(_len), symbol(_symbol)
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

struct file_desc_t {
	string file_name;
	string seq_name;
	size_t file_size;
	size_t seq_size;
	size_t n_parts;
	seq_t data;

	file_desc_t(string file_name = "", string seq_name = "", size_t file_size = 0, size_t seq_size = 0, size_t n_parts = 0) :
		file_name(file_name),
		seq_name(seq_name),
		file_size(file_size),
		seq_size(seq_size),
		n_parts(n_parts)
	{}
};


// EOF
