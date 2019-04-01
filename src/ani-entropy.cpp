// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

#include "defs.h"
#include "files.h"

#include <iostream>


seq_t s_reference;
seq_t s_data;

uint32_t n_reference;
uint32_t n_data;

uint32_t ht_size;
uint32_t ht_mask;
double ht_max_fill_factor = 0.5;
vector<uint32_t> ht;
vector<CFactor> v_parsing;

bool load_data(int argc, char **argv);
void parse();
void prepare_ht();
uint32_t my_hash(seq_t::iterator p, uint32_t len);

// ****************************************************************************
uint32_t my_hash(seq_t::iterator p, uint32_t len)
{
	uint32_t r = 0;

	for (uint32_t i = 0; i < len; ++i)
	{
		if (*(p + i) == sym_N)
			return HT_FAIL;

		r = r * 127 + *(p + i);
	}

	return r & ht_mask;
}

// ****************************************************************************
bool load_data(int argc, char **argv)
{
	if (argc < 3)
	{
		cerr << "Usage: ani-entropy <ref_name> <data_name>\n";
		return false;
	}

	if (!load_file(string(argv[1]), s_reference, n_reference))
	{
		cerr << "Cannot load " << argv[1] << endl;
		return false;
	}

	if (!load_file(string(argv[2]), s_data, n_data))
	{
		cerr << "Cannot load " << argv[2] << endl;
		return false;
	}

//	duplicate_rev_comp(s_reference);

	return true;
}

// ****************************************************************************
void parse()
{

}

// ****************************************************************************
void prepare_ht()
{
	uint32_t x = (uint32_t) (s_reference.size() / ht_max_fill_factor);

	while (x & (x - 1))
		x &= x - 1;

	ht_size = 2 * x;
	ht_mask = ht_size - 1;

	ht.clear();
	ht.resize(ht_size, HT_EMPTY);

	for (size_t i = 0; i + MIN_MATCH_LEN < s_reference.size(); ++i)
	{
		auto ht_idx = my_hash(s_reference.begin() + i, MIN_MATCH_LEN);

		if (ht_idx != HT_FAIL)
		{
			while (ht[ht_idx] != HT_EMPTY)
				ht_idx = (ht_idx + 1) & ht_mask;

			ht[ht_idx] = i;
		}

	}
}

// ****************************************************************************

int main(int argc, char **argv)
{
	if (!load_data(argc, argv))
	{
		return 0;
	}

	prepare_ht();

	cout << "End of work\n";

	return 0;
}
