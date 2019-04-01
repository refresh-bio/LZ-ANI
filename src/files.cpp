#include "files.h"
#include <cstdio>

// ****************************************************************************
bool load_file(const string &file_name, seq_t &seq, uint32_t &n_parts)
{
	seq.clear();
	
	FILE *f = fopen(file_name.c_str(), "rb");
	if (!f)
		return false;

	int c;
	bool is_comment = false;

	n_parts = 0;

	while ((c = getc(f)) != EOF)
	{
		if (c == '>')
			is_comment = true;
		else
		{
			if (c == '\n' || c == '\r')
			{
				if (is_comment)
				{
					is_comment = false;
					if(!seq.empty())
						for (int i = 0; i < close_dist; ++i)
							seq.push_back(sym_N);
					++n_parts;
				}
			}
			else if (!is_comment)
				seq.push_back((uint8_t)c);
		}
	}

	return true;
}

// ****************************************************************************
void duplicate_rev_comp(seq_t &seq)
{
	int size = (int) seq.size();

	for (int i = size - 1; i >= 0; ++i)
	{
		if (seq[i] == sym_N)
			seq.push_back(sym_N);
		else if (seq[i] == sym_A)
			seq.push_back(sym_T);
		else if (seq[i] == sym_C)
			seq.push_back(sym_G);
		else if (seq[i] == sym_G)
			seq.push_back(sym_C);
		else if (seq[i] == sym_T)
			seq.push_back(sym_A);
	}

}
