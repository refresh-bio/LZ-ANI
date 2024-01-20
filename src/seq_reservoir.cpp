#include "seq_reservoir.h"

// ****************************************************************************
void CSeqReservoir::append(const string& name, const string& seq)
{
	uint8_t* ptr = (uint8_t*) mma.allocate(seq.length());

	for (size_t i = 0; i < seq.length(); ++i)
		ptr[i] = dna_code[seq[i]];

	uint32_t id = (uint32_t)items.size();

	items.emplace_back(name, ptr, seq.length(), 1);

	seq_id_map[name] = id;
}

// ****************************************************************************
bool CSeqReservoir::load_fasta(const vector<string>& fasta_files)
{

	return true;
}

// ****************************************************************************
bool CSeqReservoir::load_multifasta(const vector<string>& fasta_files)
{
	for (const auto& fn : fasta_files)
	{
		refresh::stream_in_file sif(fn, 16 << 20, 16 << 20);
		refresh::stream_decompression sdec(&sif);

		string line;
		string name, seq;

		// !!! Doda� test czy plik uda�o si� otworzy�

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
		}

		if (!name.empty())
			append(name, seq);
	}


	return true;
}
