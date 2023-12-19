#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "data_storage.h"

// ****************************************************************************
bool CDataStorage::prepare_from_many_files(const vector<string>& file_names)
{
	if (!items.empty())
		release();

	bool all_ok = true;
	std::error_code ec;

	items.reserve(file_names.size());
	item_map.reserve(file_names.size());

	for (const auto& fn : file_names)
	{
		auto fs = std::filesystem::file_size(std::filesystem::path(fn), ec);

		if (!ec)
		{
			char* ptr = nullptr;

			if (buffered)
			{
				FILE* f = fopen(fn.c_str(), "rb");
				if (!f)
				{
					cerr << "Cannot open: " << fn << endl;
					return false;
				}

				char* ptr = new char[fs];
				fread(ptr, 1, fs, f);
				fclose(f);
			}

			item_map[fn] = items.size();
			items.emplace_back(fn, fs, ptr);
		}
		else
		{
			cerr << "Input file " << fn << " : " << ec.message() << endl;
			return false;
		}
	}

	return true;
}

// ****************************************************************************
pair<char*, size_t> CDataStorage::pack_seq(string& id, const vector<char>& vc)
{
	size_t f_size = 2 + id.size() + vc.size();

	char* ptr = new char[f_size];

	ptr[0] = '>';
	memcpy(ptr + 1, id.data(), id.size());
	ptr[id.size() + 1] = '\n';
	memcpy(ptr + id.size() + 2, vc.data(), vc.size());

	auto p = find(id.begin(), id.end(), ' ');
	if (p != id.end())
		id.erase(p, id.end());

	return make_pair(ptr, f_size);
}

// ****************************************************************************
bool CDataStorage::prepare_from_multi_fasta(const string& file_name)
{
	if (!items.empty())
		release();

	std::error_code ec;
	auto fs = std::filesystem::file_size(std::filesystem::path(file_name), ec);

	if(ec)
	{
		cerr << "Input file " << file_name << " : " << ec.message() << endl;
		return false;
	}

	FILE* f = fopen(file_name.c_str(), "rb");

	if(!f)
	{
		cerr << "Cannot open: " << file_name << endl;
		return false;
	}

	buffered = true;		// Data from multi-FASTA are always buffered

	string id;
	vector<char> vc;

	setvbuf(f, nullptr, _IOFBF, min(fs, (size_t) (16 << 20)));

	int c;
	bool is_id = false;
	size_t f_size;
	char* ptr;

	while ((c = getc(f)) != EOF)
	{
		if (c == '>')
		{
			is_id = true;

			if (!id.empty())
			{
				tie(ptr, f_size) = pack_seq(id, vc);

				item_map[id] = items.size();
				items.emplace_back(id, f_size, ptr);

				id.clear();
				vc.clear();
			}
		}
		else
		{
			if (c == '\n' || c == '\r')
				is_id = false;
			else if (is_id)
				id.push_back((char) c);
			else
				vc.emplace_back((char) c);
		}
	}

	if (!id.empty() && !vc.empty())
	{
		tie(ptr, f_size) = pack_seq(id, vc);

		item_map[id] = items.size();
		items.emplace_back(id, f_size, ptr);
	}

	fclose(f);

	items.shrink_to_fit();

	return true;
}

// ****************************************************************************
CSample* CDataStorage::open(const string& fn)
{
	auto p = item_map.find(fn);
	if (p == item_map.end())
		return nullptr;

	return _open(p->second);
}

// ****************************************************************************
CSample* CDataStorage::open(size_t fid)
{
	if (fid >= items.size())
		return nullptr;

	return _open(fid);
}

// ****************************************************************************
CSample* CDataStorage::_open(size_t fid)
{
	auto& item = items[fid];

	if (!item.data)
	{
		FILE* f = fopen(item.name.c_str(), "rb");
		if (!f)
			return nullptr;

		item.data = new char[item.size];
		fread(item.data, 1, item.size, f);

		fclose(f);
	}

	return new CSample(item.name, item.size, item.data);
}

// ****************************************************************************
void CDataStorage::close(CSample* sample)
{
	if (sample)
		delete sample;
}

// ****************************************************************************
// Release memory for file
void CDataStorage::close_hard(CSample* sample)
{
	if (!buffered)
	{
		auto p = item_map.find(sample->get_name());

		assert(p != item_map.end());

		delete[] items[p->second].data;
		items[p->second].data = nullptr;
	}

	delete sample;
}

// ****************************************************************************
void CDataStorage::free_buffered()
{
	if (buffered)
		return;

	for(auto &x : items)
		if (x.data)
		{
			delete[] x.data;
			x.data = nullptr;
		}
}

// EOF