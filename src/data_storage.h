#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <cstring>

using namespace std;

class CSample
{
	string name;
	size_t size;
	char* data;
	char* ptr;
	char* end;

public:
	CSample(const string& name, size_t size, char* data) :
		name(name),
		size(size),
		data(data),
		ptr(data),
		end(data + size)
	{}

	bool eof()
	{
		return ptr >= end;
	}

	char getc()
	{
		return *ptr++;
	}

	void restart()
	{
		ptr = data;
	}

	size_t item_size()
	{
		return size;
	}

	size_t pos()
	{
		return ptr - data;
	}

	string get_name()
	{
		return name;
	}
};

struct SDataStorageItem
{
	string name;
	size_t size;
	char* data;

	SDataStorageItem(const string& name, size_t size, char* data) :
		name(name),
		size(size),
		data(data)
	{}
};


class CDataStorage
{
	bool buffered;

	vector<SDataStorageItem> items;
	unordered_map<string, size_t> item_map;

	void release()
	{
		for (auto& x : items)
			if (x.data)
				delete[] x.data;

		items.clear();
		item_map.clear();
	}

	CSample* _open(size_t fid);
	pair<char*, size_t> pack_seq(string& id, const vector<char>& vc);

public:
	CDataStorage(bool buffered) :
		buffered(buffered)
	{}

	~CDataStorage()
	{
		release();
	}

	bool prepare_from_many_files(const vector<string>& file_names);
	bool prepare_from_multi_fasta(const string& file_name);

	CSample* open(const string& fn);
	CSample* open(size_t fid);
	void close(CSample* sample);
	void close_hard(CSample* sample);
	void free_buffered();

	void close()
	{
		release();
	}

	size_t size()
	{
		return items.size();
	}

	bool empty()
	{
		return items.empty();
	}

	vector<SDataStorageItem>::iterator begin()
	{
		return items.begin();
	}

	vector<SDataStorageItem>::const_iterator cbegin()
	{
		return items.cbegin();
	}

	vector<SDataStorageItem>::iterator end()
	{
		return items.end();
	}

	vector<SDataStorageItem>::const_iterator cend()
	{
		return items.cend();
	}
};


// EOF
