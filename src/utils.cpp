#include "utils.h"

vector<string> split(const string& str, char sep)
{
	vector<string> parts;

	string s;

	for (auto c : str)
	{
		if (c == sep)
		{
			parts.emplace_back(s);
			s.clear();
		}
		else
			s.push_back(c);
	}

	if (!s.empty())
		parts.emplace_back(s);

	return parts;
}

string remove_path_from_file(const string& file_path) 
{
	size_t pos = file_path.find_last_of("/\\");

	if (pos != string::npos)
		return file_path.substr(pos + 1);
	else
		return file_path;
}