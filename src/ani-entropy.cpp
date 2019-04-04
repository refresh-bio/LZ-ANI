// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

#include "defs.h"
#include "files.h"
#include "worker.h"

#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std::chrono;

CWorker worker;



// ****************************************************************************
int main(int argc, char **argv)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

/*	cout << "*** Loading...";
	if (!load_data(argc, argv))
	{
		return 0;
	}
	cout << " ok\n";

	cout << "*** HT preparation...";
	prepare_ht();
	cout << " ok\n";

	cout << "*** Parsing...";
	parse();
	parsing_postprocess();
	export_parsing();

	cout << " ok\n";

	calc_ani();*/

	cout << "*** End of work\n";

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Running time: " << time_span.count() << " seconds.\n";

//	getchar();

	return 0;
}
