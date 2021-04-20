/*
 * main.cpp
 *
 *  Created on: Jul 24, 2017
 *      Author: vlad
 */

#include "Particle.h"
#include "Experiment.h"
#include "argparse.hpp"
#include <unistd.h>
#include <thread>
#include <mutex>
#include <vector>

using namespace std;

int main(int argc,const char** argv) {
	time_t startTime, finishTime;
	ArgumentParser parser;
	string inputFileName;
	vector<thread> threads;

	// add some arguments to search for
	parser.addArgument("-f", "--file", 1);

	parser.parse(argc, argv);


/*	int rez=0;
	while ( (rez = getopt(argc,argv,"f:")) != -1){
		switch (rez){
		case 'f': inputFileName = optarg; break;
	       }
	}*/

	inputFileName = parser.retrieve<string>("file");
	if (inputFileName == "") inputFileName = "input.cfg";

	//cout << inputFileName << endl;

	startTime = time(NULL);

	Experiment e(inputFileName);

/*	for (int j = 0; j < e.getNumberOfPrints(); j++) {

		for(int i = 0; i < e.getNumberOfThreads(); i++) {
			thread thr([&e](){ e.runExperiment(); });
			threads.emplace_back(std::move(thr));
		}

		for(auto& thr : threads)
			thr.join();

		threads.clear();

		e.print();
	}
*/
	e.runExperiment();
	finishTime = time(NULL);
	cout<<"Evaluation time is: "<< finishTime - startTime << "sec." << endl;

	return 0;
}




