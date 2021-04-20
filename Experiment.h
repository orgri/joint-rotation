/*
 * Experiment.h
 *
 *  Created on: Sep 27, 2017
 *      Author: cuda
 */

#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include <fstream>
#include <deque>
#include <map>
#include <iterator>
#include <limits>
#include "Particle.h"
#include "ConfigFile.h"

using namespace std;

namespace Exprm {
	enum class Task{SimpleRun, PowerLoss};
	enum class Variable{Amplitude, Frequency, NotSet};
	enum class NumericalMethod{RK2, RK4, Euler};
};

struct Angles {
	double theta;
	double phi;
	double energy;
};

struct CharactersOfMotion {
	double mTheta;
	double periodTheta;
	double velocityPhi;
	double phasePhi;
	double energy;
};

class Experiment {
public:
	Experiment(const string&);
	virtual ~Experiment();

	int getNumberOfThreads() {return this->numberOfThreads;};
	void print();
	void printLast();
	void printConsole();
	void printTask();
	void printError();
	void printDebugInfo();
	void printHeader();
	void printCharacters(const CharactersOfMotion &data, string fileName);
	void runExperiment();
	void checkExtremum();
	string getVariable();
	string getTask();
	string getNumMethod();

private:
	void saveConf();
	void setParams();
	void setVar();
	void setInitialState();

	void addExtremumTo(deque<pair<double,Angles>> &data, Angles extremun);
	deque<pair<double,Angles>> processMaxima(const deque<pair<double,Angles>> &data);
	deque<pair<double,Angles>> processMinimum(const deque<pair<double,Angles>> &data);
	string findMode();

	void relaxation(int nPeriods);
	void run(int nPeriods);
	void calcPeriod();
	void calcAverageTheta(int nPeriods);
	void simpleRun();
	void powerLoss();
	void hysteresis();
	void clearExtremums();

	double getEnergyPerTime(/*double cutTime, double cutEnergy*/);
	double getEnergyPerTime(const deque<pair<double,Angles>> &data);
	CharactersOfMotion calcCharacters(const deque<pair<double,Angles>> &data);
	double calcAveragePeriod(const deque<pair<double,Angles>> &data);
	bool compare(const deque<double> &data);
	bool compare(const deque<pair<double,Angles>> &data);
	bool compareLagAngle(const deque<pair<double,Angles>> &data);
	bool compareEnergy(const deque<pair<double,Angles>> &data);
	bool relativeToleranceCompare(double x, double y);
	bool compareDiffTime(const deque<pair<double,Angles>> data);


	long long int numberOfPeriods;
	long long int maxNumberOfPeriods;
	long long int skipPeriods;
	long long int periodsForRelaxation;
	int printLastNperiods;
	unsigned int numberOfExtremums;
	int numberOfThreads;
	int outputPrecision;

	Exprm::NumericalMethod itsMethod;
	Exprm::Task itsTask;
	Exprm::Variable itsVar;

	double initialSmallTheta;
	double initialSmallPhi;
	double initialBigTheta;
	double initialBigPhi;

	double averageSmallTheta;
	double averageBigTheta;


	double cutEnergy;
	double cutTime;

	double period;
	double printStep;
	double timeStep;
	double acceptableDiff;
	double var;
	double varStep;
	double varMax;
	long long int nRuns;

	bool debug;
	bool intelegentStop;

	Particle *p;

	deque<pair<double,Angles>> maxSmallTheta;
	deque<pair<double,Angles>> maxBigTheta;
	deque<pair<double,Angles>> minSmallTheta;
	deque<pair<double,Angles>> minBigTheta;
	string inputFileName;
	string outputFileName;
	string fileDirectory;

	fstream output;
//	fstream outputTask;
	fstream outputDebug;
};

#endif /* EXPERIMENT_H_ */
