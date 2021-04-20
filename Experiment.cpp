/*
 * Experiment.cpp
 *
 *  Created on: Sep 27, 2017
 *      Author: cuda
 */

#include "Experiment.h"


Experiment::Experiment(const std::string& confFileName) {
	// TODO Auto-generated constructor stub
	inputFileName = confFileName;
	ConfigFile config(inputFileName);

	initialSmallTheta = config.read<double>("smallTheta",0.01);
	initialSmallPhi= config.read<double>("smallPhi",0.01);
	initialBigTheta = config.read<double>("bigTheta",0.01);
	initialBigPhi = config.read<double>("bigPhi",0.01);

	p = new Particle (initialSmallTheta, initialSmallPhi, initialBigTheta, initialBigPhi);

	this->setParams();

	this->output.open(fileDirectory + "log_"+ outputFileName, std::fstream::app | std::fstream::out);
	this->outputDebug.open(fileDirectory + "error.log", std::fstream::app | std::fstream::out);

	this->saveConf();
}

Experiment::~Experiment() {
	// TODO Auto-generated destructor stub
}

Exprm::Task readTask(ConfigFile& conf)
{
	using namespace Exprm;
	string readTask;
	readTask = conf.read<string>("task", "SimpleRun");

	if (readTask == "SimpleRun") return Task::SimpleRun;
	else if (readTask == "PowerLoss") return Task::PowerLoss;

	else return Task::SimpleRun;	//todo - throw exception
}

Exprm::NumericalMethod readMethod(ConfigFile& conf) {
	using namespace Exprm;
	string readMethod;
	readMethod = conf.read<string>("numericalMethod", "RK4");
	if (readMethod == "Euler") return NumericalMethod::Euler;
	else if (readMethod == "RK2") return NumericalMethod::RK2;
	else if (readMethod == "RK4") return NumericalMethod::RK4;

	else return NumericalMethod::RK4;
}

Exprm::Variable readVariable(ConfigFile& conf) {
	using namespace Exprm;
	string readVariable;
	readVariable = conf.read<string>("variable", "NotSet");

	if (readVariable == "Amplitude") return Variable::Amplitude;
	else if (readVariable == "Frequency") return Variable::Frequency;
	else if (readVariable == "NotSet") return Variable::NotSet;

	return Variable::NotSet;
}

void Experiment::setParams () {
	ConfigFile config(inputFileName);

	itsTask = readTask(config);
	itsMethod = readMethod(config);
	itsVar = readVariable(config);

	debug = config.read<bool>("debug", false);
	intelegentStop = config.read<bool>("intelegentStop", true);

	p->setMagnetization(config.read<double>("magnetization", 338.0));
	p->setAmplitude(config.read<double>("amplitude", 0.1));
	p->setFrequency(config.read<double>("frequency", 0.1));
	p->setViscosity(config.read<double>("viscosity", 0.001));
	p->setAlpha(config.read<double>("alpha", 0.1));
	p->setZAmplitude(config.read<double>("Hz", 0.0));
	p->setRo(signbit(config.read<int>("hySign", 1)) ? -1 : 1);

	p->setLinearPolarization(config.read<bool>("linearPolarization", false));
	p->setThermalBath(config.read<bool>("thermalBath", false));
	p->resetTime();

	outputFileName = config.read<string>("fileName", "output");
	fileDirectory = config.read<string>("fileDirectory", string(getenv("HOME")) + "/");

	numberOfPeriods = config.read<double>("numberOfPeriods", 1000);
	maxNumberOfPeriods = config.read<double>("maxNumberOfPeriods", 1e6);
	numberOfExtremums = config.read<double>("numberOfExtremums", 1000);
	skipPeriods = config.read<double>("skipPeriods", 0);
	periodsForRelaxation = config.read<double>("periodsForRelaxation", 0);
	printStep = config.read<double>("printStep", 1000);
	printLastNperiods = config.read<int>("printLastNperiods", 10);
	timeStep = config.read<double>("timeStep", 1000);
	acceptableDiff = config.read<double>("acceptableDiff", 1e-5);
	outputPrecision = config.read<int>("outputPrecision", 10);

	var = config.read<double>("variableMin", 0.0);
	varMax = config.read<double>("variableMax", 1.0);
	varStep = config.read<double>("variableStep", 0.1);

	cutEnergy = 0.0;
	cutTime = 0.0;
}

void Experiment::runExperiment () {
	printHeader();
	if (this->itsTask == Exprm::Task::SimpleRun)
		simpleRun();
	if (this->itsTask == Exprm::Task::PowerLoss)
		powerLoss();
}


void Experiment::run (int nPeriods) {

	switch (itsMethod) {

	case Exprm::NumericalMethod::Euler:
		for (int i = 0; i < nPeriods * nRuns; ++i) {
			p->EulerStep();
			checkExtremum();
			p->calcEnergy();
			p->applyDeltas();
			printError();
			if (i % int(printStep * nRuns) == 0){
				if(i % int(1000 * nRuns) == 0)
					printConsole();
				if (p->getTime() / period > skipPeriods)
					print();
			}
		}
		break;

	case Exprm::NumericalMethod::RK2:
		for (int i = 0; i < nPeriods * nRuns; ++i) {
			p->RungeKutta2Step();
			checkExtremum();
			p->calcEnergy();
			p->applyDeltas();
			printError();
			if (i % int(printStep * nRuns) == 0){
				if(i % int(1000 * nRuns) == 0)
					printConsole();
				if (p->getTime() / period > skipPeriods)
					print();
			}
		}
		break;

	case Exprm::NumericalMethod::RK4:
		for (long long int i = 0; i < nPeriods * nRuns; ++i) {
			p->RungeKutta4Step();
			checkExtremum();
			p->calcEnergy();
			p->applyDeltas();
			const double currentPeriod = ceil(p->getTime() / period);
			printError();
			if (i % int(printStep * nRuns) == 0){
				if(i % int(1000 * nRuns) == 0)
					printConsole();
				if (printLastNperiods == 0 && currentPeriod > skipPeriods)
					print();
			}
		}
		break;
	}
}

void Experiment::powerLoss() {
}

void Experiment::simpleRun() {
	for(; var <= varMax; var += varStep) {
		p->setInitialState(initialSmallTheta, initialSmallPhi, initialBigTheta, initialBigPhi);
		this->setVar();
		this->calcPeriod();

		if (intelegentStop) {
			bool isAcceptable = false;

			do {
				run(numberOfPeriods);
			} while (maxSmallTheta.size() < numberOfExtremums);

			do {
				run(numberOfPeriods);

				isAcceptable = this->compareEnergy(maxSmallTheta)
							&& this->compare(this->processMaxima(maxSmallTheta))
							&& this->compare(this->processMaxima(maxBigTheta));
			} while (isAcceptable == false && ceil(p->getTime() / period) < maxNumberOfPeriods);

			if (isAcceptable)
				std::cout << "Stoped due to fit in acceptable difference of energy and theta" << std::endl;

			if (p->getTime() / period > maxNumberOfPeriods)
				std::cout << "Stoped due to reach maxNumberOfPeriods" << std::endl;
		} else
			run(maxNumberOfPeriods);

		calcAverageTheta(numberOfExtremums);
		printTask();
		if (debug)	printDebugInfo();
		clearExtremums();
	}
}

void Experiment::calcPeriod() {
	period = 2.0 * M_PI / p->getFrequency();
	p->setTimeStep(period / timeStep);
	nRuns = period / p->getTimeStep();

	if (debug){
		std::cout << "period:  " << this->period;
		std::cout << "\ttimeStep:  " << p->getTimeStep() <<std::endl;
	}
}

void Experiment::calcAverageTheta(int nPeriods) {

	int steps = nRuns * nPeriods;

	for (long long int i = 0; i < steps; ++i) {
		p->RungeKutta4Step();
		p->calcEnergy();
		p->applyDeltas();
		printError();
		this->averageSmallTheta += p->getSmallTheta();
		this->averageBigTheta += p->getBigTheta();

		if (i < nRuns * printLastNperiods && i % int(printStep * nRuns) == 0) {
			print();
		}
	}

	this->averageSmallTheta = this->averageSmallTheta / steps;
	this->averageBigTheta = this->averageBigTheta / steps;

	if (debug) {
		std::cout << "==================== Average of "<< nPeriods <<" periods ===================" << std::endl;
		std::cout << "<smallTheta>:  " << this->averageSmallTheta;
		std::cout << "\t<bigTheta>:  " << this->averageBigTheta <<std::endl;
	}
}

CharactersOfMotion Experiment::calcCharacters(const deque<pair<double,Angles>> &data) {
	CharactersOfMotion params;
	params.mTheta = 0.0;

	for (auto& max : data) {
		params.mTheta += max.second.theta;
	}

	params.mTheta = params.mTheta / data.size();
	params.periodTheta = (data.back().first - data.front().first) / (data.size() - 1);
	params.velocityPhi = (data.back().second.phi - data.front().second.phi) / (data.back().first - data.front().first);
	params.phasePhi = data.back().second.phi - params.velocityPhi * data.back().first;
	params.energy = (data.back().second.energy - data.front().second.energy) / (data.back().first - data.front().first);

	return params;
}

double Experiment::calcAveragePeriod(const deque<pair<double,Angles>> &data) {
	double summ = 0;
	double diffTime = 0;
	double preTime = data.front().first;

	for (auto& max : data){
		diffTime = max.first - preTime;
		preTime = max.first;
		summ += diffTime;
	}
	return summ / (data.size() - 1);
}

void Experiment::addExtremumTo(deque<pair<double,Angles>> &data, Angles extremun) {
	data.push_back(std::make_pair(p->getTime(), extremun));
	if (data.size() > numberOfExtremums)
		data.pop_front();
}

void Experiment::checkExtremum() {
	static double prevdSmallTheta;
	static double prevdBigTheta;
	Angles small = {p->getSmallTheta(), p->getSmallPhi(), p->getEnergy()};
	Angles big = {p->getBigTheta(), p->getBigPhi(), p->getEnergy()};

	if (prevdSmallTheta > 0 && p->getDeltaSmallTheta() < 0)
		addExtremumTo(maxSmallTheta, small);
	if (prevdBigTheta > 0 && p->getDeltaBigTheta() < 0)
		addExtremumTo(maxBigTheta, big);
	if (prevdSmallTheta < 0 && p->getDeltaSmallTheta() > 0)
		addExtremumTo(minSmallTheta, small);
	if (prevdBigTheta < 0 && p->getDeltaBigTheta() > 0)
		addExtremumTo(minBigTheta, big);

	prevdSmallTheta = p->getDeltaSmallTheta();
	prevdBigTheta = p->getDeltaBigTheta();
}

bool Experiment::compare(const deque<double> &data) {
	int id = data.size() % 2;
	return ((fabs((data[id] - data.front()) / data[id]) < acceptableDiff) &&
			(fabs((data.back() - data[id]) / data.back()) < acceptableDiff) &&
			(fabs((data.front() - data.back())) / data.front() < acceptableDiff) );
}

bool Experiment::compare(const deque<pair<double,Angles>> &data) {
	int id = data.size() / 2;
	return ( (fabs((data[id].second.theta - data.front().second.theta) / data[id].second.theta) < acceptableDiff) &&
			(fabs((data.back().second.theta - data[id].second.theta) / data.back().second.theta) < acceptableDiff) &&
			(fabs((data.front().second.theta - data.back().second.theta) / data.front().second.theta) < acceptableDiff) );
}

bool Experiment::compareLagAngle(const deque<pair<double,Angles>> &data) {
	int id = data.size() / 2;
	double front = data.front().second.phi - p->getFrequency() * data.front().first;
	double middle = data[id].second.phi  - p->getFrequency() * data[id].first;
	double back = data.back().second.phi  - p->getFrequency() * data.back().first;

	return fabs((front - middle) / front) < acceptableDiff &&
			fabs((middle - back) / middle) < acceptableDiff &&
			fabs((back - front) / back) < acceptableDiff;
}

bool Experiment::compareEnergy(const deque<pair<double,Angles>> &data) {
	int id = data.size() / 2;
	double firstDiff = (data.back().second.energy - data.front().second.energy) / (data.back().first - data.front().first);
	double secondDiff = (data.back().second.energy - data[id].second.energy) / (data.back().first - data[id].first);
	double thirdDiff = (data[id].second.energy - data.front().second.energy) / (data[id].first - data.front().first);

	if (debug) {
		std::cout << "firstDiff:  "<< fabs((firstDiff - secondDiff) / firstDiff);
		std::cout << " secondDiff:  " << fabs((secondDiff - thirdDiff) / secondDiff);
		std::cout << " thirdDiff:  " << fabs((thirdDiff - firstDiff) / thirdDiff) << std::endl;
	}

	return fabs((firstDiff - secondDiff) / firstDiff) < acceptableDiff &&
			fabs((secondDiff - thirdDiff) / secondDiff) < acceptableDiff &&
			fabs((thirdDiff - firstDiff) / thirdDiff) < acceptableDiff;
}

bool Experiment::relativeToleranceCompare(double x, double y) {
	double epsilon = p->getTimeStep();
    double maxXY = std::max( std::fabs(x) , std::fabs(y) ) ;
    return std::fabs(x - y) <= epsilon * maxXY;
}

deque<pair<double,Angles>> Experiment::processMaxima(const deque<pair<double,Angles>> &data) {
	deque<pair<double,Angles>> temp = data;
	deque<pair<double,Angles>> maximas;

	do {
		double preDiff = 0;
		double preVal = 0;

		std::pair<double,Angles> preMax;

		maximas.clear();
		maximas = temp;
		temp.clear();

		for(auto& max : maximas) {
			double diff;
			if (max.first == maximas.front().first) {
				preVal = max.second.theta;
				preMax = max;
			} else {
				diff = max.second.theta - preVal;
				double error = fabs(diff / max.second.theta);
				//std::cout << "diff: " << diff << std::endl;
				if (preDiff > 0 && diff < 0 && error > acceptableDiff)
					temp.push_back(preMax);

				preVal = max.second.theta;
				preDiff = diff;
				preMax = max;
			}
		}
	} while (temp.size() > 5);

	return maximas;
}

deque<pair<double,Angles>> Experiment::processMinimum(const deque<pair<double,Angles>> &data){
	deque<pair<double,Angles>> temp = data;
	deque<pair<double,Angles>> minimums;

	do {
		double preDiff = 0;
		double preVal = 0;

		std::pair<double,Angles> preMin;

		minimums.clear();
		minimums = temp;
		temp.clear();

		for(auto& min : minimums) {
			double diff;
			if (min.first == minimums.front().first) {
				preVal = min.second.theta;
				preMin = min;
			} else {
				diff = min.second.theta - preVal;
				double error = fabs(diff / min.second.theta);
				if (preDiff < 0 && diff > 0 && error > acceptableDiff)
					temp.push_back(preMin);

				preVal = min.second.theta;
				preDiff = diff;
				preMin = min;
			}
		}
	} while (temp.size() > 5);

	return minimums;
}

string Experiment::findMode() {
	string mode = "NotDetected";
	bool isSConst = this->compare(this->processMaxima(maxSmallTheta));
	bool isBConst = this->compare(this->processMaxima(maxBigTheta));
	bool isSLagConst = this->compareLagAngle(maxSmallTheta);
	bool isBLagConst = this->compareLagAngle(maxBigTheta);
	bool isSTimeConst;
	bool isBTimeConst;

	this->acceptableDiff = 1e-8;

	for (int i = 0; i < 5; i++) {
		isBTimeConst = this->compareDiffTime(this->processMaxima(maxBigTheta));
		isSTimeConst = this->compareDiffTime(this->processMaxima(maxSmallTheta));
		if ( isBTimeConst && isSTimeConst )
			break;
		else
			this->acceptableDiff *= 10;
	}

	if (isSConst && isBConst && isSLagConst && isBLagConst) {
		mode = "P-mode";
	} else if (isSConst && isBConst && isSTimeConst && isBTimeConst) {
		mode = "Q-mode";
	}

	return mode;
}
bool Experiment::compareDiffTime(const deque<pair<double,Angles>> data) {
	double diffTime = 0;
	double preTime = 0;
	double preDiffTime = 0;
	unsigned int count = 0;

	for (auto max : data){
		if (max.first == data.front().first)
			preTime = max.first;
		else {
			diffTime = max.first - preTime;
			if (relativeToleranceCompare(diffTime, preDiffTime) && preDiffTime > 0) {
				count++;
			}
			preTime = max.first;
			preDiffTime = diffTime;
		}
	}

	if (debug) {
		std::cout << "<theta>: " << data.back().second.theta;
		std::cout << "; diffTheta: " << fabs(data.front().second.theta - data.back().second.theta);
		std::cout << "; <diffTime>: " << (data.back().first - data.front().first)/(data.size() - 1);
		std::cout << "; diffTime: " << diffTime;
		std::cout << "; countEqualTimeDiff: " << count << "; size: " << data.size() << std::endl;
	}

	return count == (data.size() - 2) && data.size() > 3;
}

void Experiment::clearExtremums() {
	maxSmallTheta.clear();
	maxBigTheta.clear();
	minSmallTheta.clear();
	minBigTheta.clear();
}

void Experiment::print() {
	using namespace std;
	output << setiosflags(ios::scientific) << setprecision(outputPrecision)
			<< (int)round(p->getTime() / period) << "  "
			<< p->getTime() << "\t" << p->getSmallTheta() << "\t" << p->getBigTheta() << "\t"
			<< p->getSmallPhi() << "\t" << p->getBigPhi() << "\t" << p->getEnergy() << "\t"
			<< getEnergyPerTime(maxSmallTheta) << "\t" << this->var;
	if (debug)
		output << "\t" << p->getDeltaSmallPhi() << "\t" << p->getDeltaSmallTheta() << "\t"
			  << p->getDeltaBigPhi() << "\t" << p->getDeltaBigTheta() << "\t"
			  << p->getDeltaEnergy() << std::endl;
	else output << endl;
}

void Experiment::printLast() {
	using namespace std;

}


void Experiment::printConsole() {
	using namespace std;
	cout << std::fixed << setprecision(8)
		<< (int)round(p->getTime() / period) << "  "
		<< p->getTime() << "  " << p->getSmallTheta() << " " << p->getBigTheta()  << "  "
		<< p->getSmallPhi() << "  " << p->getBigPhi()<< "  " << p->getEnergy() << "  "
		<< getEnergyPerTime(maxSmallTheta) << " " << p->getEnergy() / p->getTime() << " "
		<< this->var << endl;
}

void Experiment::printTask() {
	fstream outputTask;
    auto mode = this->findMode();

	auto resMaxTheta = calcCharacters(maxSmallTheta);
	auto resMinTheta = calcCharacters(minSmallTheta);
    auto firstDiffTheta = calcCharacters(processMaxima(maxSmallTheta));

    string fileName = fileDirectory + "mMoment_"+ outputFileName;
    outputTask.open(fileName.c_str(), std::fstream::app | std::fstream::out);

	outputTask << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
				<< p->getFrequency() << '\t' << p->getAmplitude() << '\t'
				<< resMaxTheta.energy << "\t"
				<< resMaxTheta.mTheta << "\t" << resMinTheta.mTheta << "\t"
				<< this->averageSmallTheta  << "\t" << resMaxTheta.periodTheta << "\t"
				<< firstDiffTheta.periodTheta << "\t"
				<< resMaxTheta.velocityPhi << "\t" << resMaxTheta.phasePhi << "\t"
				<< mode
				<< endl;
	outputTask.close();

	resMaxTheta = calcCharacters(maxBigTheta);
	resMinTheta = calcCharacters(minBigTheta);
    firstDiffTheta = calcCharacters(processMaxima(maxBigTheta));

	fileName = fileDirectory + "pBody_"+ outputFileName;
	outputTask.open(fileName.c_str(), std::fstream::app | std::fstream::out);

	outputTask << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
				<< p->getFrequency() << '\t' << p->getAmplitude() << '\t'
				<< resMaxTheta.energy << "\t"
				<< resMaxTheta.mTheta << "\t" << resMinTheta.mTheta << "\t"
				<< this->averageBigTheta  << "\t" << resMaxTheta.periodTheta << "\t"
				<< firstDiffTheta.periodTheta << "\t"
				<< resMaxTheta.velocityPhi << "\t" << resMaxTheta.phasePhi << "\t"
				<< mode
				<< endl;
	outputTask.close();

	if (debug) {
		auto maxS = processMaxima(maxSmallTheta);
		auto maxB = processMaxima(maxBigTheta);

		auto diffS = fabs(maxS.front().second.theta - maxS.back().second.theta) / maxS.back().second.theta;
		auto diffB = fabs(maxB.front().second.theta - maxB.back().second.theta) / maxB.back().second.theta;
		auto diffSP = fabs(maxS.front().second.phi - p->getFrequency() * maxS.front().first - maxS.back().second.phi  + p->getFrequency() * maxS.back().first) / (maxS.back().second.phi  - p->getFrequency() * maxS.back().first);
		auto diffBP = fabs(maxB.front().second.phi - p->getFrequency() * maxB.front().first - maxB.back().second.phi  + p->getFrequency() * maxB.back().first) / (maxB.back().second.phi  - p->getFrequency() * maxB.back().first);

		auto diffSTimeBefore = (maxSmallTheta.back().first - maxSmallTheta.front().first) / (maxSmallTheta.size() - 1);
		auto diffBTimeBefore = (maxBigTheta.back().first - maxBigTheta.front().first) / (maxBigTheta.size() - 1);
		auto diffSTime = (maxS.back().first - maxS.front().first) / (maxS.size() - 1);
		auto diffBTime = (maxB.back().first - maxB.front().first) / (maxB.size() - 1);


		fileName = fileDirectory + "mode_"+ outputFileName;
		outputTask.open(fileName.c_str(), std::fstream::app | std::fstream::out);

		outputTask << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
			<< diffS << '\t' << diffB << "\t" << diffSP << '\t' << diffBP << '\t'
			<< diffSTimeBefore << '\t' << diffSTime << "\t" << diffBTimeBefore << '\t' << diffBTime << '\t'
			<< mode << endl;

		outputTask.close();
	}
}

void Experiment::printError(){
	if (p->getError() != "") {
		outputDebug << "Error for run that located in " << outputFileName << " file" << std::endl;
		outputDebug << p->getError() << std::endl;
		p->resetError();
	}
}

void Experiment::printDebugInfo() {

	fstream outputResults;
	outputResults.open("debug.txt", std::fstream::app | std::fstream::out);

	auto resMaxSmallTheta = calcCharacters(maxSmallTheta);
	auto resMaxBigTheta = calcCharacters(maxBigTheta);


	std::cout << "=========== Before processing of maximas ===========" << std::endl;
	std::cout << "++++++++++++++++++++ smallTheta ++++++++++++++++++++" << std::endl;
	printCharacters(resMaxSmallTheta, "maxSmallTheta_" + this->outputFileName);
	std::cout << "++++++++++++++++++++ bigTheta ++++++++++++++++++++" << std::endl;
	printCharacters(resMaxBigTheta, "maxBigTheta_" + this->outputFileName);

	outputResults << "Before processing maxSmallTheta:" << std::endl;
	outputResults << "time, theta, phi" << std::endl;

    for(auto& max : maxSmallTheta) {
    	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
    				<< max.first << '\t' << max.second.theta << '\t' << max.second.phi << std::endl;
    }

    outputResults << "Before processing maxBigTheta:" << std::endl;
	outputResults << "time, theta, phi" << std::endl;

    for(auto& max : maxBigTheta) {
    	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
    				<< max.first << '\t' << max.second.theta << '\t' << max.second.phi << std::endl;
    }

    auto firstDiffSmallTheta = processMaxima(maxSmallTheta);
    auto firstDiffBigTheta = processMaxima(maxBigTheta);

	resMaxSmallTheta = calcCharacters(firstDiffSmallTheta);
	resMaxBigTheta = calcCharacters(firstDiffBigTheta);

    std::cout << "============ After processing of maximas ============" << std::endl;
    std::cout << "++++++++++++++++++++ smallTheta ++++++++++++++++++++" << std::endl;
    printCharacters(resMaxSmallTheta, "maxSmallTheta_" + this->outputFileName);
    std::cout << "++++++++++++++++++++ bigTheta ++++++++++++++++++++" << std::endl;
    printCharacters(resMaxBigTheta, "maxBigTheta_" + this->outputFileName);

    outputResults << "After proccessing of maxSmallTheta:" << std::endl;
	outputResults << "time, theta, phi" << std::endl;

    for(auto& max : firstDiffSmallTheta) {
    	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
    				<< max.first << '\t' << max.second.theta << '\t' << max.second.phi << std::endl;
    }

    outputResults << "After proccessing of maxBigTheta:" << std::endl;
	outputResults << "time, theta, phi" << std::endl;

    for(auto& max : firstDiffBigTheta) {
    	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
    				<< max.first << '\t' << max.second.theta << '\t' << max.second.phi << std::endl;
    }

	auto resMinSmallTheta = calcCharacters(minSmallTheta);
	auto resMinBigTheta = calcCharacters(minBigTheta);

	std::cout << "=========== Before processing of minimums ===========" << std::endl;
	std::cout << "++++++++++++++++++++ smallTheta ++++++++++++++++++++" << std::endl;
	printCharacters(resMinSmallTheta, "minSmallTheta_" + this->outputFileName);
	std::cout << "++++++++++++++++++++ bigTheta ++++++++++++++++++++" << std::endl;
	printCharacters(resMinBigTheta, "minBigTheta_" + this->outputFileName);

	outputResults << "Before processing minSmallTheta:" << std::endl;
	outputResults << "time, theta, phi" << std::endl;

    for(auto& min : minSmallTheta) {
    	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
    				<< min.first << '\t' << min.second.theta << '\t' << min.second.phi << std::endl;
    }

    outputResults << "Before processing minBigTheta:" << std::endl;
	outputResults << "time, theta, phi" << std::endl;

    for(auto& min : minBigTheta) {
    	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
    				<< min.first << '\t' << min.second.theta << '\t' << min.second.phi << std::endl;
    }

	firstDiffSmallTheta = processMinimum(minSmallTheta);
	firstDiffBigTheta = processMinimum(minBigTheta);

	resMinSmallTheta = calcCharacters(firstDiffSmallTheta);
	resMinBigTheta = calcCharacters(firstDiffBigTheta);

	std::cout << "=========== After processing of minimums ===========" << std::endl;
	std::cout << "++++++++++++++++++++ smallTheta ++++++++++++++++++++" << std::endl;
	printCharacters(resMinSmallTheta, "minSmallTheta_" + this->outputFileName);
	std::cout << "++++++++++++++++++++ bigTheta ++++++++++++++++++++" << std::endl;
	printCharacters(resMinBigTheta, "minBigTheta_" + this->outputFileName);

	outputResults << "After processing minSmallTheta:" << std::endl;
	outputResults << "time, theta, phi" << std::endl;

    for(auto& min : firstDiffSmallTheta) {
    	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
    				<< min.first << '\t' << min.second.theta << '\t' << min.second.phi << std::endl;
    }

    outputResults << "After processing minBigTheta:" << std::endl;
	outputResults << "time, theta, phi" << std::endl;

    for(auto& min : firstDiffBigTheta) {
    	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision)
    				<< min.first << '\t' << min.second.theta << '\t' << min.second.phi << std::endl;
    }

	outputResults.close();
}

void Experiment::printCharacters(const CharactersOfMotion &data, string fileName) {
	auto results = data;

	fstream outputResults;
	outputResults.open(fileName, std::fstream::app | std::fstream::out);

	outputResults << setiosflags(std::ios::scientific) << std::setprecision(outputPrecision);
	outputResults << p->getFrequency() << "\t";
	outputResults << results.mTheta << "\t" << results.periodTheta << "\t";
	outputResults << results.phasePhi << "\t"  << results.velocityPhi << "\t" << results.energy << endl;
	outputResults.close();

	if (debug) {
		std::cout << std::fixed << std::setprecision(8);
		std::cout << "theta: " << results.mTheta;
		std::cout << "\tvelocityPhi: " << results.velocityPhi << std::endl;
		std::cout << "periodTheta: " << results.periodTheta;
		std::cout << "\tphasePhi: " << results.phasePhi << std::endl;
		std::cout << "energy/time: " << results.energy << std::endl;
	}
}

void Experiment::printHeader() {
	using namespace std;
	using namespace Exprm;

	cout << "period  " << "time  " << "smallTheta  " << "bigTheta  " << "smallPhi  "
						<< "bigPhi  " << "energy  " << "energy/time  " << this->getVariable();
	if (debug)
		cout << "  " << "dSmallPhi  " << "dSmallTheta  " << "dBigPhi  "
			 	 	 << "dBigTheta  " << "dEnergy  "<< endl;
	else cout << endl;

	output << "period\t" << "time\t" << "smallTheta\t" << "bigTheta\t" << "smallPhi\t"
				  	  << "bigPhi\t" << "energy\t" << "energy/time\t" << this->getVariable();
	if (debug)
		output << "\t " << "dSmallPhi\t" << "dSmallTheta\t" << "dBigPhi\t"
						<< "dBigTheta\t" << "dEnergy\t"<< endl;
	else output << endl;

	fstream outputTask;

    string fileName = fileDirectory + "mMoment_"+ outputFileName;
    outputTask.open(fileName.c_str(), std::fstream::app | std::fstream::out);

	outputTask << "frequency\t" << "amplitude\t" << "energy/time\t"
				<< "maxSmallTheta\t" << "minSmallTheta\t"
				<< "averageSmallTheta\t" << "periodTheta(before)\t"
				<< "periodTheta(after)\t"
				<< "velocitySmallPhi\t" << "phaseSmallPhi\t"
				<< "mode"
				<< endl;
	outputTask.close();

    fileName = fileDirectory + "pBody_"+ outputFileName;
    outputTask.open(fileName.c_str(), std::fstream::app | std::fstream::out);

	outputTask << "frequency\t" << "amplitude\t" << "energy/time\t"
				<< "maxBigTheta\t" << "minBigTheta\t"
				<< "averageBigTheta\t" << "periodTheta(before)\t"
				<< "periodTheta(after)\t"
				<< "velocityBigPhi\t" << "phaseBigPhi\t"
				<< "mode"
				<< endl;
	outputTask.close();
}

void Experiment::relaxation(int nPeriods) {

	std::cout << "===== Skipping " << skipPeriods << " periods =====" << std::endl;
	output << "===== Skipping " << skipPeriods << " periods =====" << std::endl;

	cutEnergy = 0.0;
	cutTime = 0.0;

	if (nPeriods > 0) {
		nRuns = nPeriods * period / p->getTimeStep();
		run(1);
		nRuns = period / p->getTimeStep();
	}

	cutEnergy = p->getEnergy();
	cutTime = p->getTime();

	output << "Cut energy:	" << cutEnergy << std::endl;
	output << "Cut time:	" << cutTime << std::endl;

	std::cout << "===== Relaxation has finished =====" << std::endl;
	output << "===== Relaxation has finished =====" << std::endl;
}

double Experiment::getEnergyPerTime() {
	return (p->getEnergy() - cutEnergy) / (p->getTime() - cutTime);
}

double Experiment::getEnergyPerTime(const deque<pair<double,Angles>> &data) {
	if (data.size() > 0)
		return (data.back().second.energy - data.front().second.energy) / (data.back().first - data.front().first);
	else
		return p->getEnergy() / p->getTime();
}

void Experiment::hysteresis() {
	for (long long int i = 0; i < numberOfPeriods; ++i) {
		for (long long int j = 0; j < nRuns; ++j) {
			p->RungeKutta4Step();
			p->applyDeltas();

			if (p->isLinearPolarization() == true)
				output << "hz\t" << "Mz" << std::endl;
			else
				output << "hx+hy\t" << "Mx\t" << "My" << std::endl;

			if (i * printStep == 0) {
				if (p->isLinearPolarization() == true) {
					std::cout << p->getHz() << '\t' << p->getMz() << std::endl;
					output << p->getHz() << '\t' << p->getMz() << std::endl;
				}
				else {
					std::cout << p->getHx() + p->getHy() << '\t' << p->getMx() << '\t' << p->getMy() << std::endl;
					output << p->getHx() + p->getHy() << '\t' << p->getMx() << '\t' << p->getMy() << std::endl;
				}
			}
		}
	}
}

void Experiment::setVar() {
	switch (this->itsVar) {
		case Exprm::Variable::Amplitude:
			p->setAmplitude(this->var);
			break;
		case Exprm::Variable::Frequency:
			p->setFrequency(this->var);
			break;
		case Exprm::Variable::NotSet:
			this->var = 0;
			this->varMax = 0;
			this->varStep = 1;
			break;
		default:
			break;
	}
}

std::string Experiment::getVariable() {
	using namespace Exprm;
	std::string var = "";
	switch (itsVar) {
		case Variable::Amplitude:
			var = "Amplitude";
			break;
		case Variable::Frequency:
			var = "Frequency";
			break;
		case Variable::NotSet:
			var = "NotSet";
			break;
		default:
			var = "NotSet";
			break;
	}
	return var;
}

std::string Experiment::getTask() {
	using namespace Exprm;
	std::string task = "";
	switch (itsTask) {
		case Task::SimpleRun:
			task = "SimpleRun";
			break;
		case Task::PowerLoss:
			task = "PowerLoss";
			break;
		default:
			task = "Error";
			break;
	}
	return task;
}

std::string Experiment::getNumMethod() {
	using namespace Exprm;
	std::string method = "";
	switch (itsMethod) {
		case NumericalMethod::Euler:
			method = "Euler";
			break;
		case NumericalMethod::RK2:
			method = "RK2";
			break;
		case NumericalMethod::RK4:
			method = "RK4";
			break;
		default:
			method = "ERROR";
			break;
	}
	return method;
}

void Experiment::saveConf() {
	using namespace Exprm;

	output << "# Task\t= " << getTask() << std::endl;
	output << "# numericalMethod\t= " << getNumMethod() << std::endl;
	output << "# debug\t= " << debug << std::endl;
	output << "# intelegentStop\t= " << intelegentStop << std::endl;
	output << "# Variable\t= " << getVariable() << std::endl;

	output << "# numberOfPeriods\t= " << numberOfPeriods << std::endl;
	output << "# maxNumberOfPeriods\t= " << maxNumberOfPeriods << std::endl;
	output << "# numberOfExtremums\t= " << numberOfExtremums << std::endl;
	output << "# skipPeriods\t= " << skipPeriods << std::endl;
	output << "# printStep\t= " << printStep << std::endl;
	output << "# printLastNperiods\t= " << printLastNperiods << std::endl;
	output << "# timeStep\t= " << p->getTimeStep() << std::endl;
	output << "# outputPrecision\t= " << outputPrecision << std::endl;
	output << "# acceptableDiff\t= " << acceptableDiff << std::endl;

	output << "# alpha\t= " << p->getAlpha() << std::endl;
	output << "# amplitude\t= " << p->getAmplitude() << std::endl;
	output << "# Hz\t= " << p->getZAmplitude() << std::endl;
	output << "# frequency\t= " << p->getFrequency() << std::endl;
	output << "# magnetization\t= " << p->getMagnetization() << std::endl;
	output << "# viscosity\t= " << p->getViscosity() << std::endl;

	output << "# varMin\t= " << var << std::endl;
	output << "# varMax\t= " << varMax << std::endl;
	output << "# varStep\t= " << varStep << std::endl;

	output << "# smallTheta\t= " << p->getSmallTheta() << std::endl;
	output << "# smallPhi\t= " << p->getSmallPhi() << std::endl;
	output << "# bigTheta\t= " << p->getBigTheta ()<< std::endl;
	output << "# bigPhi\t= " << p->getBigPhi() << std::endl;

	output << "# thermalBath\t= " << p->isThermalBath() << std::endl;
	output << "# linearPolarization\t= " << p->isLinearPolarization() << std::endl;
	output << "# hySign\t= " << p->getRo() << std::endl;

	std::cout << "===== Skipping " << skipPeriods << " periods =====" << std::endl;
	output << "===== Skipping " << skipPeriods << " periods =====" << std::endl;

}
