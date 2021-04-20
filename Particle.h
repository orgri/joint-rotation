/*
 * Particle.h
 *
 *  Created on: Sep 27, 2017
 *      Author: cuda
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <random>

class Particle {

private:
	const double gamma = 1.76e+7;
	double magnetization;
	double viscosity;
	double alpha;
	double beta;
	double amplitude;
	double omega;
	double Hz;		//z-amplitude

	double alpha1;
	double beta1;

	double dt;
	double time;

	double tau1;
	double tau2;

	double smallTheta;
	double smallPhi;
	double bigTheta;
	double bigPhi;
	double energy;

	double fSmallTheta;
	double fSmallPhi;
	double fBigTheta;
	double fBigPhi;

	double dSmallTheta;
	double dSmallPhi;
	double dBigTheta;
	double dBigPhi;
	double dEnergy;


	double hx (double t);
	double hy (double t);
	double hz (double t);

	int ro;		//sign before hy

	void calcFunctions(const double &smallTheta, const double &smallPhi,const double &bigTheta,const double &bigPhi, double t);

	bool thermalBath;
	bool linearPolarization;

	std::string error;

public:
	Particle();
	Particle(const double &smallTheta, const double &smallPhi,const double &bigTheta,const double &bigPhi);
	virtual ~Particle();

	void RungeKutta2Step();
	void RungeKutta4Step();
	void EulerStep();
	void calcEnergy();
	void applyDeltas();
	void resetError() { error = ""; }

	void setAlpha(double alpha) {
		this->alpha = alpha;
		this->beta = alpha * magnetization / (6.0 * gamma * viscosity);
		this->beta1 = 1.0 + beta;
		this->alpha1 = alpha / beta1;
		this->tau1 = 1.0 / (1.0 + alpha1 * alpha1);
		this->tau2 = beta * beta1 / alpha;
	}

	double getAlpha() { return alpha; }
	double getBeta() { return beta; }
	double getTimeStep() { return dt; }
	double getTime() { return time; }
	double getAmplitude() { return amplitude; }
	double getZAmplitude() { return Hz; }
	double getMagnetization() { return magnetization; }
	double getFrequency() { return omega; }
	double getViscosity() { return viscosity; }

	double getSmallTheta() { return smallTheta; }
	double getSmallPhi() { return smallPhi; }
	double getBigTheta() { return bigTheta; }
	double getBigPhi() { return bigPhi; }
	double getEnergy() { return energy; }

	double getDeltaSmallTheta() { return dSmallTheta; }
	double getDeltaSmallPhi() { return dSmallPhi; }
	double getDeltaBigTheta() { return dBigTheta; }
	double getDeltaBigPhi() { return dBigPhi; }
	double getDeltaEnergy() { return dEnergy; }

	double getHx() { return hx(time); }
	double getHy() { return hy(time); }
	double getHz() { return hz(time); }

	int getRo() { return ro; }

	double getMx() { return sin(smallTheta) * cos(smallPhi); }
	double getMy() { return sin(smallTheta) * sin(smallPhi); }
	double getMz() { return cos(smallTheta); }

	std::string getError() { return error; }

	bool isLinearPolarization() { return linearPolarization; }
	bool isThermalBath() { return thermalBath; }

	void setTimeStep(double dt) { this->dt = dt; }
	void setAmplitude(double h) { this->amplitude = h; }
	void setZAmplitude(double hz) { this->Hz = hz; }
	void setMagnetization(double magnetization) { this->magnetization = magnetization; }
	void setFrequency(double omega) { this->omega = omega; }
	void setLinearPolarization(bool linearPolarization) { this->linearPolarization = linearPolarization; }
	void setRo(int ro) { this->ro = ro; }
	void setThermalBath(bool thermalBath) { this->thermalBath = thermalBath; }
	void setViscosity(double viscosity) { this->viscosity = viscosity; }
	void setInitialState(const double &smallTheta, const double &smallPhi,const double &bigTheta,const double &bigPhi);
	void resetEnergy() { this->energy = 0.0; }
	void resetTime() {this->time = 0.0; }

};

#endif /* PARTICLE_H_ */
