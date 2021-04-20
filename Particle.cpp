/*
 * Particle.cpp
 *
 *  Created on: Sep 27, 2017
 *      Author: cuda
 */

#include "Particle.h"

Particle::Particle() {
	// TODO Auto-generated constructor stub
	this->alpha = 1.0;
	this->time = 0.0;
	this->dt = 0.001;
	this->omega = 1.0;
	this->magnetization = 338.0;
	this->amplitude = 1.0;
	this->Hz = 0.0;

	this->beta = alpha * magnetization / (6.0 * gamma * viscosity);
	this->beta1 = 1.0 + beta;
	this->alpha1 = alpha / beta1;
	this->tau1 = 1.0 / (1.0 + alpha1 * alpha1);
	this->tau2 = beta * beta1 / alpha;

	this->thermalBath = false;
	this->linearPolarization = false;
	this->ro = 1;

	this->energy = 0.0;
	this->smallTheta = 0.1;
	this->smallPhi = 0;
	this->bigTheta = 0.1;
	this->bigPhi = 0;

	this->error = "";
}

Particle::Particle(const double &smallTheta, const double &smallPhi,const double &bigTheta,const double &bigPhi) :Particle(){
	this->smallTheta = smallTheta;
	this->smallPhi = smallPhi;
	this->bigTheta = bigTheta;
	this->bigPhi = bigPhi;
}

void Particle::setInitialState(const double &smallTheta, const double &smallPhi,const double &bigTheta,const double &bigPhi) {
	this->smallTheta = smallTheta;
	this->smallPhi = smallPhi;
	this->bigTheta = bigTheta;
	this->bigPhi = bigPhi;
	this->energy = 0.0;
	this->time = 0.0;
}

Particle::~Particle() {
	// TODO Auto-generated destructor stub
}

double Particle::hx(double t) {
	if (this->linearPolarization == true)
		return 0;
	else return amplitude * cos (omega * t);
}

double Particle::hy(double t) {
	if (this->linearPolarization == true)
		return 0;
	else return ro * amplitude * sin (omega * t);
}

double Particle::hz(double t) {
	if (this->linearPolarization == true)
		return amplitude * cos (omega * t);
	else return Hz;
}

void Particle::calcFunctions(const double& smallTheta, const double& smallPhi,
		const double& bigTheta, const double& bigPhi, double t) {

	double cosSmallTheta = cos(smallTheta);
	double cosBigTheta = cos(bigTheta);
	double cosSmallPhi = cos(smallPhi);
	double cosBigPhi = cos(bigPhi);

	double sinSmallTheta = sin(smallTheta);
	double sinBigTheta = sin(bigTheta);
	double sinSmallPhi = sin(smallPhi);
	double sinBigPhi = sin(bigPhi);

/*	double cosphi_PHI = cos(smallPhi - bigPhi);
	double sinphi_PHI = sin(smallPhi - bigPhi);

	double F = cosSmallTheta * cosBigTheta + cosphi_PHI * sinBigTheta * sinSmallTheta;
	double f1 =  - ( F * sinphi_PHI * sinBigTheta + beta1 * ( hx(t) * sinSmallPhi - hy(t) * cosSmallPhi) );
	double f2 = cosSmallTheta * ( F * cosphi_PHI * sinBigTheta + beta1 * (hx(t) * cosSmallPhi + hy(t) * sinSmallPhi) )
				- sinSmallTheta * ( F * cosBigTheta + beta1 * hz(t) );*/

	double F = cosSmallTheta * cosBigTheta + cos(smallPhi - bigPhi) * sinBigTheta * sinSmallTheta;

	double c1 = F * cos(smallPhi - bigPhi) * sinBigTheta;
	double c2 = F * sin(smallPhi - bigPhi) * sinBigTheta;

	double h1 = hx(t) * cosSmallPhi + hy(t) * sinSmallPhi;
	double h2 = hx(t) * sinSmallPhi - hy(t) * cosSmallPhi;

	double f1 = -(c2 + beta1 * h2);
	double f2 = (c1 + beta1 * h1) * cosSmallTheta - (F * cosBigTheta + beta1 * hz(t)) * sinSmallTheta;

	fSmallTheta = tau1 * (f1  + alpha1 * f2);
	fSmallPhi = tau1 * (alpha1 * f1 - f2) / sinSmallTheta;

	double wx = (fSmallTheta * cosSmallTheta * cosSmallPhi - fSmallPhi * sinSmallTheta * sinSmallPhi) / beta1
				+  hz(t) * sinSmallTheta * sinSmallPhi - hy(t) * cosSmallTheta ;
	double wy = (fSmallTheta * cosSmallTheta * sinSmallPhi + fSmallPhi * sinSmallTheta * cosSmallPhi) / beta1
				+ hx(t) * cosSmallTheta - hz(t) * sinSmallTheta * cosSmallPhi ;
	double wz = -(fSmallTheta / beta1 + h2) * sinSmallTheta;

	fBigTheta = tau2 * (wy * cosBigPhi - wx * sinBigPhi);
	fBigPhi = tau2 * (wz - cosBigTheta * (wx * cosBigPhi + wy * sinBigPhi) / sinBigTheta);

}

void Particle::EulerStep() {
	calcFunctions(smallTheta, smallPhi, bigTheta, bigPhi, time);

	dSmallTheta = dt * fSmallTheta;
	dSmallPhi = dt * fSmallPhi;
	dBigTheta = dt * fBigTheta;
	dBigPhi = dt * fBigPhi;
}

void Particle::RungeKutta2Step() {
	calcFunctions(smallTheta, smallPhi, bigTheta, bigPhi, time);

	double kSmallTheta =  fSmallTheta;
	double kSmallPhi =  fSmallPhi;
	double kBigTheta =  fBigTheta;
	double kBigPhi =  fBigPhi;

	double dkSmallTheta = dt * fSmallTheta  ;
	double dkSmallPhi = dt * fSmallPhi;
	double dkBigTheta = dt * fBigTheta  ;
	double dkBigPhi = dt * fBigPhi;

	calcFunctions(smallTheta + dkSmallTheta, smallPhi + dkSmallPhi, bigTheta + dkBigTheta, bigPhi + dkBigPhi, time + dt);

	dSmallTheta = 0.5 * dt * ( fSmallTheta +  kSmallTheta ) ;
	dSmallPhi = 0.5 * dt * ( fSmallPhi + kSmallPhi );
	dBigTheta = 0.5 * dt * ( fBigTheta +  kBigTheta ) ;
	dBigPhi = 0.5 * dt * ( fBigPhi + kBigPhi );

}

void Particle::RungeKutta4Step() {
	calcFunctions(smallTheta, smallPhi, bigTheta, bigPhi, time);

	double k1SmallTheta = dt * fSmallTheta  ;
	double k1SmallPhi = dt * fSmallPhi;
	double k1BigTheta = dt * fBigTheta  ;
	double k1BigPhi = dt * fBigPhi;

	calcFunctions(smallTheta + 0.5 * k1SmallTheta, smallPhi + 0.5 * k1SmallPhi, bigTheta + 0.5 * k1BigTheta, bigPhi + 0.5 * k1BigPhi, time + 0.5 * dt);

	double k2SmallTheta = dt * fSmallTheta;
	double k2SmallPhi = dt * fSmallPhi;
	double k2BigTheta = dt * fBigTheta;
	double k2BigPhi = dt * fBigPhi;

	calcFunctions(smallTheta + 0.5 * k2SmallTheta, smallPhi + 0.5 * k2SmallPhi, bigTheta + 0.5 * k2BigTheta, bigPhi + 0.5 * k2BigPhi, time + 0.5 * dt);

	double k3SmallTheta = dt * fSmallTheta;
	double k3SmallPhi = dt * fSmallPhi;
	double k3BigTheta = dt * fBigTheta;
	double k3BigPhi = dt * fBigPhi;

	calcFunctions(smallTheta + k3SmallTheta, smallPhi + k3SmallPhi, bigTheta + k3BigTheta, bigPhi + k3BigPhi, time + dt);

	double k4SmallTheta = dt * fSmallTheta;
	double k4SmallPhi = dt * fSmallPhi;
	double k4BigTheta = dt * fBigTheta;
	double k4BigPhi = dt * fBigPhi;

	dSmallTheta = (k1SmallTheta + 2.0 * k2SmallTheta + 2.0 * k3SmallTheta + k4SmallTheta) / 6.0;
	dSmallPhi = (k1SmallPhi + 2.0 * k2SmallPhi + 2.0 * k3SmallPhi + k4SmallPhi) / 6.0;
	dBigTheta = (k1BigTheta + 2.0 * k2BigTheta + 2.0 * k3BigTheta + k4BigTheta) / 6.0;
	dBigPhi = (k1BigPhi + 2.0 * k2BigPhi + 2.0 * k3BigPhi + k4BigPhi) / 6.0;

}

void Particle::calcEnergy() {

	double cosSmallTheta = cos(smallTheta);
	double cosBigTheta = cos(bigTheta);
	double cosSmallPhi = cos(smallPhi);

	double sinSmallTheta = sin(smallTheta);
	double sinBigTheta = sin(bigTheta);
	double sinSmallPhi = sin(smallPhi);

	double F = cosSmallTheta * cosBigTheta + cos(smallPhi - bigPhi) * sinBigTheta * sinSmallTheta;

	double h1 = hx(time) * cosSmallPhi + hy(time) * sinSmallPhi;
	double h2 = hx(time) * sinSmallPhi - hy(time) * cosSmallPhi;
	double h3 = hz(time);

	double c1 = F * cos(smallPhi - bigPhi) * sinBigTheta;
	double c2 = F * sin(smallPhi - bigPhi) * sinBigTheta;
	double c3 = F * cosBigTheta;

	dEnergy = cosSmallTheta * (h1 + c1) * dSmallTheta
			- sinSmallTheta * (h2 + c2) * dSmallPhi
			- sinSmallTheta * (h3 + c3) * dSmallTheta;
}

void Particle::applyDeltas() {
	smallTheta += dSmallTheta;
	smallPhi += dSmallPhi;
	bigTheta += dBigTheta;
	bigPhi += dBigPhi;
	energy += dEnergy;
	time += dt;
	
	if (smallTheta < 0 || smallTheta > M_PI) {
		error = "smallTheta is out of range: time=" + std::to_string(time) + " smallTheta=" + std::to_string(smallTheta) + " dSmallTheta=" + std::to_string(dSmallTheta);
		smallTheta -= 2.0 * dSmallTheta;
	}

	
	if (bigTheta < 0 || bigTheta > M_PI) {
		error = "bigTheta is out of range: time=" + std::to_string(time) + " bigTheta=" + std::to_string(bigTheta) + " dBigTheta=" + std::to_string(dBigTheta);
		bigTheta -= 2.0 * dBigTheta;
	}
/*	if (smallTheta < 0) {
		error = "smallTheta is out of range: time=" + std::to_string(time) + " smallTheta=" + std::to_string(smallTheta) + " dSmallTheta=" + std::to_string(dSmallTheta);
		smallTheta = 1e-6;
	}

	if (smallTheta > M_PI) {
		error = "smallTheta is out of range: time=" + std::to_string(time) + " smallTheta=" + std::to_string(smallTheta) + " dSmallTheta=" + std::to_string(dSmallTheta);
		smallTheta = M_PI - 1e-6;
	}

	if (bigTheta < 0) {
		error = "bigTheta is out of range: time=" + std::to_string(time) + " bigTheta=" + std::to_string(bigTheta) + " dBigTheta=" + std::to_string(dBigTheta);
		bigTheta = 1e-6;
	}

	if (bigTheta > M_PI) {
		error = "bigTheta is out of range: time=" + std::to_string(time) + " bigTheta=" + std::to_string(bigTheta) + " dBigTheta=" + std::to_string(dBigTheta);
		bigTheta = M_PI - 1e-6;
	}*/
}
/*
void Particle::checkExtremum() {
	static int prevdSmallTheta;
	static int prevdSmallPhi;
	static int prevdBigTheta;
	static int prevdBigPhi;
	
	if (prevdSmallTheta > 0 && dSmallTheta < 0)
		maxSmallTheta.append(smallTheta);
	if (prevdSmallTheta < 0 && dSmallTheta > 0)
		minSmallTheta.append(smallTheta);
	
	if (prevdSmallPhi > 0 && dSmallPhi < 0)
		maxSmallPhi.append(smallPhi);
	if (prevdSmallPhi < 0 && dSmallPhi > 0)
		minSmallPhi.append(smallPhi);
	
	if (prevdBigTheta > 0 && dBigTheta < 0)
		maxBigTheta.append(bigTheta);
	if (prevdBigTheta < 0 && dBigTheta > 0)
		minBigTheta.append(bigTheta);
	
	if (prevdBigPhi > 0 && dBigPhi < 0)
		maxBigPhi.append(bigPhi);
	if (prevdBigPhi < 0 && dBigPhi > 0)
		minBigPhi.append(bigPhi);	
	
	prevdSmallTheta = dSmallTheta;
	prevdSmallPhi = dSmallPhi;
	prevdBigTheta = dBigTheta;
	prevdBigPhi = dBigPhi;
}
*/
