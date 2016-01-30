#pragma once
#include "Polarization.h"
#include "Ray.h"
#include "Fresnel.h"
#include "AnisotropicSampler.h"
#include <vector>
#include <iostream>
#define PrintTraceInfo false
//#define PrintTraceInfo true
#define DebugMode false
#define TestNormal Normalize(Vector(1,-1,19))

using namespace std;

class RayTracing
{
private:
	vector<float> N,K,Ex,Ey;
	Vector surfaceNormal;
	Vector currentIncomingDirection;
	Vector currentOutgoingDirection;
	
	StokesVector currentStokes;
	int currentLayer;
	int NumOfLayer;
	unsigned int BounceLimit;
	unsigned int bounceNum;
	bool multiRelflectionHappen;
	void PrintStatus();
	

	void getRefractionAttenuation(
		const Vector    incomingDirecion,
		const Vector    localNormal,
		const double    n_i,
		const double    k_i,
		const double    n_t,
		const double    k_t,
		Vector  & outgoingDirecion,
		Mueller & mu
		);
	void getReflectionAttenuation(
		const Vector    incomingDirecion,
		const Vector    localNormal,
		const double    n_i,
		const double    k_i,
		const double    n_t,
		const double    k_t,
		Vector  & outgoingDirecion,
		Mueller & mu
		);
	void setFresnelMueller(Mueller &mu,double A,double B,double C,double S);
	void TraceRayUpwards();
	void TraceRayDownwards();
	Vector SampleAnisotropicNormal(int layer);
public:
	RayTracing(void);
	~RayTracing(void);
	void InitializeAndCompute(
		const Vector IncidenceDirection,
		const StokesVector   incidentceStokes,
		const vector<float>  NN,
		const vector<float>  KK,
		const vector<float>  EX,
		const vector<float>  EY,
		Vector       & outgoingDirection,
		StokesVector & sv
		);

	
};

