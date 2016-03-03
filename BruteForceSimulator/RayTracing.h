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

#define Upwards   1
#define Downwards 0

#define Polarised false

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
	
	void SamplerNormalDirection(int lightDirectionStatus,int intersectionLayerNum,Vector & sampledLocalNormal);
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
	void ChooseReflectionOrRefraction(
		const double    localtheta,
		const double    n_i,
		const double    k_i,
		const double    n_t,
		const double    k_t,
		const int       lightDirectionStatus,
			  bool    & doReflection
			  );
	void DoReflectionOrNot(
		const double    localTheta,
		const double    n_i,
		const double    k_i,
		const double    n_t,
		const double    k_t,
		const int       lightDirectionStatus,
		const int       intersectionLayer,
			  bool    & doReflection
		);
	void setFresnelMueller(Mueller &mu,double A,double B,double C,double S);
	void TraceRayUpwards();
	void TraceRayDownwards();
	void TraceRay();
	
	int GetIntersetcionLayer(int currentL,int lightDirectionStatus);
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

