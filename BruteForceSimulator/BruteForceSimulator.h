//  Created by Chi Wang on 04/12/15.
//  Copyright @ 2015 Chi Wang. All rights reserved.
//

#pragma once
#include <vector>
#include "Polarization.h"
#include "Ray.h"

#include "RayTracing.h"
//#include "windows.h"
#include "omp.h"

#include <fstream>




class BruteForceSimulator
{
public:
	vector<float> N,K,Ex,Ey;            //IOR of each layer ,  air as first layer
	StokesVector *reflection_data[91][360];
	Vector surfaceNormal;
	int NumOfLayer;
	int NumOfIteration;
	int effectiveSampleNum;
	Vector IncidenceDirection;
	StokesVector incidenceStokes;
	
	vector<Ray> outgoingRay;

	BruteForceSimulator(void);
	~BruteForceSimulator(void);

	void CollectData(Vector direction,StokesVector sv);
	void SaveToFile();
	
	void Calculate(void);
	void InitializeAndCompute();
	void AddAnisotropicLayer(float n, float k,float ex,float ey);
	void ReadAnisotropicSettingFile();
	void SetIncidence(float theta, float phi, StokesVector stokes);
	

	


};

