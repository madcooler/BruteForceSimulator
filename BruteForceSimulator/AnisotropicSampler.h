//  Created by Chi Wang on 04/12/15.
//  Copyright @ 2015 Chi Wang. All rights reserved.
//

#pragma once
#include "Vector.h"

inline double getMaxValue(double a, double b)
{
	if (a>b)
	{
		return a;
	}
	else return b;

}
inline double RandomNum()
{
	return rand()/(double)RAND_MAX;
}

inline double GaussianNoise(double variance,double mean)
{ 
	
	double Temp;
	double u1,u2;
	u1=RandomNum();
	u2=RandomNum();
	Temp=sqrt(-2*(log(u1)))*sin(2*M_PI*(u2));
	Temp=variance*Temp+mean;	
	return Temp;
}
inline void sampleFirstQuadrant(double ex,double ey,double u1, double u2, double *phi, double *costheta) 
{
		if (ex == ey)
			*phi = M_PI * u1 * 0.5f;
		else
			*phi = atanf(sqrtf((ex+1.f) / (ey+1.f)) *
                     tanf(M_PI * u1 * 0.5f));
		double cosphi = cosf(*phi), sinphi = sinf(*phi);
		*costheta = powf(u2, 1.f/(ex * cosphi * cosphi +
                              ey * sinphi * sinphi + 1));
}

inline void SampleAnisotropic(double ex,double ey, Vector &wh)
{
	double u1;
	double u2;
	u1 = RandomNum();
	u2 = RandomNum();

	double phi, costheta;

    if (u1 < .25f) 
	{
		sampleFirstQuadrant(ex,ey,4.f * u1, u2, &phi, &costheta);
    } 
	else if (u1 < .5f) 
	{
        u1 = 4.f * (.5f - u1);
        sampleFirstQuadrant(ex,ey,u1, u2, &phi, &costheta);
        phi = M_PI - phi;
    } 
	else if (u1 < .75f) 
	{
        u1 = 4.f * (u1 - .5f);
        sampleFirstQuadrant(ex,ey,u1, u2, &phi, &costheta);
        phi += M_PI;
    } 
	else {
        u1 = 4.f * (1.f - u1);
        sampleFirstQuadrant(ex,ey,u1, u2, &phi, &costheta);
        phi = 2.f * M_PI - phi;
    }
    double sintheta = sqrtf(getMaxValue(0.f, 1.f - costheta*costheta));
    wh = SphericalDirection(sintheta, costheta, phi);
}




