#include "RayTracing.h"


RayTracing::RayTracing(void)
{
	currentIncomingDirection = Vector(0,0,0);
	currentOutgoingDirection = Vector(0,0,0);
	surfaceNormal			 = Vector(0,0,1);
	currentStokes            = StokesVector(0,0,0,0);
	
	multiRelflectionHappen	 = false;
	currentLayer = 0;
	bounceNum    = 0;
}


RayTracing::~RayTracing(void)
{
}

void RayTracing::InitializeAndCompute(
			const Vector IncidenceDirection,
			const StokesVector   incidentceStokes,
			const vector<float>  NN,
			const vector<float>  KK,
			const vector<float>  EX,
			const vector<float>  EY,
				Vector       & outgoingDirection,
				StokesVector & sv
			)
{
	if(PrintTraceInfo == true)
	{
		printf("..........start tracing.......... \n");
		printf("Current layer : ");
	}
	multiRelflectionHappen = false;
	currentLayer = 0;
	bounceNum    = 0;

	currentStokes = incidentceStokes;
	currentIncomingDirection = IncidenceDirection;

	N = NN;
	K = KK;
	Ex = EX;
	Ey = EY;

	NumOfLayer = Ex.size() ;
	BounceLimit				 = NumOfLayer * 4;

	TraceRayDownwards();

	if(PrintTraceInfo == true)
		printf("\n");

	if (bounceNum < BounceLimit && multiRelflectionHappen == false)
	{
		outgoingDirection = currentOutgoingDirection;
		sv = currentStokes;
	}
	else
	{
		if(PrintTraceInfo == true)
		{
			printf("failed to sample this ray, which is due to ");
			if (bounceNum == BounceLimit)
				printf("Bounce Number meets upper limit : %d \n",bounceNum);
			if(multiRelflectionHappen == true)
				printf("multireflection \n");
		}
	}
}

void RayTracing::TraceRayDownwards()
{
	if(PrintTraceInfo == true)
		printf(" %d ", currentLayer);

	if(bounceNum == BounceLimit)
	{
		//printf("should stop %d ", bounceNum);
		return;
	}

	if (bounceNum < BounceLimit)
		bounceNum++;

	int layerNumOfTracedRay = currentLayer;

	double n_i,n_t,k_i,k_t;
	if (layerNumOfTracedRay < NumOfLayer)
	{
		n_i = N[layerNumOfTracedRay];
		k_i = K[layerNumOfTracedRay];
		n_t = N[layerNumOfTracedRay+1];
		k_t = K[layerNumOfTracedRay+1];
	}
	else
		cout<<" layerNumOfTracedRay+1 > NumOfLayer "<<endl;
	

	Vector SampledLocalNormal;
	SampledLocalNormal = SampleAnisotropicNormal(layerNumOfTracedRay+1);

	double localTheta=acos(abs(Dot(
		currentIncomingDirection,
		SampledLocalNormal)));


	double Reflectance=0,Transmittance=0;
	double Ts,Tp,TsTp,Tphase,Tfactor;
	double rs,rp,phaseS,phaseP,Fs,Fp;

	Reflectance = Fresnel_Reflection( localTheta, 
		n_i,
		k_i,
		n_t,
		k_t,
		rs,
		rp,
		phaseS,
		phaseP,
		Fs,
		Fp);
	if (k_t==0 && layerNumOfTracedRay < NumOfLayer - 1)  // next layer is not metallic or opaque AND there is not transimission at the last layer
	{
		Transmittance = Fresnel_Refraction(
			localTheta,
			n_i,
			k_i,
			n_t,
			k_t,
			Tphase,
			Ts,
			Tp,
			TsTp,
			Tfactor
			);
		if(isnan(Transmittance))   //Total Internal Refelection
			Transmittance = 0;
	}
	else Transmittance=0;

	double NomalizedReflectProbility=Reflectance/(Reflectance+Transmittance);
	double NomalizedTransmitProbility=1-NomalizedReflectProbility;

	double randomNum= rand()/(double)RAND_MAX;
	Mueller attenuation;
	if ( randomNum < NomalizedTransmitProbility )
	{
		getRefractionAttenuation(	currentIncomingDirection,
			SampledLocalNormal,
			n_i,
			k_i,
			n_t,
			k_t,
			currentOutgoingDirection,
			attenuation);
		currentStokes = attenuation * currentStokes ;

		if ( isNormalRefraction(currentIncomingDirection,currentOutgoingDirection))
		{
			currentIncomingDirection = currentOutgoingDirection;
			TraceRayDownwards();
		} 
		else
		{
			multiRelflectionHappen = true;
			return;
		}
	} 
	else
	{
		getReflectionAttenuation(	currentIncomingDirection,
			SampledLocalNormal,
			n_i,
			k_i,
			n_t,
			k_t,
			currentOutgoingDirection,
			attenuation);
		currentStokes = attenuation * currentStokes ;

		if ( isNormalReflection(currentIncomingDirection,currentOutgoingDirection))
		{
			currentIncomingDirection = currentOutgoingDirection;
			TraceRayUpwards();
		} 
		else
		{
			multiRelflectionHappen = true;
			return;
		}
	}
}

void RayTracing::TraceRayUpwards()
{
	if(PrintTraceInfo == true)
		printf(" %d ", currentLayer);

	if( currentLayer == 0  && currentOutgoingDirection.z > 0 )
		return;

	if(bounceNum == BounceLimit)
	{	
		return;
	}

	if (bounceNum < BounceLimit)
		bounceNum++;

	int layerNumOfTracedRay = currentLayer;
	double n_i,n_t,k_i,k_t;
	if (layerNumOfTracedRay > 0)
	{
		n_i = N[layerNumOfTracedRay];
		k_i = K[layerNumOfTracedRay];
		n_t = N[layerNumOfTracedRay-1];
		k_t = K[layerNumOfTracedRay-1];
	} 
	else
	{
		cout<<"layerNumOfTracedRay-1 < 0"<<endl;
	}


	Vector SampledLocalNormal;
	SampledLocalNormal = SampleAnisotropicNormal(layerNumOfTracedRay); // ??

	double localTheta=acos(abs(Dot(
		currentIncomingDirection,
		SampledLocalNormal)));


	double Reflectance=0,Transmittance=0;
	double Ts,Tp,TsTp,Tphase,Tfactor;
	double rs,rp,phaseS,phaseP,Fs,Fp;

	Reflectance = Fresnel_Reflection( localTheta, 
		n_i,
		k_i,
		n_t,
		k_t,
		rs,
		rp,
		phaseS,
		phaseP,
		Fs,
		Fp);

	if (k_t==0)  // next layer is not metallic or opaque AND there is not transimission at the last layer
	{
		Transmittance = Fresnel_Refraction(
			localTheta,
			n_i,
			k_i,
			n_t,
			k_t,
			Tphase,
			Ts,
			Tp,
			TsTp,
			Tfactor
			);
		if(isnan(Transmittance))   //Total Internal Refelection
			Transmittance = 0;
	}
	else Transmittance=0;

	double NomalizedReflectProbility=Reflectance/(Reflectance+Transmittance);
	double NomalizedTransmitProbility=1-NomalizedReflectProbility;

	double randomNum= rand()/(double)RAND_MAX;
	Mueller attenuation;
	if ( randomNum < NomalizedTransmitProbility )
	{
		getRefractionAttenuation(	currentIncomingDirection,
			SampledLocalNormal,
			n_i,
			k_i,
			n_t,
			k_t,
			currentOutgoingDirection,
			attenuation);
		currentStokes = attenuation * currentStokes ;

		if ( isNormalRefraction(currentIncomingDirection,currentOutgoingDirection))
		{
			currentIncomingDirection = currentOutgoingDirection;
			TraceRayUpwards();
		} 
		else
		{
			multiRelflectionHappen = true;
			return;
		}
	} 
	else
	{
		getReflectionAttenuation(	currentIncomingDirection,
			SampledLocalNormal,
			n_i,
			k_i,
			n_t,
			k_t,
			currentOutgoingDirection,
			attenuation);
		currentStokes = attenuation * currentStokes ;

		if ( isNormalReflection(currentIncomingDirection,currentOutgoingDirection))
		{
			currentIncomingDirection = currentOutgoingDirection;
			TraceRayDownwards();
		} 
		else
		{
			multiRelflectionHappen = true;
			return;
		}
	}
}

void RayTracing::getRefractionAttenuation(
	const Vector    incomingDirecion,
	const Vector    localNormal,
	const double    n_i,
	const double    k_i,
	const double    n_t,
	const double    k_t,
	Vector  & outgoingDirecion,
	Mueller & mu
	)
{
	if ( incomingDirecion.z < 0)
	{
		currentLayer++;
		if (currentLayer > NumOfLayer)
		{
			cout<<"Error.....currentLayer > TotalLayerNum"<<endl;
		}
	}
	else
	{
		currentLayer--;
	}
	bool TIR;
	GetSpecularRefractionDirection(incomingDirecion,
		localNormal,
		n_i,
		n_t,
		TIR,
		outgoingDirecion);
	double cosi=abs( Dot( incomingDirecion ,  localNormal ));
	double theta_i=acos(cosi);
	double Ts,Tp,TsTp,Tphase,Tfactor;
	Fresnel_Refraction(
		theta_i,
		n_i,
		k_i,
		n_t,
		k_t,
		Tphase,
		Ts,
		Tp,
		TsTp,
		Tfactor
		);

	double A,B,C,S;
	A=Tfactor*(Ts+Tp)/2;
	B=Tfactor*(Ts-Tp)/2;
	C=Tfactor*TsTp;
	S=Tfactor*0;

	double localToglobal_angle;
	double globalTolocal_angle;

	rotationAngle(
		incomingDirecion,
		outgoingDirecion,
		surfaceNormal,
		localNormal,
		globalTolocal_angle,
		localToglobal_angle
		);

	JonesMat R1=Jones_Rotator(globalTolocal_angle),R2=Jones_Rotator(localToglobal_angle);

	Mueller mR1(R1),mR2(R2);
	setFresnelMueller(mu,A,B,C,S);
	mu=mR2*mu;
	mu=mu*mR1;

}

void RayTracing::getReflectionAttenuation(
	const Vector    incomingDirecion,
	const Vector    localNormal,
	const double    n_i,
	const double    k_i,
	const double    n_t,
	const double    k_t,
	Vector  & outgoingDirecion,
	Mueller & mu
	)
{
	GetSpecularReflectionDirection(
		incomingDirecion,
		localNormal,
		n_i,
		n_t,
		outgoingDirecion);

	double cosi=abs( Dot( incomingDirecion ,  localNormal ));
	double theta_i=acos(cosi);

	double rs,rp,phaseS,phaseP,Fs,Fp;
	Fresnel_Reflection(
		theta_i,
		n_i,
		k_i,
		n_t,
		k_t,
		rs,
		rp,
		phaseS,
		phaseP,
		Fs,
		Fp
		);

	double A,B,C,S;
	A=(Fs+Fp)/2;
	B=(Fs-Fp)/2;
	C=cos(phaseS-phaseP)*sqrt(Fs*Fp);
	S=sin(phaseS-phaseP)*sqrt(Fs*Fp);

	double localToglobal_angle;
	double globalTolocal_angle;

	rotationAngle(
		incomingDirecion,
		outgoingDirecion,
		surfaceNormal,
		localNormal,
		globalTolocal_angle,
		localToglobal_angle
		);

	JonesMat R1=Jones_Rotator(globalTolocal_angle),R2=Jones_Rotator(localToglobal_angle);

	Mueller mR1(R1),mR2(R2);
	setFresnelMueller(mu,A,B,C,S);
	mu=mR2*mu;
	mu=mu*mR1;

}

void RayTracing::setFresnelMueller(Mueller &mu,double A,double B,double C,double S)
{
	mu.SetValue(0,0,A);
	mu.SetValue(1,1,A);
	mu.SetValue(0,1,B);
	mu.SetValue(1,0,B);
	mu.SetValue(2,2,C);
	mu.SetValue(3,3,C);
	mu.SetValue(2,3,S);
	mu.SetValue(3,2,-S);
}

Vector RayTracing::SampleAnisotropicNormal(int layer)
{
	int index = layer-1;
	float ex = Ex[index];
	float ey = Ey[index];

	Vector normal;
	SampleAnisotropic(ex,ey,normal);
	return normal;
}
