#include "BruteForceSimulator.h"


BruteForceSimulator::BruteForceSimulator(void)
{
	N.push_back(1);  // Top layer: AIR
	K.push_back(0);

	NumOfLayer         = 0;
	NumOfIteration     = 0;
	effectiveSampleNum = 0;

	surfaceNormal			 = Vector(0,0,1);

	for (int i=0;i<91;i++)
	{
		for(int j=0;j<360;j++)
			reflection_data[i][j] = new StokesVector();
	}

	//ReadSettingFile();
	ReadAnisotropicSettingFile();
	srand(1);

	
}


BruteForceSimulator::~BruteForceSimulator(void)
{
}

void BruteForceSimulator::Calculate(void)
{
	#pragma omp parallel for
	for (int i=0; i < NumOfIteration; i++)
	{
		RayTracing rt = RayTracing();
		Vector outgoingDirection;
		StokesVector outgoingStokes;
		rt.InitializeAndCompute(IncidenceDirection,incidenceStokes,N,K,Ex,Ey,outgoingDirection,outgoingStokes);
		CollectData(outgoingDirection,outgoingStokes);
	}

	SaveToFile();
}


void BruteForceSimulator::CollectData(Vector direction,StokesVector sv)
{
	Vector d = direction;
	int thetaDeg = 0 , phiDeg = 0;

	float x = d.x , y = d.y , z = d.z;

	if (z>0)
	{

		float costheta = CosBetweenVector ( d , surfaceNormal );
		float temp = acosf ( costheta );
		temp = temp * MATH_RAD_TO_DEG;
		//thetaDeg=(int)temp;
		//int temp1=(int)temp;
		thetaDeg = ( int ) ( temp + 0.5 );  // >=temp?temp1+1:temp1;
		thetaDeg = Clamp ( thetaDeg , 0 , 90 );


		float absphi = abs ( atanf ( y / x ) )*MATH_RAD_TO_DEG;
		//int temp2=(int)absphi;
		phiDeg = ( int ) ( absphi + 0.5 );//(2*temp2+1)>=2*absphi?temp2+1:temp2;

		if (x > 0 && y >= 0)
			phiDeg = phiDeg;
		if (x > 0 && y <= 0)
			phiDeg = 360 - phiDeg;
		if (x < 0 && y <= 0)
			phiDeg = 180 + phiDeg;
		if (x < 0 && y >= 0)
			phiDeg = 180 - phiDeg;

		phiDeg=phiDeg%360;
		phiDeg=Clamp(phiDeg,0,359);

		*reflection_data[thetaDeg][phiDeg]=*reflection_data[thetaDeg][phiDeg]+sv;
	}

}


void BruteForceSimulator::SaveToFile()
{
	char I[]="I.txt";
	char Q[]="Q.txt";
	char U[]="U.txt";
	char V[]="V.txt";
	char DOP[]="DOP.txt";
	char AOP[]="AOP.txt";

	ofstream Ifs(I);
	ofstream Qfs(Q),Ufs(U),Vfs(V),DOPfs(DOP),AOPfs(AOP);

	for (int i=0;i<91;i++)
	{
		for(int j=0;j<360;j++)
		{
			StokesVector s=*reflection_data[i][j];
			Ifs<<s.I()<<" ";
			Qfs<<s.Q()<<" ";
			Ufs<<s.U()<<" ";
			Vfs<<s.V()<<" ";
			DOPfs<<s.DOP()<<" ";
			AOPfs<<s.eta()<<" ";
		}
		Ifs<<endl;
		Qfs<<endl;
		Ufs<<endl;
		Vfs<<endl;
		DOPfs<<endl;
		AOPfs<<endl;
	}

}


void BruteForceSimulator::SetIncidence(float theta, float phi, StokesVector stokes)
{
	theta	= theta*MATH_DEG_TO_RAD;
	phi		= phi*MATH_DEG_TO_RAD;
	float sintheta=sin(theta),costheta=cos(theta);

	IncidenceDirection=Vector(sintheta * cosf(phi),sintheta * sinf(phi),-costheta);
	incidenceStokes=stokes;
	
}

void BruteForceSimulator::AddAnisotropicLayer(float n, float k,float ex,float ey)
{
	N.push_back(n);
	K.push_back(k);
	Ex.push_back(ex);
	Ey.push_back(ey);

	NumOfLayer++;
}

void BruteForceSimulator::ReadAnisotropicSettingFile()
{
	float incidentAngle,I,Q,U,V;

	ifstream x("setting.txt");
	if (!x.is_open())
	{
		cout<<"Couldn't find SETTING files"<<endl;
	}

	x>>incidentAngle;

	x>>I>>Q>>U>>V;

	this->SetIncidence(incidentAngle,0,StokesVector(I,Q,U,V));

	x>>NumOfIteration;

	cout<<"Incident Angle : "<<incidentAngle<<endl;
	cout<<"Incident Stokes : "<<I<<" "<<Q<<" "<<U<<" "<<V<<endl;
	cout<<"Number Of Iteration : "<<NumOfIteration<<endl;

	int c;
	x>>c;

	cout<<"Number of Layer : "<<c<<endl;
	cout<<"Layer 0 (AIR)  IOR: ( 1, 0 )"<<endl;
	for (int j=0;j<c;j++)
	{
		float n,k,ex,ey;
		x>>n>>k>>ex>>ey;

		this->AddAnisotropicLayer(n,k,ex,ey);
		cout<<"Layer "<<j+1<< " IOR: ( "<<n<<", "<<k<<" )"<<endl;
	}

	x.close();
}





