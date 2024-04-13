#pragma once
#include "Point.h"
#include <fstream>
using namespace std;
class CSalp{
public:
	CPoint pself;
	CPoint pnearest; // the closest point of other swarm or local optimum if it is within their radius
	double fitness;
	int m_nPopulationNo;

	//For DE
	double ctrlParamMutation, ctrlParam;
	int salpAge;

public:
	static double *vmax;
	static double c1,c2, c3; //accelerators 
	static double wMax, accMinionParam; //accelerators 
	static double aAlpha, bBeta, r, u; //accelerators
	static double r1, r2, r3; //accelerators
	static double m;     // inertia weight
public:
	CSalp();
	CSalp(CSalp &p);
	~CSalp();
	void Initialization(ofstream &);
	void calcSalpFitness(ofstream &);
	CSalp & operator =(const CSalp &p);
	double Velocity();
	int Comparison(const CSalp &p);
	void moveCAQSSA(const CSalp &, const CSalp &, const CSalp &, const CSalp &, bool, double, double, double, ofstream &);
};
