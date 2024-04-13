#include "Salpchain.h"
#include "Salp.h"
#include "Global.h"
#include <iostream>

#include "Composition_DBG.h"
#include "Rotation_DBG.h"

using namespace std;

double *CSalp::vmax = 0;
double CSalp::c1 = 1.4960;
double CSalp::c2 = 1.4960;
double CSalp::c3;
double CSalp::accMinionParam = 0;
double CSalp::wMax = 1.5960;
double CSalp::r1;
double CSalp::r2;
double CSalp::r3;
double CSalp::m = 0.7298f;
double CSalp::aAlpha;
double CSalp::bBeta;
double CSalp::r;
double CSalp::u;

CSalp::CSalp(){
	
}
CSalp::~CSalp(){

}

void CSalp::Initialization(ofstream &traceLog){
	int i;
	Real rTemp;
	double u,l;
	for (i = 0; i < Global::g_dbg->Get_Dimension(); i++) {

		u = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[i].upper;
		l = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[i].lower;

		rTemp = Global::uniform.Next();
		pself.x[i] = l + (u - l)*rTemp;

		if (u - l > 2 * vmax[i])
		{
			rTemp = Global::uniform.Next();
			pself.v[i] = -vmax[i] + 2 * vmax[i] * rTemp;
		}
		else
		{
			rTemp = Global::uniform.Next();
			pself.v[i] = l + (u - l)*rTemp;
		}
	}
	calcSalpFitness(traceLog);
	ctrlParamMutation = 0.5;
	ctrlParam = 0.9;
	salpAge = 0;
}
void CSalp::calcSalpFitness(ofstream &traceLog){

	//Calculating the fitness
	fitness=static_cast<Real_DBG*>(Global::g_dbg)->Evaluation(pself.x);
	//traceLog << "CSalp::calcSalpFitness() : " << fitness << endl;

}
CSalp & CSalp::operator=(const CSalp &p){
	if(this==&p) return *this;
	pself=p.pself;
	fitness=p.fitness;
	
	pnearest=p.pnearest;

	ctrlParamMutation = p.ctrlParamMutation;
	ctrlParam = p.ctrlParam;
	salpAge = p.salpAge;
	
	return *this;
}
CSalp::CSalp(CSalp &p){
	pself=p.pself;
	fitness=p.fitness;
	pnearest=p.pnearest;
	ctrlParamMutation = p.ctrlParamMutation;
	ctrlParam = p.ctrlParam;
	salpAge = p.salpAge;
}

double CSalp::Velocity(){
	int i;
	double ve=0;
	
	for( i=0;i<Global::g_dbg->Get_Dimension();i++)
		ve+=pself.v[i]*pself.v[i];
	if(ve==0.0) return 0;
	return sqrt(ve);
}

int CSalp::Comparison(const CSalp &p){
	int flag;
	switch(Global::optimization_type){
	case MIN:
		if(fitness<p.fitness) flag=1;
		else if(fitness==p.fitness) flag= 0;
		else flag= -1;
		break;
	case MAX:
		if(fitness<p.fitness) flag= -1;
		else if(fitness==p.fitness) flag= 0;
		else flag= 1;
		break;
	}
	return flag;
}

void CSalp::moveCAQSSA(const CSalp & salp1, const CSalp &salp2, const CSalp &salp3, const CSalp &gbest, bool bAge, double CR, double dBestMean, double dBorderValue, ofstream &traceLog) {
	
	int nIndexjDE = (int)(Global::uniform.Next() * Global::g_dbg->Get_Dimension());

	for (int j = 0; j < Global::g_dbg->Get_Dimension(); j++) {

		double lowerCoordinate = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].lower;
		double upperCoordinate = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].upper;

		//Change at least one parameter in case of AGE
		if ((Global::uniform.Next() < CR) || j == (Global::g_dbg->Get_Dimension() - 1))
		{
			double lBound = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[nIndexjDE].lower;
			double uBound = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[nIndexjDE].upper;
			r1 = Global::uniform.Next();
			r2 = Global::uniform.Next();
			r3 = Global::uniform.Next();
			aAlpha = ((r1 * salp1.pself.x[nIndexjDE]) + (r2 * salp2.pself.x[nIndexjDE]) + (r3 * salp3.pself.x[nIndexjDE])) / (r1 + r2 + r3);
			
			//With Novel equation to update only followers
			r = Global::uniform.Next();
			u = 3 * 0.96 * (1 - 0.96) * Global::uniform.Next();
			c3 = Global::uniform.Next();
			if (c3 >= 0.5){
				pself.x[nIndexjDE] = aAlpha + bBeta * (gbest.pself.x[nIndexjDE] - pself.x[nIndexjDE]) * log(r / u);
			}
			else if (c3 < 0.5) {
				pself.x[nIndexjDE] = aAlpha - bBeta * (gbest.pself.x[nIndexjDE] - pself.x[nIndexjDE]) * log(r / u);
			}
			pself.x[nIndexjDE] = salp1.pself.x[nIndexjDE] + accMinionParam * (aAlpha - pself.x[nIndexjDE]) + m * (salp2.pself.x[nIndexjDE] - salp3.pself.x[nIndexjDE]);
			
			if (bAge) {  //Re-initializing the i-th individual
				pself.x[nIndexjDE] = uBound + (lBound - uBound) * Global::uniform.Next();
			}

			//Borders (check new values)
			if (dBorderValue < 0.5) {
				if (pself.x[nIndexjDE] > uBound) {
					pself.x[nIndexjDE] = 2 * uBound - pself.x[nIndexjDE];
				}
				if (pself.x[nIndexjDE] < lBound) {
					pself.x[nIndexjDE] = 2 * lBound - pself.x[nIndexjDE];
				}
			}
			else {
				if (pself.x[nIndexjDE] > uBound) {
					pself.x[nIndexjDE] = uBound;
				}
				if (pself.x[nIndexjDE] < lBound) {
					pself.x[nIndexjDE] = lBound + (uBound - lBound)*Global::uniform.Next();
				}
			}
		}
		nIndexjDE = (nIndexjDE + 1) % Global::g_dbg->Get_Dimension();

		if (pself.x[j] > upperCoordinate || pself.x[j] < lowerCoordinate)
		{
			pself.x[j] = lowerCoordinate + (upperCoordinate - lowerCoordinate)*Global::uniform.Next();;
		}
	}
}