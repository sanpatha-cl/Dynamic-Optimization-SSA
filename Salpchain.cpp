#include "Salpchain.h"

#include <iostream>
#include "Global.h"

#include "Composition_DBG.h"
#include "Rotation_DBG.h"
#include <random>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>

using namespace std;

CSalp *CSalpchain::globalBest;
CSalp *CSalpchain::prevGlobalBest;

int CSalpchain::swarmCount=0;
int CSalpchain::popSize = 0;
CSalpchain *CSalpchain::subSwarm=0;
const int CSalpchain::subPopSize = 10;

CSalpchain::CSalpchain(void){
	
}
CSalpchain::CSalpchain(const CSalpchain &s){
	int i;
	popSize=s.popSize ;
	population = new CSalp[popSize];
	populationBest = new CSalp[popSize];

	for( i=0;i<popSize;i++){
		population[i]=s.population[i];
		populationBest[i]=s.populationBest[i];
	}
	
	bestPosition =s.bestPosition;
	
	avgFitness=s.avgFitness;
	stdError=s.stdError;

}
CSalpchain::CSalpchain(const int pop_size){
	popSize=pop_size;
	population = new CSalp[pop_size];
	populationBest = new CSalp[pop_size];
}
CSalpchain &CSalpchain::operator= (const CSalpchain &s){
	if(this==&s)return *this;

	int i;
	popSize=s.popSize ;
	for( i=0;i<popSize;i++){
		population[i]=s.population[i];
		populationBest[i]=s.populationBest[i];
	}
	bestPosition =s.bestPosition;
	avgFitness=s.avgFitness;
	stdError=s.stdError;
	return *this;
}
CSalpchain::~CSalpchain(void){
	if(population){
		delete [] population;
		population =0;
	}
	if(populationBest) {
		delete [] populationBest;
		populationBest =0;
	}
}

void CSalpchain::Initial(const int pop_size, ofstream &traceLog){
	if (Global::isLogEnable)
		traceLog << "Start ->CSalpchain::Initial" << endl;
	int i;

	popSize=pop_size;

	population = new CSalp[pop_size];
	for (i = 0; i < pop_size; i++) {
		population[i].Initialization(traceLog);
		population[i].m_nPopulationNo = i + 1;
	}
		
	populationBest = new CSalp[pop_size];
	for( i=0;i<pop_size;i++)
		populationBest[i]= population[i];

	sortPopulation();
	bestPosition = population[0];
	CSalpchain::globalBest = &bestPosition;
	CSalpchain::prevGlobalBest = CSalpchain::globalBest;

	if (Global::isLogEnable)
		traceLog << "End ->CSalpchain::Initial" << endl;
}

//Sort the Salp Population in Ascending order or Descending order depends on the MIN/MAX problem
void CSalpchain::sortPopulation()
{
	int minPos = 0;
	int nSalpNo = 0;
	CSalp tempParticle;
	switch (Global::optimization_type) {
	case MIN:
		for (int i = 1; i < popSize; i++)
		{
			for (int j = 0; j < popSize - i; j++)
			{
				if (population[j].fitness > population[j + 1].fitness)
				{
					tempParticle = population[j];
					nSalpNo = population[j].m_nPopulationNo;
					population[j] = population[j + 1];
					population[j].m_nPopulationNo = population[j + 1].m_nPopulationNo;
					population[j + 1] = tempParticle;
					population[j + 1].m_nPopulationNo = nSalpNo;
				}
			}
		}
		break;
	case MAX:
		for (int i = 1; i < popSize; i++)
		{
			for (int j = 0; j < popSize - i; j++)
			{
				if (population[j].fitness < population[j + 1].fitness)
				{
					tempParticle = population[j];
					nSalpNo = population[j].m_nPopulationNo;
					population[j] = population[j + 1];
					population[j].m_nPopulationNo = population[j + 1].m_nPopulationNo;
					population[j + 1] = tempParticle;
					population[j + 1].m_nPopulationNo = nSalpNo;
				}
			}
		}
		break;
	}
}

const int CSalpchain::findBest(void)const
{
	int index_best=0;
	for(int i=1;i<popSize;i++){
		if(population[i].Comparison(population[index_best])==1){
			index_best=i;
		}
	}
	return index_best;
}

void CSalpchain::printPopulationSalp(ofstream &traceLog)
{
	for (int i = 0; i < popSize; i++) {
		for (int j = 0; j < Global::g_dbg->Get_Dimension(); j++) {
			traceLog << population[i].pself.x[j] << "	";
		}
		traceLog << "( " << population[i].fitness << " )";
		traceLog << endl;
	}
	if (Global::isLogEnable)
		traceLog << "Food Position ( " << bestPosition.fitness << " )";
}

void CSalpchain::transposePopulationSalp(ofstream &traceLog)
{
	double tempVal = 0.0;
	int nDim = Global::g_dbg->Get_Dimension();
	int iIndex = 0, jIndex = 0;

	for (int i = 0; i < popSize; i++) {
		for (int j = i+1; j < Global::g_dbg->Get_Dimension(); j++) {
			tempVal = population[i].pself.x[j];
			population[i].pself.x[j] = population[j].pself.x[i];
			population[j].pself.x[i] = tempVal;
			//iIndex = i * nDim;
			//jIndex = j * nDim;
			//tempVal = population[iIndex].pself.x[j];
			//population[iIndex].pself.x[j] = population[jIndex].pself.x[i];
			//population[jIndex].pself.x[i] = tempVal;
		}
	}
}

void CSalpchain::Statistic(void)
{
	calcAverageFitness();
	stdError=0;
	for(int i=0;i<popSize;i++)
		stdError+=(population[i].fitness-avgFitness)*(population[i].fitness-avgFitness);
	stdError=sqrt(stdError/popSize);
}

void CSalpchain::calcAverageFitness(void){
	avgFitness=0;
	// calculate the average fitness of its Salp Swarm
	for(int i=0;i<popSize;i++)
		avgFitness+= population[i].fitness;
	avgFitness/=popSize;
}

void CSalpchain::personalBestUpdate(ofstream &traceLog){

	int b=findBest();
	if (population[b].Comparison(bestPosition) == 1)
		bestPosition = population[b];

	// update the best one of all swarms
	prevGlobalBest = globalBest;
	if (bestPosition.Comparison(*globalBest) == 1)
		globalBest = &bestPosition;
}
void CSalpchain::createSalpchain(CSalpchain & w, ofstream &traceLog){
	if (Global::isLogEnable)
		traceLog << "Start ->CSalpchain::createSalpchain" << endl;
	CSalpchain *sw=new CSalpchain[swarmCount+1];
	int i=0;
	for( i=0;i<swarmCount;i++){
		sw[i].Initial(subSwarm[i].popSize, traceLog);
		sw[i]=subSwarm[i];
	}
	sw[i].Initial(w.popSize, traceLog);
	sw[i]=w;
	delete [] subSwarm;
	swarmCount++;
	subSwarm= sw;

	if (Global::isLogEnable)
		traceLog << "End ->CSalpchain::createSalpchain" << endl;
}
void CSalpchain::deletePopulation(int index, ofstream &traceLog){
	if (Global::isLogEnable)
		traceLog << "Start ->CSalpchain::deletePopulation" << endl;
	if(swarmCount==1) {
		delete  [] subSwarm;
		subSwarm=0;
		swarmCount--;
		return ;
	}
	CSalpchain *sw=new CSalpchain[swarmCount-1];
	int i,j;
	for( i=0,j=0;j<swarmCount;i++,j++){
		if(j==index){
			i--;
			continue;	
		}
		sw[i].Initial(subSwarm[j].popSize, traceLog);
		sw[i]=subSwarm[j];
	}
	delete [] subSwarm;
	swarmCount--;
	subSwarm=sw;
	if (Global::isLogEnable)
		traceLog << "End ->CSalpchain::deletePopulation" << endl;
}
void CSalpchain::addSalp(CSalp &add){
	CSalp *sw=new CSalp[popSize+1];
	CSalp *swbest=new CSalp[popSize+1];
	int i;
	for( i=0;i<popSize;i++){
		sw[i]= population[i];
		swbest[i]= populationBest[i];
	}
	sw[i]=add;
	swbest[i]=add;
	
	delete [] population;
	population =0;
	delete [] populationBest;
	populationBest =0;
	population =sw;
	populationBest =swbest;
	popSize++;
}
void CSalpchain::deleteSalp(int del){
	if(popSize==1) {
		delete  [] population;
		delete []populationBest;
		populationBest =0;
		population =0;
		popSize--;
		return;
	}
	CSalp *sw=new CSalp[popSize-1];
	CSalp *swbest=new CSalp[popSize-1];
	int i,j;
	for( i=0,j=0;j<popSize;i++,j++){
		if(j==del){
			i--;
			continue;	
		}
		sw[i]= population[j];
		swbest[i]= populationBest[j];
	}
	delete [] population;
	population =0;
	delete []populationBest;
	populationBest =0;
	population =sw;
	populationBest =swbest;
	popSize--;
}

void CSalpchain::ReinitializeChain(CSalp* pSalpParticleArchive, int nArchiveNum, ofstream &traceLog) {
	
	int nIndexFoodPosition = findBest();
	CSalp tBest = population[nIndexFoodPosition];
	for (int i = 0; i < popSize; i++) {
		if (i >= 2 * subPopSize) { //Reinitializing the subPopulations more than 2*subPopSize
			population[i].Initialization(traceLog);
		}
		else { //Existing knowledge from archive for first (two) subPopulation using
			bool bIsSubPopNormal = false;
			if (i % subPopSize < subPopSize / 2) {
				bIsSubPopNormal = true;
			}
			if (i < subPopSize && nArchiveNum > 0 && i < nArchiveNum) {
				population[i] = pSalpParticleArchive[(int)(nArchiveNum*Global::uniform.Next())]; // archive
				population[i].calcSalpFitness(traceLog);
			}
			else {
				population[i] = tBest;
				population[i].calcSalpFitness(traceLog);
			}
		}
	}
	for (int i = 0; i < popSize; i++)
		populationBest[i] = population[i];
	bestPosition = population[findBest()];
}

void CSalpchain::subBest(int sizeSub)
{
	CSalp tmp;
	int n = sizeSub; // 
	for (int k = 0; k < popSize; k = k + n) {
		int ibest = k;
		for (int i = k + 1; i < k + n; i++) {
			if (population[i].Comparison(population[ibest]) == 1) {
				ibest = i;
			}
		}
		if (ibest != k) {
			tmp = population[k];
			population[k] = population[ibest];
			population[ibest] = tmp;
		}
	}
}

double CSalpchain::bestMean(int nIndex)
{
	double dSum = 0.0, dBestMean = 0.0;
	for (int i = nIndex; i < popSize; i++) {
		dSum = dSum + population[i].fitness;
	}
	dBestMean = dSum / popSize;
	return dBestMean;
}

void CSalpchain::SalpSwarm_Locomate(double globalOptima, int &fit_eva, double &r_value, CSalp* pSalpParticleArchive, int nArchiveNum, ofstream &debugInfo, ofstream &traceLog) {
	int r1, r2, r3;
	int nBestIndex = findBest(), nLocalBestIndex = 0;
	bool bSalpAge = false;
	int NP = popSize, nCurSalpAge = 0, nSubPopReinitialize = 0, nIndexjDE = 0;
	int nSubPop = NP / subPopSize; // number of subPopulations
	double dBorderValue = Global::uniform.Next();
	double F = 0.5, CR = 0.9;
	CSalp tempSalp;
	double dBestMean = 0.0;

	//Find the local best in Sub Population
	subBest(popSize / nSubPop);

	//SSA Implementation Main loop for the Population Size
	for (int i = 0; i < popSize; i++) {

		int tmpNP;
		tmpNP = NP / nSubPop;
		int pomoc = (i / tmpNP)*tmpNP;

		do                        /* Pick a random population member */
		{                         /* Endless loop for NP < 2 !!!     */ 
			r1 = (int)(Global::uniform.Next() * tmpNP + pomoc);
		} while (r1 == i);

		do                        /* Pick a random population member */
		{                         /* Endless loop for NP < 3 !!!     */
			r2 = (int)(Global::uniform.Next() * tmpNP + pomoc);
		} while ((r2 == i) || (r2 == r1));

		do                        /* Pick a random population member */
		{                         /* Endless loop for NP < 4 !!!     */
			r3 = (int)(Global::uniform.Next() * tmpNP + pomoc);
		} while ((r3 == i) || (r3 == r1) || (r3 == r2));

		tempSalp = population[i];
		CSalp::accMinionParam = CSalp::accMinionParam * Global::uniform.Next();

		// F jDE
		if (Global::uniform.Next() < 0.1) {
			CSalp::m = tempSalp.ctrlParamMutation = 0.46 + Global::uniform.Next() * (1.0 - 0.46);
		}
		else {
			CSalp::m = tempSalp.ctrlParamMutation;
		}

		// CR jDE
		if (Global::uniform.Next() < 0.1) {
			CR = tempSalp.ctrlParam = Global::uniform.Next()*1.0;
		}
		else {
			CR = tempSalp.ctrlParam;
		}

		nLocalBestIndex = (i / (NP / nSubPop)) * (NP / nSubPop);

		//Get the current salp age and reset the bSalpAge variable
		nCurSalpAge = tempSalp.salpAge;
		bSalpAge = false;

		//If the current salp age is greate than 25 and neither local best not global best then re-initialize i-th individual (Ref.370 Alg.1 line 4)
		if (true && (i != nLocalBestIndex) && (i != nBestIndex) && nCurSalpAge > 25 && Global::uniform.Next() < 0.1 && (i != nBestIndex)) {
			bSalpAge = true;
		}

		//If the current salp age is greate than 30 and is local best but not global best then re-initialize subPopulation that contain i (Ref. 370 Alg.1 line 2)  
		if (i == nLocalBestIndex && (i != nBestIndex) && nCurSalpAge > 30 && Global::uniform.Next() < 0.1) {
			nSubPopReinitialize = subPopSize;
		}//If the current salp is local best but not global best and two local best close to each other then re-initialize subPopulation that contain i (Ref. 370 Alg.2)
		else if (i == nLocalBestIndex && (i != nBestIndex)) {
			for (int k = 0; k < NP; k += subPopSize) {
				if (i == k /*|| indexOfBest == k */)
					continue;
				if (tempSalp.pself.Distance(population[k].pself.x) < 0.05) // 0.01 
					nSubPopReinitialize = subPopSize;
			}
		}

		//If the current salp is neither local best, nor global best and does not belong to best subpopulation and distance between i and local best is small then re-initialize i-th individual (Ref. 370 Alg.3)
		if (i != nLocalBestIndex && (i != nBestIndex) && (i / subPopSize != nBestIndex / subPopSize) && tempSalp.pself.Distance(population[nLocalBestIndex].pself.x) < 0.001 && tempSalp.salpAge > 15) {
			bSalpAge = true;
		}
		CSalp salp1 = population[r1];
		CSalp salp2 = population[r2];
		CSalp salp3 = population[r3];
		dBestMean = bestMean(i);
		tempSalp.moveCAQSSA(salp1, salp2, salp3, bestPosition, bSalpAge, CR, dBestMean, dBorderValue, traceLog);
		
		if (bSalpAge) {//Use existing knowledge for the small movements
			if (Global::uniform.Next() < 0.5 && i < subPopSize && nArchiveNum > 0 && i < nArchiveNum) {
				tempSalp = pSalpParticleArchive[(int)(nArchiveNum*Global::uniform.Next())];
				tempSalp.calcSalpFitness(traceLog);
			}
		}

		if (nSubPopReinitialize > 0) {
			nSubPopReinitialize--;

			if (true && Global::uniform.Next() < 0.5 && i < subPopSize && nArchiveNum > 0 && i < nArchiveNum) {
				//Use existing knowledge for the small movements
				tempSalp = pSalpParticleArchive[(int)(nArchiveNum*Global::uniform.Next())];
				tempSalp.calcSalpFitness(traceLog);
			}
			else
			{
				tempSalp.Initialization(traceLog);
			}
			bSalpAge = true;
		}
		else {
			//Re-calculate the fitness of each salp
			tempSalp.calcSalpFitness(traceLog);
		}

		if (tempSalp.Comparison(population[i]) == 1 || bSalpAge) {
			if (bSalpAge) {
				population[i] = tempSalp;
				population[i].salpAge = 0;
			}
			else {
				// Ref. 370 Alg.4 Improvement and aging
				if (tempSalp.pself.Distance(population[i].pself.x) < 0.01 || fabs(tempSalp.fitness - population[i].fitness) < 0.1) {// if small improvement 
					population[i] = tempSalp;
					if (population[i].salpAge > 20) {
						population[i].salpAge = 20;
					}
				}
				else {
					population[i] = tempSalp;
					if (population[i].salpAge > 5) {
						population[i].salpAge = 5;
					}
				}
			}
		}

		//Increment the age
		population[i].salpAge++;

		//Compare it with previous best and update
		if (population[i].Comparison(populationBest[i]) == 1) {
			populationBest[i] = population[i];
			nBestIndex = findBest();
			if (i == nBestIndex) personalBestUpdate(traceLog);
		}

		//Recording the r-value for report
		fit_eva++;
		if (fit_eva%Global::sample_frequency == 0) {

			if (Global::optimization_type == MIN)
				r_value += (1 - globalOptima / CSalpchain::globalBest->fitness);
			else
				r_value += (1 - CSalpchain::globalBest->fitness / globalOptima);
		}
	}
}

void CSalpchain::BuildSalp_Chain(double globalOptima, int& fit_eva, double& r_value, CSalp* pSalpParticleArchive, int nArchiveNum, ofstream& debugInfo, ofstream& traceLog) {
	//Random numbers 
	double lowerCoordinate, upperCoordinate;

	for (int i = 0; i < popSize; i++) {

		if (i < popSize / 2) //First Salp Leader
		{
			for (int j = 0; j < Global::g_dbg->Get_Dimension(); j++) {

				lowerCoordinate = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].lower;
				upperCoordinate = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].upper;

				CSalp::c2 = Global::uniform.Next();
				CSalp::c3 = Global::uniform.Next();

				if (CSalp::c3 < 0.5) // Eq. 3.1 from SSA
				{
					population[i].pself.x[j] = bestPosition.pself.x[j] + (CSalp::c1 * ((upperCoordinate - lowerCoordinate) * CSalp::c2 + lowerCoordinate));
				}
				else
				{
					population[i].pself.x[j] = bestPosition.pself.x[j] - (CSalp::c1 * ((upperCoordinate - lowerCoordinate) * CSalp::c2 + lowerCoordinate));
				}

				if (population[i].pself.x[j] > upperCoordinate || population[i].pself.x[j] < lowerCoordinate)
				{
					population[i].pself.x[j] = lowerCoordinate + (upperCoordinate - lowerCoordinate);
				}
			}
		}
		else if (i >= popSize / 2 && i < popSize + 1) //Equation 3.4 from SSA paper for Follower SALP
		{
			for (int j = 0; j < Global::g_dbg->Get_Dimension(); j++) {

				double lowerCoordinate = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].lower;
				double upperCoordinate = static_cast<Real_DBG*>(Global::g_dbg)->Get_Boundary()[j].upper;

				population[i].pself.x[j] = (population[i - 1].pself.x[j] + population[i].pself.x[j]) / 2;

				if (population[i].pself.x[j] > upperCoordinate || population[i].pself.x[j] < lowerCoordinate)
				{
					population[i].pself.x[j] = lowerCoordinate + (upperCoordinate - lowerCoordinate);
				}
			}
		}

		//Re-calculate the fitness of each salp
		population[i].calcSalpFitness(traceLog);

		fit_eva++;

		if (fit_eva % Global::sample_frequency == 0) {

			if (Global::optimization_type == MIN)
				r_value += (1 - globalOptima / CSalpchain::globalBest->fitness);
			else
				r_value += (1 - CSalpchain::globalBest->fitness / globalOptima);
		}

		if (population[i].Comparison(populationBest[i]) == 1) {
			populationBest[i] = population[i];
		}
	}
}