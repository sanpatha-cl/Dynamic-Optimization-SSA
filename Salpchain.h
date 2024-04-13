#pragma once

#include "Salp.h"
#include "Global.h"
#include <fstream>
using namespace std;

class CSalpchain
{
public:
	static int popSize;  // swarm size  
	const static int subPopSize;
	
	double avgFitness; // agerage fitness of whole swarm
	double stdError;   //standard deviation of whole swarm
	
	CSalp bestPosition;  // the best one of swarm
	
	CSalp *population;	 // the population
	CSalp *populationBest; //population of previous best position
	int m_nPopulationNo;
	
public:

	static CSalp *globalBest; // global best one found by LPSO
	static CSalp *prevGlobalBest; // global best one found by LPSO

	// the total number of swarms
	static int swarmCount; 
	static CSalpchain *subSwarm;

public:
	CSalpchain(void);
	CSalpchain(const CSalpchain &s);
	CSalpchain(const int pop_size);
	~CSalpchain(void);
	void Initial(const int pop_size, ofstream &);
	void ReinitializeChain(CSalp*, int, ofstream &);
	const int findBest(void)const;

	void SalpSwarm_Locomate(double best, int&, double &, CSalp*, int, ofstream &, ofstream &);
	void BuildSalp_Chain(double best, int&, double&, CSalp*, int, ofstream&, ofstream&);
	double bestMean(int nIndex);
	void subBest(int);  // local best v sub-populacijah

public:
	void Statistic(void);
public:
	void calcAverageFitness(void);
	void personalBestUpdate(ofstream &);
	static void createSalpchain(CSalpchain &w, ofstream &);
	void sortPopulation();
	static void deletePopulation(int index, ofstream &);
	void addSalp(CSalp &add);
	void deleteSalp(int del);
	CSalpchain &operator= (const CSalpchain &s);
	void printPopulationSalp(ofstream &);
	void transposePopulationSalp(ofstream &);
};
