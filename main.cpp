#include "stdafx.h"
#include "Global.h"
#include "Composition_DBG.h"
#include "Rotation_DBG.h"
#include "Salpchain.h"
#include <fstream>
#include <iostream>
#include <conio.h>
#include <random>
using namespace std;


void Initial_Change_Counter() {
	for (int i = 0; i < Num_Change_Type; i++) {
		Global::change_type[i].type = (Change_type)i;
		Global::change_type[i].counter = 0;
	}
}
void SalpSwarmAlgo_Initial(ofstream &traceLog) {
	if (Global::isLogEnable)
		traceLog << "Start ->Main : SalpSwarmAlgo_Initial" << endl;
	//particle velocity initial
	if (CSalp::vmax != 0) {
		delete[] CSalp::vmax;
		CSalp::vmax = 0;
	}

	double x = 0.0;
	CSalp::vmax = new double[Global::num_dim];
	for (int i = 0; i < Global::num_dim; i++) {
		x = (Global::boundary.upper - Global::boundary.lower) / 2.;
		CSalp::vmax[i] = x;
	}
	if (Global::isLogEnable)
		traceLog << "End ->Main : SalpSwarmAlgo_Initial" << endl;
}
void System_Initial(int t, double seed, const int num_dim, ofstream &traceLog) {
	
	if (Global::isLogEnable)
		traceLog << "Start ->Main : System_Initial" << endl;
	
	static LGM_mixed urng(seed);
	Random::Set(urng);

	Global::min_dimension = 5;
	Global::max_dimension = 15;
	Global::num_dim = num_dim;
	while (Global::num_dim < 5 || Global::num_dim>15) {
		cout << "please input dimensions within [5,15]:";
		cin >> Global::num_dim;
		getchar();
	}
	Global::boundary.Set_Boundary(-5, 5);
	if (t == Num_Change_Type) {
		Global::flag_dimension_change = true;
		Global::change = (Change_type)(u_random);
	}
	else {
		Global::flag_dimension_change = false;
		Global::change = (Change_type)(t);
	}
	if (t == recurrent || t == recurrent_noisy)
		Global::periodicity = 12;
	else
		Global::periodicity = 0;

	Initial_Change_Counter();

	SalpSwarmAlgo_Initial(traceLog);
	
	if (Global::isLogEnable)
		traceLog << "End ->Main : System_Initial" << endl;
}
void Continuous_Setting(Real_DBG *p) {
	p->Set_Boundary(Global::boundary);
	double *t = new double[Global::num_peakorfun];

	for (int i = 0; i < Global::num_peakorfun; i++) {
		if (p->Get_Change_Type() == chaotic)
			t[i] = Global::min_height + (Global::max_height - Global::min_height)*Global::uniform.Next();
		else
			t[i] = 50;
	}
	p->Set_Height(t);
	p->Set_Height_Severity(5);

	double **position;
	position = new double*[Global::num_peakorfun];
	for (int i = 0; i < Global::num_peakorfun; i++)
		position[i] = new double[Global::num_dim];
	for (int i = 0; i < Global::num_peakorfun; i++)
		for (int j = 0; j < Global::num_dim; j++) {
			position[i][j] = p->Get_Boundary()[j].lower + (p->Get_Boundary()[j].upper - p->Get_Boundary()[j].lower)*Global::uniform.Next();
		}
	p->Set_Position(const_cast<const double **>(position));
	delete[] t;
	for (int i = 0; i < Global::num_peakorfun; i++)
		delete[]position[i];
	delete[] position;
}
void Rotation_DBG_Setting(Rotation_DBG *p) {

	Continuous_Setting(p);
	p->Set_Width_Severity(0.5);
	p->Set_Weight(5);		// between (1,10)
	p->Calculate_Global_Optima();
}
void Composition_DBG_Setting(Composition_DBG * p, const int f) {

	Continuous_Setting(p);

	Fun_Name *basic_fun = new Fun_Name[Global::num_peakorfun];
	switch (f) {
	case 1:
		for (int i = 0; i < Global::num_peakorfun; i++) basic_fun[i] = Sphere;
		break;
	case 2:
		for (int i = 0; i < Global::num_peakorfun; i++) basic_fun[i] = Rastrigin;
		break;
	case 3:
		for (int i = 0; i < Global::num_peakorfun; i++) basic_fun[i] = Griewank;
		break;
	case 4:
		for (int i = 0; i < Global::num_peakorfun; i++) basic_fun[i] = Ackley;
		break;
	case 5:
		basic_fun[0] = Sphere;		basic_fun[1] = Sphere;
		basic_fun[2] = Rastrigin;		basic_fun[3] = Rastrigin;
		basic_fun[4] = Weierstrass;	basic_fun[5] = Weierstrass;
		basic_fun[6] = Griewank;		basic_fun[7] = Griewank;
		basic_fun[8] = Ackley;		basic_fun[9] = Ackley;
		break;
	}
	p->Set_Basic_Function(basic_fun);
	double *t = new double[Global::num_peakorfun];
	for (int i = 0; i < Global::num_peakorfun; i++)t[i] = 1.;
	p->Set_Coverge_Sevrity(t);
	p->Set_Stretch_Severity();
	p->Set_Rotation_Matrix();
	p->Calculate_Global_Optima();
	delete[]basic_fun;
	delete[] t;

}


void System_Setting(General_DBG * g_dbg, const int f) {

	g_dbg->Set_Change_Frequency(Global::change_frequency);
	g_dbg->Set_Change_Type(Global::change);

	g_dbg->Set_Periodicity(Global::periodicity);

	g_dbg->Set_Dimension_Change(Global::flag_dimension_change);


	if (g_dbg->Get_Change_Type() == chaotic) {
		Global::chaotic_constant = 3.67;
		while (Global::chaotic_constant > 4 || Global::chaotic_constant < 1) {
			cout << "invalid value of chaotic_constant,reset please" << endl;
			cin >> Global::chaotic_constant;
			getchar();
		}
	}

	if (g_dbg->Get_Change_Type() == recurrent_noisy)
		g_dbg->Set_RR_Severity(0.8f);

	if (Composition_DBG * p = dynamic_cast<Composition_DBG*>(g_dbg)) {
		Composition_DBG_Setting(p, f);
	}
	else if (Rotation_DBG * p = dynamic_cast<Rotation_DBG*>(g_dbg)) {
		Rotation_DBG_Setting(p);
	}
}
void Update_Global_Best(ofstream &traceLog) {
	
	//If the fitness value is better than the best fitness value(pBest) in history, set current value as the new pBest
	CSalpchain::prevGlobalBest = CSalpchain::globalBest;
	CSalpchain::globalBest = &CSalpchain::subSwarm[0].bestPosition;
	
	for (int k = 1; k < CSalpchain::swarmCount; k++)
		if (CSalpchain::globalBest->Comparison(CSalpchain::subSwarm[k].bestPosition) == -1) {
			CSalpchain::prevGlobalBest = CSalpchain::globalBest;
			CSalpchain::globalBest = &CSalpchain::subSwarm[k].bestPosition;
		}
}

void printPopulation(ofstream &traceLog)
{
	for (int k = 0; k < CSalpchain::swarmCount; k++) {
		for (int i = 0; i < Global::max_popsize; i++)
		{
			traceLog << "P" << CSalpchain::subSwarm[k].population[i].m_nPopulationNo <<" :  ";
			for (int j = 0; j < Global::g_dbg->Get_Dimension(); j++)
			{
				traceLog << CSalpchain::subSwarm[k].population[i].pself.x[j] << "	";
			}
			traceLog << "( " << CSalpchain::subSwarm[k].population[i].fitness << " )";
			traceLog << endl;
		}
	}
}

void Generate_SalpSwarm(ofstream &traceLog) {
	if (Global::isLogEnable)
		traceLog << "Start ->Main : Generate_SalpSwarm" << endl;
	
	for (int j = 0; j < Global::max_popnum; j++) {

		CSalpchain *newp = new CSalpchain(Global::max_popsize);
		//Initializing and calculating the fitness of each particle
		newp->Initial(newp->popSize, traceLog);
		CSalpchain::createSalpchain(*newp, traceLog);
	}

	if (Global::isLogEnable)
		traceLog << "End ->Main : Generate_SalpSwarm" << endl;
}
void Delete_SalpSwarm(ofstream &traceLog) {
	if (Global::isLogEnable)
		traceLog << "Start ->Main : Delete_SalpSwarm" << endl;
	for (int j = 0; j < CSalpchain::swarmCount; j++) {
		CSalpchain::deletePopulation(j, traceLog);
		j--;
	}
	if (Global::isLogEnable)
		traceLog << "End ->Main : Delete_SalpSwarm" << endl;
}

void printPopulationBest(ofstream &traceLog)
{
	for (int k = 0; k < CSalpchain::swarmCount; k++) {
		for (int i = 0; i < Global::max_popsize; i++)
		{
			for (int j = 0; j < Global::g_dbg->Get_Dimension(); j++)
			{
				traceLog << CSalpchain::subSwarm[k].populationBest[i].pself.x[j] << "	";
			}
			traceLog << "( " << CSalpchain::subSwarm[k].populationBest[i].fitness << " )";
			traceLog << endl;
		}
	}
}

void SalpSwarmAlgo(double **best_rel_fit, double **best_abs_fit, double **fit, double **relative, General_DBG * g_dbg, int *number_dimension, const int num_run, ofstream &debugInfo, ofstream &traceLog) {
	
	Generate_SalpSwarm(traceLog); //Initialize and calculate the fitness of each SALP

	CSalp salpParticleArchive[100];
	int nArchiveIndexNum = 0;

	if (!Global::isLogEnable) {
		traceLog << "Initial Popolation : " << endl;
		printPopulation(traceLog);
	}

	int max_gen;
	int fit_evas;
	double r_value;
	int iChangeCounter = 1;
	int numChange = 1;
	for (int i = 0; i < Global::num_change; i++) {

		max_gen = Global::change_frequency*Global::num_dim / Global::max_popsize;
		
		number_dimension[i] = Global::num_dim;
		
		if (num_run == 0) {// allocate memory to fit and relative for the first run
			
			fit[i] = new double[max_gen];
			
			relative[i] = new double[max_gen];
			
			for (int j = 0; j < max_gen; j++) {
				fit[i][j] = 0;
				relative[i][j] = 0;
			}
		}
		
		fit_evas = 0;
		r_value = 0;
		
		debugInfo << endl << "PASS : " << numChange++ << endl;

		Global::gen_Number = 0;
		
		//Relative calculation for Framework
		for (int l = 0; l < max_gen ; l++) {

			cout << i << " " << l << " : " << CSalpchain::globalBest->fitness << " " << static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima() << endl;

			fit[i][l] += CSalpchain::globalBest->fitness;
			if (Global::optimization_type == MIN)
				relative[i][l] += static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima() / CSalpchain::globalBest->fitness;
			else
				relative[i][l] += CSalpchain::globalBest->fitness / static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima();

			for (int k = 0; k < CSalpchain::swarmCount; k++) {

				if (0 == l) {
					CSalp::c1 = 2 * exp(-(4 * l / max_gen) ^ 2); //*From SSA Equation No. 3.2
					CSalpchain::subSwarm[k].BuildSalp_Chain(static_cast<Real_DBG*>(g_dbg)->Get_Global_Optima(), fit_evas, r_value, salpParticleArchive, nArchiveIndexNum, debugInfo, traceLog);
				}
				else{
					CSalp::bBeta = (0.5 * (max_gen - l)) / (l + 0.5); //*From QSSA Equation No. 10
					CSalp::accMinionParam = (0.75 * sin( PI/4 * (1 - (l / max_gen)) * Global::uniform.Next()));
					CSalpchain::subSwarm[k].SalpSwarm_Locomate(static_cast<Real_DBG*>(g_dbg)->Get_Global_Optima(), fit_evas, r_value, salpParticleArchive, nArchiveIndexNum, debugInfo, traceLog);
				}
				CSalpchain::subSwarm[k].personalBestUpdate(traceLog);
			}
			traceLog << "Gen : " << Global::gen_Number << " fit[i][l] : " << fit[i][l] << "	Global Optima : " << static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima() << "	Global Best :" << CSalpchain::globalBest->fitness;
			traceLog << endl;
			Global::gen_Number++;
		}
		//traceLog << "Generation :" << Global::gen_Number << "	Global Optima : " << static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima() << "	Global Best :" << CSalpchain::globalBest->fitness << "	Food Position :" << CSalpchain::subSwarm[0].bestPosition.fitness << endl;
		if (!Global::isLogEnable) {
			traceLog << "Population in Change : " << numChange - 1 << endl;
			printPopulation(traceLog);
		}

		if (Global::optimization_type == MIN)
			best_rel_fit[num_run][i] = static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima() / CSalpchain::globalBest->fitness;
		else
			best_rel_fit[num_run][i] = CSalpchain::globalBest->fitness / static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima();

		best_abs_fit[num_run][i] = fabs(CSalpchain::globalBest->fitness - static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima());
		best_rel_fit[num_run][i] = best_rel_fit[num_run][i] / (1 + r_value / (Global::change_frequency*Global::num_dim / Global::sample_frequency));

		g_dbg->Change();

		if (g_dbg->Get_Dimension_Change_Flag() == true) {
			
			g_dbg->Dimension_Change();
			
			if (Composition_DBG * p = dynamic_cast<Composition_DBG*>(g_dbg)) {
				p->Delete_Composition_DBG();
			}
			else if (Rotation_DBG * p = dynamic_cast<Rotation_DBG*>(g_dbg)) {
				p->Delete_Rotation_DBG();
			}

			g_dbg = Global::g_dbg;
		}
		
		if (g_dbg->Get_Dimension_Change_Flag() == false) {

			salpParticleArchive[nArchiveIndexNum++] = CSalpchain::subSwarm[0].bestPosition;

			for (int k = 0; k < CSalpchain::swarmCount; k++) 
				CSalpchain::subSwarm[k].ReinitializeChain(salpParticleArchive, nArchiveIndexNum, traceLog);
			Update_Global_Best(traceLog);
		}
		else {
			Delete_SalpSwarm(traceLog);
			SalpSwarmAlgo_Initial(traceLog);
			Generate_SalpSwarm(traceLog);
		}

		traceLog << endl << "Final Popolation in Change :" << i + 1 << endl;
		printPopulation(traceLog);
		iChangeCounter++;
	}
	if (Global::isLogEnable) {
		traceLog << endl;
		traceLog << "Population best after Generation :" << Global::gen_Number << " and Change Type :" << iChangeCounter << endl;
		printPopulationBest(traceLog);
		traceLog << "Generation :" << Global::gen_Number << "	Global Optima : " << static_cast<Real_DBG *>(g_dbg)->Get_Global_Optima() << "	Global Best :" << CSalpchain::globalBest->fitness << "	Food Position :" << CSalpchain::subSwarm[0].bestPosition.fitness << endl;
	}

	Delete_SalpSwarm(traceLog);
	if (Global::isLogEnable)
		traceLog << "End ->Main::SalpSwarmAlgo()" << endl;
}

double marking(const int f, const int t) {
	double mark;
	if (f == 0) {
		if (t == 6) mark = 0.2*0.1*0.5;
		else mark = 0.2*0.15*0.5;
	}
	else {
		if (t == 6) mark = 0.16*0.1;
		else mark = 0.16*0.15;
	}
	return mark;
}
void Output_result(char * file, double ** fit, double ** relative, double **best_rel_fit, double **best_abs_fit, int *number_dimension, ofstream & perfor, const int f, const int t, double &total_mark) {

	char name[50], chName[50];
	strcpy(chName, "SSA-TracelogInfo\\");
	strcpy(name, file);
	strcat(name, "_fit.txt");
	strcat(chName, name);
	ofstream ofit(chName);

	strcpy(chName, "SSA-TracelogInfo\\");
	strcpy(name, file);
	strcat(name, "_relative.txt");
	strcat(chName, name);
	ofstream orel(chName);

	strcpy(chName, "SSA-TracelogInfo\\");
	strcpy(name, file);
	strcat(name, "_statistic.txt");
	strcat(chName, name);
	ofstream osta(chName);

	for (int i = 0; i < Global::num_change; i++) {
		for (int j = 0; j < Global::change_frequency*number_dimension[i] / Global::max_popsize; j++) {
			ofit << fit[i][j] / Num_Run << endl;
			orel << relative[i][j] / Num_Run << endl;
		}
	}
	double min, max, std = 0, avg_rel = 0, avg_abs = 0;
	double avg_min = 0, avg_max = 0;
	for (int i = 0; i < Num_Run; i++) {
		min = max = best_abs_fit[i][0];
		for (int j = 0; j < Global::num_change; j++) {
			if (min > best_abs_fit[i][j])
				min = best_abs_fit[i][j];
			if (max < best_abs_fit[i][j])
				max = best_abs_fit[i][j];
			avg_rel += best_rel_fit[i][j];
			avg_abs += best_abs_fit[i][j];
		}
		avg_min += min;
		avg_max += max;
	}
	avg_min /= Num_Run;
	avg_max /= Num_Run;
	avg_abs = avg_abs / (Num_Run*Global::num_change);
	avg_rel = avg_rel / (Num_Run*Global::num_change);

	for (int i = 0; i < Num_Run; i++) {
		for (int j = 0; j < Global::num_change; j++)
			std += (best_abs_fit[i][j] - avg_abs)*(best_abs_fit[i][j] - avg_abs);
	}
	std = sqrt(std / (Num_Run*Global::num_change - 1));
	osta << "avg_min: " << avg_min << endl;
	osta << "avg_max: " << avg_max << endl;
	osta << "avg_mean: " << avg_abs << endl;
	osta << "std: " << std << endl;

	double node_mark = marking(f, t);
	perfor <<"	" << f + 1 << "	|	" << t + 1 << "		  |	" << avg_abs << "		|	 " << avg_rel << "	|	 " << node_mark << "	|	" << node_mark * avg_rel << endl;
	total_mark += node_mark * avg_rel;

	ofit.close();
	orel.close();
	osta.close();


	for (int i = 0; i < Global::num_change; i++) {
		delete[]fit[i];
		delete[]relative[i];
	}
	delete[] fit;
	delete[] relative;
}
void Generate_file_name(char * file, const int f, const int t, const int d, const int m) {
	// the name returned by file, f: function, t: change type, d: dimension, m: the number of peaks or basci functions
	char t1[50], t2[50];
	sprintf(t1, "%d", m);
	strcpy(t2, "_peak");
	strcat(t2, t1);
	sprintf(t1, "%d", d);
	strcat(t1, t2);
	strcpy(t2, "_D");
	strcat(t2, t1);
	sprintf(t1, "%d", t + 1);
	strcat(t1, t2);
	strcpy(t2, "_T");
	strcat(t2, t1);
	sprintf(t1, "%d", f + 1);
	strcat(t1, t2);
	strcpy(t2, "F");
	strcat(t2, t1);
	strcpy(file, t2);
}
void run(const int f, const int t, const int num_dim, const int num_peak, double *seed, ofstream &perfor, double & total_mark, ofstream &debugInfo, ofstream &traceLog) {
	//run(function, change_type, number_dimensions, number_peaks, seed, output file, performance mark)
	int *number_dimension = new int[Global::num_change]; // save the number of dimensions of each change
	double **relative = new double*[Global::num_change];
	double **fit = new double*[Global::num_change];
	double **best_rel_fit = new double *[Num_Run];
	double **best_abs_fit = new double *[Num_Run];
	for (int i = 0; i < Num_Run; i++) {
		best_rel_fit[i] = new double[Global::num_change];
		best_abs_fit[i] = new double[Global::num_change];
	}
	
	if (Global::isLogEnable)
		traceLog << "Start-> Main::Run()" <<endl;

	if (f == 0) {
		for (int i = 0; i < Num_Run; i++) {
			debugInfo << endl << "Function : " << f + 1 << " Change : " << t + 1<< endl;
			
			if (Global::isLogEnable)
				traceLog << "Function : " << f+1 << " Change : " << t+1 << " Num-Dim : " << num_dim << " Num_Run : " << i+1 <<endl;
			
			System_Initial(t, seed[i], num_dim, traceLog);
			Global::optimization_type = MAX;
			Rotation_DBG *rot_dbg;
			Global::num_peakorfun = num_peak;
			Global::max_width = 10.;
			Global::min_width = 1;
			Global::min_height = 10;
			Global::max_height = 100;
			rot_dbg = Rotation_DBG::Get_Rotation_DBG();
			System_Setting(rot_dbg, f);
			Global::g_dbg = rot_dbg;
			SalpSwarmAlgo(best_rel_fit, best_abs_fit, fit, relative, rot_dbg, number_dimension, i, debugInfo, traceLog); // call Salp Swarm algorithm
			rot_dbg->Delete_Rotation_DBG();
		}
		char file[50];
		Generate_file_name(file, f, t, num_dim, Global::num_peakorfun); // generate output file name
		Output_result(file, fit, relative, best_rel_fit, best_abs_fit, number_dimension, perfor, f, t, total_mark);

	}
	else {
		debugInfo << endl << "Function : " << f + 1 << " Change : " << t + 1 << endl;
		for (int i = 0; i < Num_Run; i++) {
			
			if (Global::isLogEnable)
				traceLog << "Function : " << f+1 << " Change : " << t+1 << " Num-Dim : " << num_dim << " Num_Run : " << i << endl;
			
			System_Initial(t, seed[i], num_dim, traceLog);
			Global::optimization_type = MIN;
			Composition_DBG * com_dbg;
			Global::num_peakorfun = 10;
			Global::hn_s = 2000.;
			Global::min_height = 10;
			Global::max_height = 100;

			com_dbg = Composition_DBG::Get_Composition_DBG();
			System_Setting(com_dbg, f);
			Global::g_dbg = com_dbg;

			SalpSwarmAlgo(best_rel_fit, best_abs_fit, fit, relative, com_dbg, number_dimension, i, debugInfo, traceLog); // call Salp Swarm algorithm
			com_dbg->Delete_Composition_DBG();
		}
		char file[50];
		Generate_file_name(file, f, t, num_dim, Global::num_peakorfun);// generate output file name
		Output_result(file, fit, relative, best_rel_fit, best_abs_fit, number_dimension, perfor, f, t, total_mark);
	}
	for (int i = 0; i < Num_Run; i++) {
		delete[] best_rel_fit[i];
		delete[] best_abs_fit[i];
	}

	delete[] best_rel_fit;
	delete[] best_abs_fit;

	delete[] number_dimension;
	
	if (Global::isLogEnable)
		traceLog << "End-> Main::Run()" << endl;
}
int main() {
	srand(17);
	double seed[Num_Run];
	for (int i = 0; i < Num_Run; i++)
		seed[i] = (double)rand() / RAND_MAX;
	int num_peak;

	ofstream perfor("SSA-TracelogInfo\\performance.txt");
	double total_mark = 0;
	perfor << "Problem : Change_type :	Avg_abs		:	Avg_rel		: Case_mark :	Score	" << endl;

	ofstream ofDebugInfo("SSA-DebuggingInfo\\PerformanceLog.txt");
	ofstream ofTraceLog("SSA-DebuggingInfo\\TraceLog.txt");
	ofDebugInfo << "Generation | Parameter(s1)| CalcFit(F)      | CalcVal(R)    | Global Optima | FoodPosition | GlobalBest" << endl;

	if(Global::isLogEnable)
		ofTraceLog << "Start -> Main::main()"<<endl;

	for (int f = 0; f < 6; f++) {// 6 test problems
		for (int t = 0; t < (Num_Change_Type + 1); t++) {// 7 (Num_Change_Type + 1)change types
			if (f == 0) {
				for (int p = 0; p < 2; p++) {// 2 test number of peak: 10 and 50
					switch (p) {
					case 0:
						num_peak = 10;
						break;
					case 1:
						num_peak = 50;
						break;
					}
					run(f, t, 10, num_peak, seed, perfor, total_mark, ofDebugInfo, ofTraceLog);
				}
			}
			else	run(f, t, 10, 10, seed, perfor, total_mark, ofDebugInfo, ofTraceLog);

		}
	}
	delete[] CSalp::vmax;
	perfor << "	Total mark (100*sum(score)):	" << 100 * total_mark << endl;
	
	if (Global::isLogEnable)
		ofTraceLog << "End -> Main::main()" << endl;
	
	perfor.close();
	ofDebugInfo.close();
	ofTraceLog.close();
	return 1;
}
