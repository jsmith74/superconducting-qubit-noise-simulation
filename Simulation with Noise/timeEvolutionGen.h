#ifndef TIMEEVOLUTIONGEN_H_INCLUDED
#define TIMEEVOLUTIONGEN_H_INCLUDED

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <mpi.h>
#include <string>
#include <sstream>
#include<iomanip>


class importH{

public:

	importH(int samplingPoints,int numberoffluxparams,int TEorder,int truncHSD);
	Eigen::MatrixXcd at(Eigen::MatrixXcd& Htot,Eigen::VectorXd& noiseLocation);
	void set(Eigen::MatrixXcd H[]);

private:

	std::complex<double> I;
	Eigen::MatrixXi partialDerivativeIndexMatrix;
	int numbPoints,numbFluxParameters,subHSDimension,taylorExpansionOrder,numbPartialDerivMatrices,numbPrevMPIProcesses;
	int genNumbPartialDerivMatrices(int& taylorExpansionOrder,int& numbFluxParameters);
	double doublefactorial(int x);
	inline int g(int n,int m);
	Eigen::MatrixXi generateSubBasisVector(int subPhotons, int subModes);
	void setThisH(Eigen::MatrixXcd& H,std::ifstream& HimportStream);
	int genNumPrevMPIProcesses();
	void setFilename(std::string& filename,int fileNumber);
	double taylorArgument(Eigen::VectorXd& testNoiseVector,Eigen::VectorXi partialVector);
	Eigen::VectorXi numbPartialDerivMatricesAt;
	void fixMPIPhaseDiscontinuity(Eigen::MatrixXcd& H,Eigen::VectorXcd& beginPhaseDiscontinuity);
};

class timeEvolutionOperatorGen{

    public:
        void generate(Eigen::MatrixXcd Htot[],importH& taylorH);
        timeEvolutionOperatorGen(double xStart,double xEnd,int numSamplingPoints,int truncHSD,int numberofFluxP);
        double gateFidelity(Eigen::MatrixXcd& U1,Eigen::MatrixXcd& U2);
        void printTEOperatorRect();
        void printTEOperatorPolar();

    private:
        int subHSDimension,numbPoints,numbFluxParameters,numbHMatrices,numbPartialDerivMatrices;
        double tInit,tFin,hb,dt,t;
        std::complex<double> I;
        Eigen::MatrixXcd matrixExp(Eigen::MatrixXcd X);
        Eigen::VectorXd fluxNoise;
        Eigen::MatrixXcd Uop;
        void updateFluxNoise(double& t);

};


#endif // TIMEEVOLUTIONGEN_H_INCLUDED
