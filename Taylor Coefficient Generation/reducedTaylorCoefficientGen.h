/*! \file reducedTaylorCoefficientGen.h
 *  \brief header file for reducedTaylorCoefficientGen.cpp
 *
 *
 *      This is the class and function definitions for reducedTaylorCoefficientGen.cpp
 *
 *
 * \author Jake Smith <jsmith74@tulane.edu>
 */



#ifndef NONADIABHGEN_H_INCLUDED
#define NONADIABHGEN_H_INCLUDED

#include "finiteDifferenceFuncs.h"





class RTCGen{

	public:

		RTCGen(double xStart,double xEnd,int numSamplingPoints,int initialHSD,int truncHSD,int NderivPoints,int NFluxParams,int TEorder,int RdagDerivApprox,double FDI,int procNumb);
		void generate(Eigen::MatrixXcd& Htot);

	private:

        Eigen::VectorXd fluxVec;
        void updateFluxVec(double& x);
        double OPTIMAL_DERIVATIVE_H;
		int numbPartialDerivMatrices,numbFluxParameters,numbPoints, pointsPerProcess, pointsLeftover, beginandendbuffer,numbRdagDerivativePoints, processNumber;
		int taylorExpansionOrder,numbHMatrices,startIndex,endIndex,myid,numprocs, numbIndices, subHSDimension,totalHSDimension,numbDerivativePoints,derivGridMidpoint;
		double dx,x,xFinal,xInit,dtTaylorProj;
		std::complex<double> I;
		double normCompare(Eigen::VectorXcd vec1, Eigen::VectorXcd vec2);

		void initializeNDSamplingInterval(Eigen::MatrixXd eigenValues[],Eigen::MatrixXcd Rdag[],Eigen::MatrixXcd& RdagPrime,Eigen::MatrixXd& eigenValuesProj,Eigen::MatrixXi& eVDegenAddress,std::string filename);
		void mainGeneratorLoop(Eigen::MatrixXcd& Htot,Eigen::MatrixXd eigenValues[],Eigen::MatrixXcd Rdag[],Eigen::MatrixXcd& RdagPrime);
		inline void setRdagPrime(Eigen::MatrixXcd Rdag[],Eigen::MatrixXcd& RdagPrime,int& oldestIndex,finiteDifferenceFormula& numericalDerivativeRdag);
		inline void addEigenvaluestoDiagonals(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Htot);
		Eigen::VectorXd zeroVec;
        void setFluxPerturbatedEVandRdagFirstOrderInit(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress);
		void setFluxPerturbatedEVandRdagFirstOrder(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress);
		inline void setHMatrices(Eigen::MatrixXcd& Htot,Eigen::MatrixXcd& Rdag,Eigen::MatrixXcd& RdagPrime,Eigen::MatrixXd& eigenValues,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXi partialDerivativeIndexMatrix[]);
		void correctInitialEigenbasisPhase(Eigen::MatrixXcd& Rdag);
		int genNumbPartialDerivMatrices(int& taylorExpansionOrder,int& numbFluxParameters);
		double doublefactorial(int x);
		void setFluxPerturbatedEVandRdagHigherOrderInit(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filename,Eigen::MatrixXi& eVDegenAddress);
		void setFluxPerturbatedEVandRdagHigherOrder(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress);
		int HIConvert(Eigen::VectorXi& gridPosition);
		void iterate(Eigen::VectorXi& gridPosition);
		Eigen::MatrixXi generateSubBasisVector(int subPhotons, int subModes);
		inline int g(int n,int m);
		void buildHtotFirstOrder(Eigen::MatrixXcd& Htot,Eigen::MatrixXcd& HtotTemp,finiteDifferenceFormula& numericalDerivative);
		void buildHtotHigherOrder(Eigen::MatrixXcd& Htot,Eigen::MatrixXcd& HtotTemp,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXi partialDerivativeIndexMatrix[]);
		inline void takePartialDerivative(Eigen::MatrixXcd& Htot,Eigen::MatrixXcd& HtotTemp,Eigen::VectorXi v,int& k,finiteDifferenceFormula numericalDerivative[]);
		inline void collapseDim(int& dim,int& derivOrder,Eigen::VectorXi& gridPosition,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXcd& reducedHyperCube,int& loc,Eigen::MatrixXcd& HtotTemp);
		inline void setReducedHyperCubeInitial(Eigen::VectorXi& v,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXcd& reducedHyperCube,Eigen::MatrixXcd& HtotTemp);
		inline void setReducedHyperCube(Eigen::VectorXi& v,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXcd& reducedHyperCubeIn,Eigen::MatrixXcd& reducedHyperCubeOut);
        void pullEigensystemfromMathematica(int placingAddress,Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,Eigen::MatrixXi& eVDegenAddress);
        void printMatrixToMathematica(Eigen::MatrixXcd M);
        void outputComplexMtoFile(Eigen::MatrixXcd& M,std::ofstream& outfile);
        void setEigenValuesProj(Eigen::MatrixXd eigenValues[],Eigen::MatrixXd& eigenValuesProj,int& oldestIndex,finiteDifferenceFormula numericalDerivativeEv[]);
        void seteVDegenAddress(Eigen::MatrixXd& eigenValuesProj,Eigen::MatrixXi& eVDegenAddress);
        bool isCorrectOrder(Eigen::MatrixXd& eigenValuesProj,Eigen::MatrixXi& eVDegenAddress);
        void initializeEVDegenAddress(Eigen::MatrixXi& eVDegenAddress);
        void checkDegenEvPerturbFirstOrder(Eigen::MatrixXcd& Rdag,Eigen::MatrixXd& eigenValues,std::string& filename,Eigen::MatrixXi& eVDegenAddress);
        void correctInitialEigenDegeneracy(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,Eigen::MatrixXi& eVDegenAddress);
        void popVecEnt(Eigen::VectorXi& v,int i);
        void checkDegenEvPerturbHigherOrder(Eigen::MatrixXcd& Rdag,Eigen::MatrixXd& eigenValues,std::string& filename,Eigen::MatrixXi& eVDegenAddress);
        void iterateDegenEv(Eigen::VectorXi& gridPosition,bool dir[]);
        void correctInitialNDSamplingDegeneracy(Eigen::MatrixXcd& RdagPrev,Eigen::MatrixXcd& RdagCurr,Eigen::MatrixXd& eigenValues,Eigen::MatrixXi& eVDegenAddress);
        void removeOldFiles(std::string& filename, std::string filenameEv);
        void printFinalBasis(Eigen::MatrixXcd& Rdag);
        void printMathematicaForm(Eigen::MatrixXcd M,std::ofstream& outfile);

        void setFluxPerturbatedEVandRdagZerothOrderInit(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress);
        void setFluxPerturbatedEVandRdagZerothOrder(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress);

};


#endif
