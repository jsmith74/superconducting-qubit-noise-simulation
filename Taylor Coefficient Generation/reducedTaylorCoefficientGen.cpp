/*! \file reducedTaylorCoefficientGen.cpp
 *  \brief functions used to generate taylor series coefficients
 *
 *       These are the functions for importing a large Hamiltonian from Mathematica, converting it into the
 *      instantaneous eigenbasis, and then generating the taylor expansion coefficients to allow later consideration of flux noise See the user guide.
 *
 *      In order to change the output format of the taylor coefficients, change the first function commented below.
 *
 * \author Jake Smith <jsmith74@tulane.edu>
 */



#include "reducedTaylorCoefficientGen.h"

#include <iomanip>
#include <fstream>

#define PTOL 1.0e-10


//#define DEBUG_ENABLED


/** ====== THIS FUNCTION CONTROLS THE OUTPUT OF THE TAYLOR COEFFICIENTS. TO CHANGE THE WAY THESE MATRICES ARE EXPORTED, YOU ONLY NEED TO EDIT THIS FUNCTION ========================================================================================================== */



void RTCGen::outputComplexMtoFile(Eigen::MatrixXcd& M,std::ofstream& outfile){

    double current_Interval_Midpoint_Time = x-(beginandendbuffer+1)*dx;

    outfile << std::setprecision(16) << current_Interval_Midpoint_Time << std::endl << std::endl;

    for(int i=0;i<numbPartialDerivMatrices;i++){

       outfile << M.block(0,i*subHSDimension,subHSDimension,subHSDimension) << std::endl << std::endl;
       //printMathematicaForm(M.block(0,i*subHSDimension,subHSDimension,subHSDimension),outfile);

    }

    outfile << std::endl << std::endl << std::endl << std::endl;

    return;

}



/** ================================================================================================================================================================================================================================================================= */




void RTCGen::generate(Eigen::MatrixXcd& Htot){

	Eigen::MatrixXcd Rdag[numbRdagDerivativePoints];

	Eigen::MatrixXd eigenValues[numbRdagDerivativePoints];

	Eigen::MatrixXcd RdagPrime(totalHSDimension,subHSDimension*numbHMatrices);

	mainGeneratorLoop(Htot,eigenValues,Rdag,RdagPrime);

	return;

}




void RTCGen::mainGeneratorLoop(Eigen::MatrixXcd& Htot,Eigen::MatrixXd eigenValues[],Eigen::MatrixXcd Rdag[],Eigen::MatrixXcd& RdagPrime){

    std::stringstream ss;
    ss << processNumber;
    std::string procNumb;
    ss >> procNumb;
    std::string filename,filenameEv;
    filename =  "RTC.dat";//"H" + procNumb + ".dat";
    filenameEv = "evH" + procNumb + ".dat";

    removeOldFiles(filename,filenameEv);

	Eigen::MatrixXd eigenValuesProj(subHSDimension,numbHMatrices);

    dtTaylorProj = dx * ((numbRdagDerivativePoints-1)/2 +1);

    Eigen::MatrixXi eVDegenAddress(subHSDimension,numbHMatrices);

    initializeEVDegenAddress(eVDegenAddress);


	initializeNDSamplingInterval(eigenValues,Rdag,RdagPrime,eigenValuesProj,eVDegenAddress,filenameEv);


	finiteDifferenceFormula numericalDerivative[taylorExpansionOrder];
	for(int i=0;i<taylorExpansionOrder;i++)	numericalDerivative[i].setfiniteDifferenceFormula(numbDerivativePoints,i+1,OPTIMAL_DERIVATIVE_H);

	Eigen::MatrixXi partialDerivativeIndexMatrix[taylorExpansionOrder];
	for(int i=0;i<taylorExpansionOrder;i++) partialDerivativeIndexMatrix[i] = generateSubBasisVector(i+1,numbFluxParameters);

	std::ofstream PDIMOut("Taylor_Coefficient_Order.dat");
    for(int i=0;i<=taylorExpansionOrder;i++) PDIMOut << generateSubBasisVector(i,numbFluxParameters) << std::endl << std::endl;
	PDIMOut.close();

	finiteDifferenceFormula numericalDerivativeRdag;

	numericalDerivativeRdag.setfiniteDifferenceFormula(numbRdagDerivativePoints,1,dx);

	finiteDifferenceFormula numericalDerivativeEv[numbRdagDerivativePoints-1];

	for(int i=0;i<numbRdagDerivativePoints-1;i++) numericalDerivativeEv[i].setfiniteDifferenceFormula(numbRdagDerivativePoints,i+1,dx);

	Htot.resize(subHSDimension,numbPartialDerivMatrices*subHSDimension);


	for(int i=0;i<numbIndices;i++){


		int oldestIndex = i % numbRdagDerivativePoints;
		int centralIndex = (oldestIndex + beginandendbuffer) % numbRdagDerivativePoints;
		int prevIndex = oldestIndex - 1;
		if(prevIndex < 0) prevIndex = numbRdagDerivativePoints - 1;

		setRdagPrime(Rdag,RdagPrime,oldestIndex,numericalDerivativeRdag);

        setEigenValuesProj(eigenValues,eigenValuesProj,oldestIndex,numericalDerivativeEv);

        seteVDegenAddress(eigenValuesProj,eVDegenAddress);

		setHMatrices(Htot,Rdag[centralIndex],RdagPrime,eigenValues[centralIndex],numericalDerivative,partialDerivativeIndexMatrix);

        std::ofstream Hout(filename.c_str(),std::ofstream::app);
        outputComplexMtoFile(Htot,Hout);
        Hout.close();

        if(taylorExpansionOrder==0 && i<numbIndices-1) setFluxPerturbatedEVandRdagZerothOrder(eigenValues[oldestIndex],Rdag[oldestIndex],filenameEv,eVDegenAddress);
		else if(taylorExpansionOrder == 1 && i<numbIndices-1) setFluxPerturbatedEVandRdagFirstOrder(eigenValues[oldestIndex],Rdag[oldestIndex],filenameEv,eVDegenAddress);
		else if(i<numbIndices-1) setFluxPerturbatedEVandRdagHigherOrder(eigenValues[oldestIndex],Rdag[oldestIndex],filenameEv,eVDegenAddress);


		for(int j=0;j<subHSDimension*numbHMatrices;j++){

            double phaseArg = std::arg(Rdag[oldestIndex](0,j));
            Rdag[oldestIndex].col(j) *= exp(-I*phaseArg);

		}

		x += dx;

		if(i==numbIndices-1) printFinalBasis(Rdag[centralIndex]);

	}


}




void RTCGen::printMathematicaForm(Eigen::MatrixXcd M,std::ofstream& outfile){

	outfile << "{";

	for(int i=0;i < M.rows();i++){

		outfile << "{";

		for(int j=0;j<M.cols();j++){

			outfile << std::setprecision(16) << std::real(M(i,j));
			if(std::imag(M(i,j))<0) outfile << std::imag(M(i,j)) << "*I";
			else outfile << "+" << std::imag(M(i,j)) << "*I";
			if(j<M.cols()-1) outfile << ",";

		}

		outfile << "}";
		if(i<M.rows()-1) outfile << ",";

	}

	outfile << "}";

	outfile << std::endl << std::endl;

	return;

}




void RTCGen::initializeEVDegenAddress(Eigen::MatrixXi& eVDegenAddress){


    for(int i=0;i<subHSDimension;i++){

        for(int j=0;j<numbHMatrices;j++) eVDegenAddress(i,j) = i;

    }

    return;

}




void RTCGen::updateFluxVec(double& x){

    WSPutFunction(stdlink,"EvaluatePacket",1);
    WSPutFunction(stdlink,"fluxVec",1);
    WSPutReal64(stdlink,x);
    WSEndPacket(stdlink);
    long size;
    WSCheckFunction(stdlink,"ReturnPacket",&size);
    WSCheckFunction(stdlink,"List",&size);
    for(int i=0;i<numbFluxParameters;i++) WSGetReal64(stdlink,&fluxVec(i));
    return;

}




RTCGen::RTCGen(double xStart,double xEnd,int numSamplingPoints,int initialHSD,int truncHSD,int NderivPoints,int NFluxParams,int TEorder,int RdagDerivApprox,double FDI,int procNumb){


    processNumber = procNumb;
    OPTIMAL_DERIVATIVE_H = FDI;
	taylorExpansionOrder = TEorder;
	numbFluxParameters = NFluxParams;

	numbPartialDerivMatrices = genNumbPartialDerivMatrices(taylorExpansionOrder,numbFluxParameters);

	numbHMatrices = std::pow(NderivPoints,numbFluxParameters);

	if(taylorExpansionOrder == 0) numbHMatrices = 1;
	else if(taylorExpansionOrder == 1) numbHMatrices = (NderivPoints-1) * numbFluxParameters + 1;
	else numbHMatrices = std::pow(NderivPoints,numbFluxParameters);

	zeroVec = Eigen::VectorXd::Zero(numbFluxParameters);

	numbPoints = numSamplingPoints;
	subHSDimension = truncHSD;


	numbDerivativePoints = NderivPoints;
	numbRdagDerivativePoints = RdagDerivApprox;

	derivGridMidpoint = (numbDerivativePoints-1)/2;

	myid =0;
	numprocs = 1;


	pointsPerProcess = numbPoints/numprocs;
	pointsLeftover = numbPoints%numprocs;
	beginandendbuffer = (numbRdagDerivativePoints-1)/2;

	if(myid<pointsLeftover){
		startIndex = myid*(pointsPerProcess+1);
		endIndex = startIndex + pointsPerProcess;
	}

	else{
		startIndex = pointsLeftover * (pointsPerProcess + 1) + (myid-pointsLeftover) * pointsPerProcess;
		endIndex = startIndex + pointsPerProcess - 1;
	}

	numbIndices = endIndex - startIndex + 1;

	xFinal = xEnd;
	xInit = xStart;

	dx  = (xFinal-xInit)/(numbPoints-1);

	printf("Grid size: %f \n",dx);

	x = xInit + startIndex*dx - beginandendbuffer*dx;

	std::complex<double> Igen(0.0,1.0);
	I = Igen;

	totalHSDimension = initialHSD;

    fluxVec.resize(numbFluxParameters);

}




double RTCGen::normCompare(Eigen::VectorXcd vec1, Eigen::VectorXcd vec2){
    std::complex<double> output = vec1.conjugate().transpose() * vec2;
    return std::norm(output);

}



void RTCGen::pullEigensystemfromMathematica(int placingAddress,Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,Eigen::MatrixXi& eVDegenAddress){

    long size;


    WSPutFunction(stdlink,"EvaluatePacket",1);
    WSPutFunction(stdlink,"es",1);
    WSPutFunction(stdlink,"List",numbFluxParameters);
    for(int i=0;i<numbFluxParameters;i++){
        WSPutReal64(stdlink,fluxVec(i));
    }
    WSEndPacket(stdlink);


    WSCheckFunction(stdlink,"ReturnPacket",&size);

    WSCheckFunction(stdlink,"List",&size);

    WSCheckFunction(stdlink,"List",&size);

    for(int i=0;i<subHSDimension;i++){
        if(!WSCheckFunction(stdlink,"Complex",&size)){
            WSClearError(stdlink);
            WSGetReal64(stdlink,&eigenValues(eVDegenAddress(i,placingAddress),placingAddress));
            continue;
        }
        double imagThrowAway;
        WSGetReal64(stdlink,&eigenValues(eVDegenAddress(i,placingAddress),placingAddress));
        WSGetReal64(stdlink,&imagThrowAway);
    }


    WSCheckFunction(stdlink,"List",&size);

    for(int i=0;i<subHSDimension;i++){
        WSCheckFunction(stdlink,"List",&size);

        for(int j=0;j<totalHSDimension;j++){
            if(!WSCheckFunction(stdlink,"Complex",&size)){
                WSClearError(stdlink);
                double realPart;
                WSGetReal64(stdlink,&realPart);
                Rdag(j,eVDegenAddress(i,placingAddress)+placingAddress*subHSDimension) = realPart;
                continue;
            }

            double realPart,imagPart;
            WSGetReal64(stdlink,&realPart);
            WSGetReal64(stdlink,&imagPart);
            Rdag(j,eVDegenAddress(i,placingAddress)+placingAddress*subHSDimension) = realPart + I * imagPart;
        }
    }

    return;
}




void RTCGen::correctInitialEigenbasisPhase(Eigen::MatrixXcd& Rdag){


    for(int i=0;i<numbHMatrices*subHSDimension;i++){

        double phaseArg = std::arg(Rdag(0,i));
        Rdag.col(i) *= exp(-I*phaseArg);

    }


	return;

}




int RTCGen::genNumbPartialDerivMatrices(int& taylorExpansionOrder,int& numbFluxParameters){

	int output = 0;
	for(int i=0;i<=taylorExpansionOrder;i++){

		output += doublefactorial(i+numbFluxParameters-1)/(doublefactorial(i)*doublefactorial(numbFluxParameters-1));

	}

	return output;

}




void RTCGen::setEigenValuesProj(Eigen::MatrixXd eigenValues[],Eigen::MatrixXd& eigenValuesProj,int& oldestIndex,finiteDifferenceFormula numericalDerivativeEv[]){

	Eigen::VectorXi ix(numbRdagDerivativePoints);
	for(int i=0;i<numbRdagDerivativePoints;i++){
		ix(i) = (oldestIndex+i) % numbRdagDerivativePoints;
	}

	int centralIndex = (numbRdagDerivativePoints - 1) / 2;

    eigenValuesProj = eigenValues[ix(centralIndex)];


    for(int i=0;i<numbRdagDerivativePoints-1;i++){

        for(int j=0;j<numbRdagDerivativePoints;j++){

            eigenValuesProj += (1.0/doublefactorial(i+1)) * std::pow(dtTaylorProj,i+1) * eigenValues[ix(j)] * numericalDerivativeEv[i].weights(j);

        }

    }

#ifdef DEBUG_ENABLED

    for(int i=0;i<numbHMatrices;i++){

        std::stringstream ss;
        ss << i;
        std::string s;
        ss >> s;

        std::string filenameProj = "eV"  + s + "Proj.dat";

        std::ofstream eigenvalueCheck(filenameProj.c_str(),std::ofstream::app);
        eigenvalueCheck << std::setprecision(16) << x << "\t";
        for(int ii=0;ii<subHSDimension;ii++) eigenvalueCheck << eigenValuesProj(ii,i) << "\t";
        eigenvalueCheck << std::endl;
        eigenvalueCheck.close();

    }

#endif

    return;

}




inline void RTCGen::setRdagPrime(Eigen::MatrixXcd Rdag[],Eigen::MatrixXcd& RdagPrime,int& oldestIndex,finiteDifferenceFormula& numericalDerivativeRdag){

	Eigen::VectorXi ix(numbRdagDerivativePoints);
	for(int i=0;i<numbRdagDerivativePoints;i++){
		ix(i) = (oldestIndex+i) % numbRdagDerivativePoints;
	}

	RdagPrime = numericalDerivativeRdag.weights(0) * Rdag[ix(0)];

	for(int i=1;i<numbRdagDerivativePoints;i++){
		RdagPrime += numericalDerivativeRdag.weights(i) * Rdag[ix(i)];
	}

	return;
}




bool RTCGen::isCorrectOrder(Eigen::MatrixXd& eigenValuesProj,Eigen::MatrixXi& eVDegenAddress){

    for(int j=0;j<numbHMatrices;j++){

        for(int i=0;i<subHSDimension-1;i++){

            if(eigenValuesProj(eVDegenAddress(i,j),j) > eigenValuesProj(eVDegenAddress(i+1,j),j)){

                int tempSwap = eVDegenAddress(i,j);
                eVDegenAddress(i,j) = eVDegenAddress(i+1,j);
                eVDegenAddress(i+1,j) = tempSwap;

                return false;
            }

        }
    }

    return true;

}




void RTCGen::seteVDegenAddress(Eigen::MatrixXd& eigenValuesProj,Eigen::MatrixXi& eVDegenAddress){

    while(true){

        if(isCorrectOrder(eigenValuesProj,eVDegenAddress)) break;

    }

    return;

}



void RTCGen::setFluxPerturbatedEVandRdagZerothOrder(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress){

    updateFluxVec(x);

    pullEigensystemfromMathematica(0,eigenValues,Rdag,eVDegenAddress);

    std::ofstream eigenvalueCheck(filenameEv.c_str(),std::ofstream::app);
	eigenvalueCheck << std::setprecision(16) << x << "\t";
	for(int i=0;i<subHSDimension;i++) eigenvalueCheck << eigenValues(i,0) << "\t";
	eigenvalueCheck << std::endl;
	eigenvalueCheck.close();


	return;

}




void RTCGen::setFluxPerturbatedEVandRdagFirstOrder(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress){

    updateFluxVec(x);

    pullEigensystemfromMathematica(0,eigenValues,Rdag,eVDegenAddress);

    std::ofstream eigenvalueCheck(filenameEv.c_str(),std::ofstream::app);
	eigenvalueCheck << std::setprecision(16) << x << "\t";
	for(int i=0;i<subHSDimension;i++) eigenvalueCheck << eigenValues(i,0) << "\t";
	eigenvalueCheck << std::endl;
	eigenvalueCheck.close();

	for(int i=0;i<numbFluxParameters;i++){

		Eigen::VectorXd fluxPerturbations = Eigen::VectorXd::Zero(numbFluxParameters);

		int k = 0;
		for(int j=0;j<numbDerivativePoints-1;j++){

			if(j==derivGridMidpoint) k=1;
			fluxPerturbations(i) = (-1*derivGridMidpoint + j + k)*OPTIMAL_DERIVATIVE_H;

            fluxVec += fluxPerturbations;

            pullEigensystemfromMathematica((numbDerivativePoints-1)*i+1+j,eigenValues,Rdag,eVDegenAddress);

            fluxVec -= fluxPerturbations;
		}

	}


#ifdef DEBUG_ENABLED

	for(int i=0;i<numbHMatrices;i++){

        std::stringstream ss;
        ss << i;
        std::string s;
        ss >> s;
        std::string filenamePert = "eV" + s + "Pert.dat";
        std::ofstream eigenvalueCheckt(filenamePert.c_str(),std::ofstream::app);
        eigenvalueCheckt << std::setprecision(16) << x << "\t";
        for(int k=0;k<subHSDimension;k++) eigenvalueCheckt << eigenValues(k,i) << "\t";
        eigenvalueCheckt << std::endl;

        eigenvalueCheckt.close();

	}

#endif

	return;

}




void RTCGen::setFluxPerturbatedEVandRdagHigherOrder(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress){


    updateFluxVec(x);

	Eigen::VectorXi gridPosition(numbFluxParameters);

	for(int i=0;i<numbFluxParameters;i++) gridPosition(i) = -1*derivGridMidpoint;

	gridPosition(0)--;

	Eigen::VectorXd fluxPerturbations(numbFluxParameters);

	for(int i=0;i<numbHMatrices;i++){

		iterate(gridPosition);

		for(int j=0;j<numbFluxParameters;j++) fluxPerturbations(j) = OPTIMAL_DERIVATIVE_H * gridPosition(j);

        fluxVec += fluxPerturbations;

		pullEigensystemfromMathematica(HIConvert(gridPosition),eigenValues,Rdag,eVDegenAddress);

        if(i==(numbHMatrices-1)/2){

            std::ofstream eigenvalueCheck(filenameEv.c_str(),std::ofstream::app);
            eigenvalueCheck << std::setprecision(16) << x << "\t";
            for(int ii=0;ii<subHSDimension;ii++) eigenvalueCheck << eigenValues(ii,i) << "\t";
            eigenvalueCheck << std::endl;
            eigenvalueCheck.close();

		}


        fluxVec -= fluxPerturbations;

	}

#ifdef DEBUG_ENABLED

	for(int i=0;i<numbHMatrices;i++){

        std::stringstream ss;
        ss << i;
        std::string s;
        ss >> s;
        std::string filenamePert = "eV" + s + "Pert.dat";
        std::ofstream eigenvalueCheckt(filenamePert.c_str(),std::ofstream::app);
        eigenvalueCheckt << std::setprecision(16) << x << "\t";
        for(int k=0;k<subHSDimension;k++) eigenvalueCheckt << eigenValues(k,i) << "\t";
        eigenvalueCheckt << std::endl;

        eigenvalueCheckt.close();

	}

#endif

	return;

}

void RTCGen::setFluxPerturbatedEVandRdagZerothOrderInit(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress){

    updateFluxVec(x);

    pullEigensystemfromMathematica(0,eigenValues,Rdag,eVDegenAddress);

    std::ofstream eigenvalueCheck(filenameEv.c_str(),std::ofstream::app);
	eigenvalueCheck << std::setprecision(16) << x << "\t";
	for(int i=0;i<subHSDimension;i++) eigenvalueCheck << eigenValues(i,0) << "\t";
	eigenvalueCheck << std::endl;
	eigenvalueCheck.close();


	return;

}



void RTCGen::setFluxPerturbatedEVandRdagFirstOrderInit(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress){

    updateFluxVec(x);

    pullEigensystemfromMathematica(0,eigenValues,Rdag,eVDegenAddress);

    std::ofstream eigenvalueCheck(filenameEv.c_str(),std::ofstream::app);
	eigenvalueCheck << std::setprecision(16) << x << "\t";
	for(int i=0;i<subHSDimension;i++) eigenvalueCheck << eigenValues(i,0) << "\t";
	eigenvalueCheck << std::endl;
	eigenvalueCheck.close();

	for(int i=0;i<numbFluxParameters;i++){

		Eigen::VectorXd fluxPerturbations = Eigen::VectorXd::Zero(numbFluxParameters);

		int k = 0;
		for(int j=0;j<numbDerivativePoints-1;j++){

			if(j==derivGridMidpoint) k=1;
			fluxPerturbations(i) = (-1*derivGridMidpoint + j + k)*OPTIMAL_DERIVATIVE_H;

			fluxVec += fluxPerturbations;

            pullEigensystemfromMathematica((numbDerivativePoints-1)*i+1+j,eigenValues,Rdag,eVDegenAddress);

            fluxVec -= fluxPerturbations;
		}

	}

#ifdef DEBUG_ENABLED

	for(int i=0;i<numbHMatrices;i++){

        std::stringstream ss;
        ss << i;
        std::string s;
        ss >> s;
        std::string filenamePert = "eV" + s + "Pert.dat";
        std::ofstream eigenvalueCheckt(filenamePert.c_str(),std::ofstream::app);
        eigenvalueCheckt << std::setprecision(16) << x << "\t";
        for(int k=0;k<subHSDimension;k++) eigenvalueCheckt << eigenValues(k,i) << "\t";
        eigenvalueCheckt << std::endl;

        eigenvalueCheckt.close();

	}

#endif

	return;

}




void RTCGen::setFluxPerturbatedEVandRdagHigherOrderInit(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,std::string& filenameEv,Eigen::MatrixXi& eVDegenAddress){

	Eigen::VectorXi gridPosition(numbFluxParameters);

	for(int i=0;i<numbFluxParameters;i++) gridPosition(i) = -1*derivGridMidpoint;

	gridPosition(0)--;

	Eigen::VectorXd fluxPerturbations(numbFluxParameters);

	updateFluxVec(x);

	for(int i=0;i<numbHMatrices;i++){

		iterate(gridPosition);

		for(int j=0;j<numbFluxParameters;j++) fluxPerturbations(j) = OPTIMAL_DERIVATIVE_H * gridPosition(j);

		fluxVec += fluxPerturbations;

        pullEigensystemfromMathematica(HIConvert(gridPosition),eigenValues,Rdag,eVDegenAddress);

        if(i==(numbHMatrices-1)/2){

            std::ofstream eigenvalueCheck(filenameEv.c_str(),std::ofstream::app);
            eigenvalueCheck << std::setprecision(16) << x << "\t";
            for(int ii=0;ii<subHSDimension;ii++) eigenvalueCheck << eigenValues(ii,i) << "\t";
            eigenvalueCheck << std::endl;
            eigenvalueCheck.close();

		}

        fluxVec -= fluxPerturbations;

	}

#ifdef DEBUG_ENABLED

	for(int i=0;i<numbHMatrices;i++){

        std::stringstream ss;
        ss << i;
        std::string s;
        ss >> s;
        std::string filenamePert = "eV" + s + "Pert.dat";
        std::ofstream eigenvalueCheckt(filenamePert.c_str(),std::ofstream::app);
        eigenvalueCheckt << std::setprecision(16) << x << "\t";
        for(int k=0;k<subHSDimension;k++) eigenvalueCheckt << eigenValues(k,i) << "\t";
        eigenvalueCheckt << std::endl;

        eigenvalueCheckt.close();

	}

#endif

	return;

}




void RTCGen::popVecEnt(Eigen::VectorXi& v,int i){

    Eigen::VectorXi tempV(v.size()-1);
    int k = 0;
    for(int j=0;j<v.size();j++){

        if(j!=i) tempV(k) = v(j);
        else continue;
        k++;

    }

    v.resize(v.size()-1);
    v = tempV;

    return;

}




void RTCGen::iterate(Eigen::VectorXi& gridPosition){

	gridPosition(0)++;

	for(int i=0;i<gridPosition.size();i++){

		if(gridPosition(i)>derivGridMidpoint){

			gridPosition(i) = -derivGridMidpoint;
			gridPosition(i+1)++;

		}

	}

	return;
}




void RTCGen::iterateDegenEv(Eigen::VectorXi& gridPosition,bool dir[]){


	if(dir[0] == true) gridPosition(0)++;
	else gridPosition(0)--;

	for(int i=0;i<gridPosition.size();i++){

		if(std::abs(gridPosition(i)) > derivGridMidpoint){

			if(dir[i] == true){
                gridPosition(i)--;
                dir[i] = false;
			}
			else{
                gridPosition(i)++;
                dir[i] = true;
			}
			if(dir[i+1] == true && i<gridPosition.size()-1) gridPosition(i+1)++;
			else if(i<gridPosition.size()-1) gridPosition(i+1)--;

		}

	}

    return;

}




void RTCGen::checkDegenEvPerturbHigherOrder(Eigen::MatrixXcd& Rdag,Eigen::MatrixXd& eigenValues,std::string& filename,Eigen::MatrixXi& eVDegenAddress){

#ifdef DEBUG_ENABLED
    std::ofstream evOverDerivIntervalCheck(filename.c_str());
#endif
    Eigen::VectorXi gridPosition1 = Eigen::VectorXi::Zero(numbFluxParameters);
    bool dir[numbFluxParameters];

    for(int i=0;i<numbFluxParameters;i++){

        gridPosition1(i) = -1*(numbDerivativePoints-1)/2;
        dir[i] = true;

    }



    for(int i=0;i<std::pow(numbDerivativePoints,numbFluxParameters)-1;i++){

        Eigen::VectorXi gridPosition2 = gridPosition1;
        bool dir2[numbFluxParameters];
        for(int j=0;j<numbFluxParameters;j++) dir2[j] = dir[j];
        iterateDegenEv(gridPosition2,dir2);

        int position1 = HIConvert(gridPosition1);
        int position2 = HIConvert(gridPosition2);
#ifdef DEBUG_ENABLED
        evOverDerivIntervalCheck << "GP\n" << gridPosition1 << "\n\n" << gridPosition2 << "\n";
#endif
        iterateDegenEv(gridPosition1,dir);

        Eigen::VectorXi swappedEigenvalues;
        for(int k=0;k<subHSDimension;k++){
            double Norm = normCompare(Rdag.col(position1*subHSDimension+k),Rdag.col(position2*subHSDimension+k));
#ifdef DEBUG_ENABLED
            evOverDerivIntervalCheck << Norm << "\t";
#endif
            if(Norm < PTOL) {
                swappedEigenvalues.conservativeResize(swappedEigenvalues.size()+1);
                swappedEigenvalues(swappedEigenvalues.size()-1) = k;

            }
        }
#ifdef DEBUG_ENABLED
        evOverDerivIntervalCheck << std::endl << swappedEigenvalues << std::endl;
#endif
        while(swappedEigenvalues.size() > 1){
            for(int kk=1;kk<swappedEigenvalues.size();kk++){
                double Norm = normCompare(Rdag.col(position1*subHSDimension+swappedEigenvalues(0)),Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));
                if(Norm > PTOL){
                    int tempInt = eVDegenAddress(swappedEigenvalues(0),position2);
                    eVDegenAddress(swappedEigenvalues(0),position2) = eVDegenAddress(swappedEigenvalues(kk),position2);
                    eVDegenAddress(swappedEigenvalues(kk),position2) = tempInt;

                    Rdag.col(position2*subHSDimension+swappedEigenvalues(0)).swap(Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));
                    double tempDouble = eigenValues(swappedEigenvalues(0),position2);
                    eigenValues(swappedEigenvalues(0),position2) = eigenValues(swappedEigenvalues(kk),position2);
                    eigenValues(swappedEigenvalues(kk),position2) = tempDouble;

                    double NormP = normCompare(Rdag.col(position1*subHSDimension+swappedEigenvalues(kk)),Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));
                    if(NormP > PTOL) popVecEnt(swappedEigenvalues,kk);
                    popVecEnt(swappedEigenvalues,0);

                    break;
                }
            }
        }
#ifdef DEBUG_ENABLED
        for(int k=0;k<subHSDimension;k++){
            double Norm = normCompare(Rdag.col(position1*subHSDimension+k),Rdag.col(position2*subHSDimension+k));
            evOverDerivIntervalCheck << Norm << "\t";
        }
        evOverDerivIntervalCheck <<  "\n"<< eVDegenAddress;
        evOverDerivIntervalCheck << "\n\n";
#endif

    }


#ifdef DEBUG_ENABLED
    evOverDerivIntervalCheck.close();
#endif
    return;

}




void RTCGen::checkDegenEvPerturbFirstOrder(Eigen::MatrixXcd& Rdag,Eigen::MatrixXd& eigenValues,std::string& filename,Eigen::MatrixXi& eVDegenAddress){

#ifdef DEBUG_ENABLED
    std::ofstream evOverDerivIntervalCheck(filename.c_str());
#endif
    int branchLength = (numbDerivativePoints-1)/2;

    for(int i=0;i<numbFluxParameters;i++){

        for(int j=0;j<branchLength;j++){

            int position1,position2;
            position1 = (numbDerivativePoints-1)*i + branchLength + 0 + j;
            position2 = (numbDerivativePoints-1)*i + branchLength + 1 + j;
            if(j==0) position1 = 0;
#ifdef DEBUG_ENABLED
            evOverDerivIntervalCheck << position1 << "\t" << position2 << "\n";
#endif
            Eigen::VectorXi swappedEigenvalues;
            for(int k=0;k<subHSDimension;k++){
                double Norm = normCompare(Rdag.col(position1*subHSDimension+k),Rdag.col(position2*subHSDimension+k));
#ifdef DEBUG_ENABLED
                evOverDerivIntervalCheck << Norm << "\t";
#endif
                if(Norm < PTOL) {
                    swappedEigenvalues.conservativeResize(swappedEigenvalues.size()+1);
                    swappedEigenvalues(swappedEigenvalues.size()-1) = k;

                }
            }
#ifdef DEBUG_ENABLED
            evOverDerivIntervalCheck << std::endl << swappedEigenvalues << std::endl;
#endif
            while(swappedEigenvalues.size() > 1){
                for(int kk=1;kk<swappedEigenvalues.size();kk++){
                    double Norm = normCompare(Rdag.col(position1*subHSDimension+swappedEigenvalues(0)),Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));
                    if(Norm > PTOL){
                        int tempInt = eVDegenAddress(swappedEigenvalues(0),position2);
                        eVDegenAddress(swappedEigenvalues(0),position2) = eVDegenAddress(swappedEigenvalues(kk),position2);
                        eVDegenAddress(swappedEigenvalues(kk),position2) = tempInt;

                        Rdag.col(position2*subHSDimension+swappedEigenvalues(0)).swap(Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));
                        double tempDouble = eigenValues(swappedEigenvalues(0),position2);
                        eigenValues(swappedEigenvalues(0),position2) = eigenValues(swappedEigenvalues(kk),position2);
                        eigenValues(swappedEigenvalues(kk),position2) = tempDouble;


                        double NormP = normCompare(Rdag.col(position1*subHSDimension+swappedEigenvalues(kk)),Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));
                        if(NormP > PTOL) popVecEnt(swappedEigenvalues,kk);
                        popVecEnt(swappedEigenvalues,0);

                        break;
                    }
                }
            }
#ifdef DEBUG_ENABLED
            for(int k=0;k<subHSDimension;k++){
                double Norm = normCompare(Rdag.col(position1*subHSDimension+k),Rdag.col(position2*subHSDimension+k));
                evOverDerivIntervalCheck << Norm << "\t";
            }
            evOverDerivIntervalCheck <<  "\n"<< eVDegenAddress;
            evOverDerivIntervalCheck << "\n\n";
#endif

            position1 = (numbDerivativePoints-1)*i + branchLength + 1 - j;
            position2 = (numbDerivativePoints-1)*i + branchLength + 0 - j;
            if(j==0) position1 = 0;

#ifdef DEBUG_ENABLED
            evOverDerivIntervalCheck << position1 << "\t" << position2 << "\n";
#endif
            swappedEigenvalues.resize(0);
            for(int k=0;k<subHSDimension;k++){
                double Norm = normCompare(Rdag.col(position1*subHSDimension+k),Rdag.col(position2*subHSDimension+k));
#ifdef DEBUG_ENABLED
                evOverDerivIntervalCheck << Norm << "\t";
#endif
                if(Norm < PTOL) {
                    swappedEigenvalues.conservativeResize(swappedEigenvalues.size()+1);
                    swappedEigenvalues(swappedEigenvalues.size()-1) = k;

                }
            }
#ifdef DEBUG_ENABLED
            evOverDerivIntervalCheck << std::endl << swappedEigenvalues << std::endl;
#endif
            while(swappedEigenvalues.size() > 1){
                for(int kk=1;kk<swappedEigenvalues.size();kk++){
                    double Norm = normCompare(Rdag.col(position1*subHSDimension+swappedEigenvalues(0)),Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));
                    if(Norm > PTOL){
                        int tempInt = eVDegenAddress(swappedEigenvalues(0),position2);
                        eVDegenAddress(swappedEigenvalues(0),position2) = eVDegenAddress(swappedEigenvalues(kk),position2);
                        eVDegenAddress(swappedEigenvalues(kk),position2) = tempInt;

                        Rdag.col(position2*subHSDimension+swappedEigenvalues(0)).swap(Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));

                        double tempDouble = eigenValues(swappedEigenvalues(0),position2);
                        eigenValues(swappedEigenvalues(0),position2) = eigenValues(swappedEigenvalues(kk),position2);
                        eigenValues(swappedEigenvalues(kk),position2) = tempDouble;

                        double NormP = normCompare(Rdag.col(position1*subHSDimension+swappedEigenvalues(kk)),Rdag.col(position2*subHSDimension+swappedEigenvalues(kk)));
                        if(NormP > PTOL) popVecEnt(swappedEigenvalues,kk);
                        popVecEnt(swappedEigenvalues,0);

                        break;
                    }
                }
            }
#ifdef DEBUG_ENABLED
            for(int k=0;k<subHSDimension;k++){
                double Norm = normCompare(Rdag.col(position1*subHSDimension+k),Rdag.col(position2*subHSDimension+k));
                evOverDerivIntervalCheck << Norm << "\t";
            }
            evOverDerivIntervalCheck <<  "\n"<< eVDegenAddress;
            evOverDerivIntervalCheck << "\n\n";
#endif
        }

    }
#ifdef DEBUG_ENABLED
    evOverDerivIntervalCheck.close();
#endif
    return;
}




void RTCGen::correctInitialEigenDegeneracy(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Rdag,Eigen::MatrixXi& eVDegenAddress){

    std::string filename = "eVInitialDegenCheck.dat";

    if(taylorExpansionOrder==0) assert(2>1);
    else if(taylorExpansionOrder==1) checkDegenEvPerturbFirstOrder(Rdag,eigenValues,filename,eVDegenAddress);
    else checkDegenEvPerturbHigherOrder(Rdag,eigenValues,filename,eVDegenAddress);
    return;

}




void RTCGen::correctInitialNDSamplingDegeneracy(Eigen::MatrixXcd& RdagPrev,Eigen::MatrixXcd& RdagCurr,Eigen::MatrixXd& eigenValues,Eigen::MatrixXi& eVDegenAddress){

#ifdef DEBUG_ENABLED
    std::ofstream eVDegenCheck("EvDegenNDSamplingCheck.dat",std::ofstream::app);
#endif
    for(int i=0;i<numbHMatrices;i++){

        Eigen::VectorXi swappingEV;

        for(int k=0;k<subHSDimension;k++){

            double Norm = normCompare(RdagPrev.col(i*subHSDimension+k),RdagCurr.col(i*subHSDimension+k));
#ifdef DEBUG_ENABLED
            eVDegenCheck << Norm << "\t";
#endif
            if(Norm < PTOL){

                swappingEV.conservativeResize(swappingEV.size()+1);
                swappingEV(swappingEV.size()-1) = k;

            }

        }
#ifdef DEBUG_ENABLED
        eVDegenCheck << std::endl << swappingEV << std::endl;
#endif
        while(swappingEV.size() > 1){
            for(int kk=1;kk<swappingEV.size();kk++){

                double Norm = normCompare(RdagPrev.col(i*subHSDimension+swappingEV(0)),RdagCurr.col(i*subHSDimension+swappingEV(kk)));

                if(Norm > PTOL){

                    int tempInt = eVDegenAddress(swappingEV(0),i);
                    eVDegenAddress(swappingEV(0),i) = eVDegenAddress(swappingEV(kk),i);
                    eVDegenAddress(swappingEV(kk),i) = tempInt;

                    RdagCurr.col(i*subHSDimension+swappingEV(0)).swap(RdagCurr.col(i*subHSDimension+swappingEV(kk)));

                    double tempDouble = eigenValues(swappingEV(0),i);
                    eigenValues(swappingEV(0),i) = eigenValues(swappingEV(kk),i);
                    eigenValues(swappingEV(kk),i) = tempDouble;

                    double NormP = normCompare(RdagPrev.col(i*subHSDimension+swappingEV(kk)),RdagCurr.col(i*subHSDimension+swappingEV(kk)));
                    if(NormP > PTOL) popVecEnt(swappingEV,kk);
                    popVecEnt(swappingEV,0);

                    break;
                }
            }
        }

#ifdef DEBUG_ENABLED
        for(int k=0;k<subHSDimension;k++){
            double Norm = normCompare(RdagPrev.col(i*subHSDimension+k),RdagCurr.col(i*subHSDimension+k));
            eVDegenCheck << Norm << "\t";
        }

        eVDegenCheck <<  "\n"<< eVDegenAddress;
        eVDegenCheck << "\n\n";

        eVDegenCheck << "\n";
#endif
    }
#ifdef DEBUG_ENABLED
    eVDegenCheck.close();
#endif
    return;

}




void RTCGen::initializeNDSamplingInterval(Eigen::MatrixXd eigenValues[],Eigen::MatrixXcd Rdag[],Eigen::MatrixXcd& RdagPrime,Eigen::MatrixXd& eigenValuesProj,Eigen::MatrixXi& eVDegenAddress,std::string filename){

	eigenValues[0].resize(subHSDimension,numbHMatrices);

	Rdag[0].resize(totalHSDimension,subHSDimension*numbHMatrices);

	if(taylorExpansionOrder == 0) setFluxPerturbatedEVandRdagZerothOrderInit(eigenValues[0],Rdag[0],filename,eVDegenAddress);
	else if(taylorExpansionOrder == 1) setFluxPerturbatedEVandRdagFirstOrderInit(eigenValues[0],Rdag[0],filename,eVDegenAddress);
	else setFluxPerturbatedEVandRdagHigherOrderInit(eigenValues[0],Rdag[0],filename,eVDegenAddress);

    correctInitialEigenDegeneracy(eigenValues[0],Rdag[0],eVDegenAddress);

	correctInitialEigenbasisPhase(Rdag[0]);

	x += dx;


	for(int i=1;i<numbRdagDerivativePoints;i++){

		eigenValues[i].resize(subHSDimension,numbHMatrices);

		Rdag[i].resize(totalHSDimension,subHSDimension*numbHMatrices);

		if(taylorExpansionOrder == 0) setFluxPerturbatedEVandRdagZerothOrderInit(eigenValues[i],Rdag[i],filename,eVDegenAddress);
		else if(taylorExpansionOrder == 1) setFluxPerturbatedEVandRdagFirstOrderInit(eigenValues[i],Rdag[i],filename,eVDegenAddress);
		else setFluxPerturbatedEVandRdagHigherOrderInit(eigenValues[i],Rdag[i],filename,eVDegenAddress);

		correctInitialNDSamplingDegeneracy(Rdag[i-1],Rdag[i],eigenValues[i],eVDegenAddress);

		for(int j=0;j<subHSDimension*numbHMatrices;j++){

			double phaseArg = std::arg(Rdag[i](0,j));
            Rdag[i].col(j) *= exp(-I*phaseArg);

		}

		x += dx;

	}

	return;
}




void RTCGen::removeOldFiles(std::string& filename, std::string filenameEv){

	remove(filename.c_str());
	remove(filenameEv.c_str());
    remove("EvDegenNDSamplingCheck.dat");
    remove("eVInitialDegenCheck.dat");
    remove("Taylor_Coefficient_Order.dat");
    remove("finalEb.dat");

    int numbFilestoDelete = 0;
    while(true){

        std::stringstream ss;
        ss << numbFilestoDelete;
        std::string s;
        ss >> s;
        std::string filenamePert = "eV" + s + "Pert.dat";
        std::ifstream fileCheck(filenamePert.c_str());
        if(!fileCheck.is_open()) break;
        fileCheck.close();
        numbFilestoDelete++;

    }

	for(int i=0;i<numbFilestoDelete;i++){

        std::stringstream ss;
        ss << i;
        std::string s;
        ss >> s;
        std::string filenamePert = "eV" + s + "Pert.dat";
        remove(filenamePert.c_str());
        std::string filenameProj = "eV"  + s + "Proj.dat";
        remove(filenameProj.c_str());


	}

    return;

}




void RTCGen::printFinalBasis(Eigen::MatrixXcd& Rdag){

    std::ofstream finalBasisOut("finalEb.dat");
    Eigen::VectorXi zeroVecInt = Eigen::VectorXi::Zero(numbFluxParameters);
    if(taylorExpansionOrder==0)  finalBasisOut << Rdag.block(0,0,totalHSDimension,subHSDimension) << std::endl;
    else if(taylorExpansionOrder==1) finalBasisOut << Rdag.block(0,0,totalHSDimension,subHSDimension) << std::endl;
    else finalBasisOut << Rdag.block(0,HIConvert(zeroVecInt)*subHSDimension,totalHSDimension,subHSDimension) << std::endl;
    finalBasisOut.close();

    return;
}




inline void RTCGen::addEigenvaluestoDiagonals(Eigen::MatrixXd& eigenValues,Eigen::MatrixXcd& Htot){

	for(int i=0;i<subHSDimension;i++){
		for(int j=0;j<numbHMatrices;j++){
			Htot(i,i+subHSDimension*j) += eigenValues(i,j);
		}
	}

	return;

}




inline void RTCGen::buildHtotFirstOrder(Eigen::MatrixXcd& Htot,Eigen::MatrixXcd& HtotTemp,finiteDifferenceFormula& numericalDerivative){

	for(int j=0;j<numbFluxParameters;j++){

		Htot.block(0,subHSDimension*(j+1),subHSDimension,subHSDimension) = numericalDerivative.weights(derivGridMidpoint) * HtotTemp.block(0,0,subHSDimension,subHSDimension);

		int k=0;
		for(int i=0;i<numbDerivativePoints-1;i++){
			if(i==derivGridMidpoint) k++;
			Htot.block(0,subHSDimension*(j+1),subHSDimension,subHSDimension) += numericalDerivative.weights(i+k) * HtotTemp.block(0,subHSDimension*(1+i+j*(numbDerivativePoints-1)),subHSDimension,subHSDimension);
		}
	}

	return;

}




int RTCGen::HIConvert(Eigen::VectorXi& gridPosition){

	int output = 0;

	for(int i=0;i<gridPosition.size();i++){
		output += (gridPosition(i) + derivGridMidpoint) * std::pow(numbDerivativePoints,i);
	}

	return output;
}




inline void RTCGen::collapseDim(int& dim,int& derivOrder,Eigen::VectorXi& gridPosition,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXcd& reducedHyperCube,int& loc,Eigen::MatrixXcd& HtotTemp){

	gridPosition(dim) = -derivGridMidpoint;

	reducedHyperCube.block(0,loc*subHSDimension,subHSDimension,subHSDimension) = numericalDerivative[derivOrder-1].weights(0) * HtotTemp.block(0,HIConvert(gridPosition)*subHSDimension,subHSDimension,subHSDimension);

	for(int i=1;i<numbDerivativePoints;i++){
		gridPosition(dim)++;
		reducedHyperCube.block(0,loc*subHSDimension,subHSDimension,subHSDimension) += numericalDerivative[derivOrder-1].weights(i) * HtotTemp.block(0,HIConvert(gridPosition)*subHSDimension,subHSDimension,subHSDimension);
	}


	return;
}




inline void RTCGen::setReducedHyperCubeInitial(Eigen::VectorXi& v,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXcd& reducedHyperCube,Eigen::MatrixXcd& HtotTemp){

	Eigen::VectorXi nonZeros;
	for(int i=0;i<v.size();i++){
		if(v(i) != 0){
			nonZeros.conservativeResize(nonZeros.size()+1);
			nonZeros(nonZeros.size()-1) = i;
		}
	}

	Eigen::VectorXi gridPosition = Eigen::VectorXi::Zero(v.size());
	Eigen::VectorXi subGridPosition(nonZeros.size()-1);

	for(int i=1;i<nonZeros.size();i++) subGridPosition(i-1) = -derivGridMidpoint;
	for(int j=1;j<nonZeros.size();j++) gridPosition(nonZeros(j)) = subGridPosition(j-1);

	int upperLim = std::pow(numbDerivativePoints,nonZeros.size()-1);
	for(int i=0;i<upperLim;i++){

		collapseDim(nonZeros(0),v(nonZeros(0)),gridPosition,numericalDerivative,reducedHyperCube,i,HtotTemp);
		if(subGridPosition.size()>0 && i<upperLim-1) iterate(subGridPosition);
		for(int j=1;j<nonZeros.size();j++) gridPosition(nonZeros(j)) = subGridPosition(j-1);

	}

	for(int i=1;i<nonZeros.size();i++) v(i-1) = v(nonZeros(i));
	v.conservativeResize(nonZeros.size()-1);

	return;
}




inline void RTCGen::setReducedHyperCube(Eigen::VectorXi& v,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXcd& reducedHyperCubeIn,Eigen::MatrixXcd& reducedHyperCubeOut){

	int zero = 0;
	Eigen::VectorXi gridPosition(v.size());
	Eigen::VectorXi subGridPosition(v.size()-1);

	for(int i=1;i<gridPosition.size();i++) subGridPosition(i-1) = -derivGridMidpoint;
	for(int j=1;j<gridPosition.size();j++) gridPosition(j) = subGridPosition(j-1);

	reducedHyperCubeOut.resize(subHSDimension,subHSDimension*std::pow(numbDerivativePoints,gridPosition.size()-1));

	int upperLim = std::pow(numbDerivativePoints,gridPosition.size()-1);
	for(int i=0;i<upperLim;i++){

		collapseDim(zero,v(0),gridPosition,numericalDerivative,reducedHyperCubeOut,i,reducedHyperCubeIn);
		if(subGridPosition.size()>0 && i<upperLim-1) iterate(subGridPosition);
		for(int j=1;j<gridPosition.size();j++) gridPosition(j) = subGridPosition(j-1);

	}

	for(int i=1;i<gridPosition.size();i++) v(i-1) = v(i);
	v.conservativeResize(gridPosition.size()-1);
	return;
}




inline void RTCGen::takePartialDerivative(Eigen::MatrixXcd& Htot,Eigen::MatrixXcd& HtotTemp,Eigen::VectorXi v,int& k,finiteDifferenceFormula numericalDerivative[]){

	int newDimension = 1;

	for(int i=0;i<numbFluxParameters;i++){
		newDimension *= std::pow(numbDerivativePoints,std::min(v(i),1));
	}

	newDimension /= numbDerivativePoints;

	Eigen::MatrixXcd reducedHyperCube1(subHSDimension,subHSDimension*newDimension);
	Eigen::MatrixXcd reducedHyperCube2;

	setReducedHyperCubeInitial(v,numericalDerivative,reducedHyperCube1,HtotTemp);

	bool takeRHC1toRHC2 = false;

	while(v.size()>0){

		takeRHC1toRHC2 = !takeRHC1toRHC2;
		if(takeRHC1toRHC2) setReducedHyperCube(v,numericalDerivative,reducedHyperCube1,reducedHyperCube2);
		else setReducedHyperCube(v,numericalDerivative,reducedHyperCube2,reducedHyperCube1);

	}


	if(takeRHC1toRHC2) Htot.block(0,subHSDimension*k,subHSDimension,subHSDimension) = reducedHyperCube2;
	else Htot.block(0,subHSDimension*k,subHSDimension,subHSDimension) = reducedHyperCube1;
	return;
}




inline void RTCGen::buildHtotHigherOrder(Eigen::MatrixXcd& Htot,Eigen::MatrixXcd& HtotTemp,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXi partialDerivativeIndexMatrix[]){

	int k=1;
	for(int i=0;i<taylorExpansionOrder;i++){
		for(int j=0;j<partialDerivativeIndexMatrix[i].rows();j++){


			takePartialDerivative(Htot,HtotTemp,partialDerivativeIndexMatrix[i].row(j),k,numericalDerivative);

			k++;
		}
	}

	return;
}




inline void RTCGen::setHMatrices(Eigen::MatrixXcd& Htot,Eigen::MatrixXcd& Rdag,Eigen::MatrixXcd& RdagPrime,Eigen::MatrixXd& eigenValues,finiteDifferenceFormula numericalDerivative[],Eigen::MatrixXi partialDerivativeIndexMatrix[]){


	Eigen::MatrixXcd HtotTemp(subHSDimension,subHSDimension*numbHMatrices);

	for(int i=0;i<numbHMatrices;i++){

		HtotTemp.block(0,i*subHSDimension,subHSDimension,subHSDimension) = -I/(2*M_PI) * (Rdag.block(0,i*subHSDimension,totalHSDimension,subHSDimension).conjugate().transpose() * RdagPrime.block(0,i*subHSDimension,totalHSDimension,subHSDimension));

	}


	addEigenvaluestoDiagonals(eigenValues,HtotTemp);


	if(taylorExpansionOrder==0) Htot.block(0,0,subHSDimension,subHSDimension) = HtotTemp.block(0,0,subHSDimension,subHSDimension);
	else if(taylorExpansionOrder==1) Htot.block(0,0,subHSDimension,subHSDimension) = HtotTemp.block(0,0,subHSDimension,subHSDimension);
	else {

		Eigen::VectorXi zeroVecInt = Eigen::VectorXi::Zero(numbFluxParameters);

		Htot.block(0,0,subHSDimension,subHSDimension) = HtotTemp.block(0,subHSDimension * HIConvert(zeroVecInt),subHSDimension,subHSDimension);

	}

    if(taylorExpansionOrder==0) assert(2>1);
	else if(taylorExpansionOrder==1) buildHtotFirstOrder(Htot,HtotTemp,numericalDerivative[0]);
	else buildHtotHigherOrder(Htot,HtotTemp,numericalDerivative,partialDerivativeIndexMatrix);


	return;

}




void RTCGen::printMatrixToMathematica(Eigen::MatrixXcd M){

    int Mcols,Mrows;
    Mcols = M.cols();
    Mrows = M.rows();
    WSPutFunction(stdlink,"List",Mrows);
    for(int i=0;i<Mrows;i++){
        WSPutFunction(stdlink,"List",Mcols);
        for(int j=0;j<Mcols;j++){
            WSPutFunction(stdlink,"Complex",2);
            WSPutReal64(stdlink,std::real(M(i,j)));
            WSPutReal64(stdlink,std::imag(M(i,j)));
        }
    }

    return;
}




Eigen::MatrixXi RTCGen::generateSubBasisVector(int subPhotons, int subModes){
    int markers = subPhotons + subModes - 1;
    int myints[markers];
    int i = 0;
    while(i<subPhotons){
        myints[i]=1;
        i++;
    }
    while(i<markers){
        myints[i]=0;
        i++;
    }
    Eigen::MatrixXi nv = Eigen::MatrixXi::Zero(g(subPhotons,subModes),subModes);
    i = 0;
    int j,k = 0;
    do {
        j = 0;
        k = 0;
        while(k<markers){
        if(myints[k]==1){
            nv(i,j)=nv(i,j)+1;
        }
        else if(myints[k]==0){
            j++;
        }

        k++;
        }
        i++;
    } while ( std::prev_permutation(myints,myints+markers) );
    return nv;
}




inline int RTCGen::g(int n,int m){
    if(n==0 && m==0){
        return 0;
    }
    else if(n==0 && m>0){
        return 1;
    }

    else{
        return doublefactorial(n+m-1)/(doublefactorial(n)*doublefactorial(m-1));
    }
}




double RTCGen::doublefactorial(int x){

    double total=1.0;
    if (x>=0){
        for(int i=x;i>0;i--){
            total=i*total;
        }
    }
    else{
        std::cout << "invalid factorial" << std::endl;
        total=-1;
    }
    return total;

}

