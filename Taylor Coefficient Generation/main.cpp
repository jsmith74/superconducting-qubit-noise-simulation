/*! \file main.cpp
 *  \brief The main function is hijacked by WSTP
 *
 *
 *      These Functions are for loading this program into a Mathematica notebook and executing the Taylor coefficient generation.
 *
 *
 *
 * \author Jake Smith <jsmith74@tulane.edu>
 */



#include "reducedTaylorCoefficientGen.h"


void GenerateTaylorCoefficients(int initialHSD,int truncHSD,int numberoffluxparams,double tInit,double tFin,int taylorExpansionOrder,int samplingPoints,int RdagDerivApprox,int derivApprox,double finiteDifferenceInterval){

	if(derivApprox <= taylorExpansionOrder) WSPutString(stdlink,"Error: The number of grid point used to sample the numerical derivative, dHdphiSamplingPoints, must be strictly larger than the highest order partial derivatives in the Taylor series");

	int processNumber = 0;

	RTCGen HGen(tInit,tFin,samplingPoints,initialHSD,truncHSD,derivApprox,numberoffluxparams,taylorExpansionOrder,RdagDerivApprox,finiteDifferenceInterval,processNumber);

	Eigen::MatrixXcd H;

	HGen.generate(H);

	WSPutString(stdlink,"Successfully exported Reduced Hamiltonian Taylor Coefficients to file.");

    return;

}


int main(int argc, char *argv[]){

	return WSMain(argc, argv);

}
