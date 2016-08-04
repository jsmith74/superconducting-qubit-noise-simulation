#include "timeEvolutionGen.h"


int main(){


	const int samplingPoints = 300;
	const double tInit = 0.0;
	const double tFin = 0.2;
	const int numberoffluxparams = 2;
	const int taylorExpansionOrder = 4;
	const int truncHSD = 10;

	Eigen::MatrixXcd H[samplingPoints];

	importH taylorH(samplingPoints,numberoffluxparams,taylorExpansionOrder,truncHSD);

	taylorH.set(H);

    timeEvolutionOperatorGen U(tInit,tFin,samplingPoints,truncHSD,numberoffluxparams);

    U.generate(H,taylorH);

    std::cout << "Resulting Time evolution operator: " << std::endl << std::endl;

    U.printTEOperatorRect();

    U.printTEOperatorPolar();

	return 0;
}

