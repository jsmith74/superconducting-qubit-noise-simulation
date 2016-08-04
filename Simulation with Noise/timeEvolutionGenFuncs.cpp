#include "timeEvolutionGen.h"




 /** FOR NOW MAKE A FUNCTION OF NO FLUX NOISE =========================================== */

void timeEvolutionOperatorGen::updateFluxNoise(double& t){

    fluxNoise(0) = 0.0;
    fluxNoise(1) = 0.0;

    return;

}

/** ===================================================================================== */

timeEvolutionOperatorGen::timeEvolutionOperatorGen(double xStart,double xEnd,int numSamplingPoints,int truncHSD,int numberofFluxP){

    subHSDimension = truncHSD;
    numbPoints = numSamplingPoints;
    tInit = xStart;
    tFin = xEnd;
    dt = (tFin-tInit)/(numbPoints-1);
    std::complex<double> Igen(0.0,1.0);
    I = Igen;
    numbFluxParameters = numberofFluxP;
    numbHMatrices = 4*numbFluxParameters + 1;
    numbPartialDerivMatrices = numbFluxParameters + 1;
    hb = 1.0/(2.0*M_PI);
    fluxNoise.resize(numbFluxParameters);
    Uop.resize(subHSDimension,subHSDimension);

    t = tInit;

    updateFluxNoise(t);

}




void timeEvolutionOperatorGen::generate(Eigen::MatrixXcd Htot[],importH& taylorH){

    Eigen::MatrixXcd H(subHSDimension,subHSDimension);

    H = taylorH.at(Htot[0],fluxNoise);

    Uop = matrixExp(H*dt/(2.0*hb));

    t += dt;

    for(int i=1;i<numbPoints-1;i++){

        updateFluxNoise(t);

        H = taylorH.at(Htot[i],fluxNoise);

        Uop = matrixExp(H*dt/hb) * Uop;

        t += dt;

    }

    updateFluxNoise(t);

    H = taylorH.at(Htot[numbPoints-1],fluxNoise);

    Uop = matrixExp(H*dt/(2.0*hb)) * Uop;

	return;

}

void timeEvolutionOperatorGen::printTEOperatorRect(){

    std::cout << "In rectangular coordinates\n" << Uop << std::endl << std::endl;

    std::cout << "Verify Unitarity: \n" << Uop.conjugate().transpose() * Uop << std::endl << std::endl;
    return;
}

void timeEvolutionOperatorGen::printTEOperatorPolar(){

    std::cout << "In polar coordinates: \n";
	for(int i=0;i<Uop.rows();i++){
		for(int j=0;j<Uop.cols();j++){
			std::cout << "(" << abs(Uop(i,j)) << "," << arg(Uop(i,j)) << ")" << "\t";
		}
		std::cout << std::endl;
	}
    std::cout << std::endl;
    return;
}

Eigen::MatrixXcd timeEvolutionOperatorGen::matrixExp(Eigen::MatrixXcd X){

    int matrixSize;
    matrixSize = X.rows();


    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(X);
    Eigen::VectorXd evalues=ces.eigenvalues();
    Eigen::MatrixXcd evectors=(ces.eigenvectors());
    Eigen::MatrixXcd cevectors=evectors.conjugate();
    Eigen::MatrixXcd sylvester[matrixSize];

    for(int i=0;i < matrixSize;i++){
        sylvester[i].resize(matrixSize,matrixSize);
        for(int m=0; m<matrixSize;m++){
            for(int n=0;n<matrixSize;n++){
                sylvester[i](n,m)=evectors(n,i)*cevectors(m,i);
            }
        }
    }

    Eigen::MatrixXcd result(matrixSize,matrixSize);
    result = exp(-I*evalues(0))*sylvester[0];
    for(int j=1;j<matrixSize;j++){
        result=result+exp(-I*evalues(j))*sylvester[j];
    }

    return result;
}
