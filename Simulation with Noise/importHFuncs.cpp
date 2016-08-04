#include "timeEvolutionGen.h"



importH::importH(int samplingPoints,int numberoffluxparams,int TEorder,int truncHSD){

	numbPoints = samplingPoints;
	numbFluxParameters = numberoffluxparams;
	taylorExpansionOrder = TEorder;
	subHSDimension = truncHSD;
	numbPartialDerivMatrices = genNumbPartialDerivMatrices(taylorExpansionOrder,numbFluxParameters);

	partialDerivativeIndexMatrix.resize(numbPartialDerivMatrices-1,numbFluxParameters);

	Eigen::MatrixXi partialDerivativeIndexMatrixImport[taylorExpansionOrder];
	int k = 0;
	for(int i=0;i<taylorExpansionOrder;i++){
		partialDerivativeIndexMatrixImport[i] = generateSubBasisVector(i+1,numbFluxParameters);
		partialDerivativeIndexMatrix.block(k,0,partialDerivativeIndexMatrixImport[i].rows(),numbFluxParameters) = partialDerivativeIndexMatrixImport[i];
		k += partialDerivativeIndexMatrixImport[i].rows();
	}

	std::complex<double> Igen(0.0,1.0);
	I = Igen;
	numbPrevMPIProcesses = genNumPrevMPIProcesses();
	numbPartialDerivMatricesAt.resize(taylorExpansionOrder);
	for(int i=1;i<=taylorExpansionOrder;i++) numbPartialDerivMatricesAt(i-1) = doublefactorial(i+numbFluxParameters-1)/(doublefactorial(i)*doublefactorial(numbFluxParameters-1));

}

int importH::genNumPrevMPIProcesses(){

	int output = 0;

	while(true){

		std::string filename;
		setFilename(filename,output);
		std::ifstream HCount(filename.c_str());
		if(!HCount.is_open()) return output;
		output++;

	}

}

void importH::setFilename(std::string& filename,int fileNumber){

	std::stringstream ss;
	std::string s;
	ss << fileNumber;
	ss >> s;
	filename = "H" + s + ".dat";
	return;
}

void importH::setThisH(Eigen::MatrixXcd& H,std::ifstream& HimportStream){

    int Hcols,Hrows;
    HimportStream >> Hrows;
    HimportStream >> Hcols;
    H.resize(Hrows,Hcols);

	for(int j=0;j<Hcols;j++){
		for(int i=0;i<Hrows;i++){
			HimportStream >> H(i,j);
			double imagPart;
			HimportStream >> imagPart;
			H(i,j) += imagPart*I;
		}
	}


	return;
}

void importH::fixMPIPhaseDiscontinuity(Eigen::MatrixXcd& H,Eigen::VectorXcd& beginPhaseDiscontinuity){


    for(int i=0;i<H.rows();i++){
        H.row(i) *= beginPhaseDiscontinuity(i);
    }

    for(int k=0;k<numbPartialDerivMatrices;k++){

        for(int j=0;j<subHSDimension;j++){

            H.block(0,k*subHSDimension,subHSDimension,subHSDimension).col(j) *= conj(beginPhaseDiscontinuity(j));

        }

    }


    return;

}

void importH::set(Eigen::MatrixXcd H[]){

    Eigen::VectorXcd beginPhaseDiscontinuity(subHSDimension);
    for(int j=0;j<subHSDimension;j++) beginPhaseDiscontinuity(j) = 1.0;


	int k=0;
	for(int i=0;i<numbPrevMPIProcesses;i++){

		std::string filename,filename2;
		setFilename(filename,i);
		setFilename(filename2,i-1);
        filename2 = "endPhase" + filename2;
		std::ifstream HimportStream(filename.c_str());

		for(int j=0;j<numbPoints;j++){

			setThisH(H[k],HimportStream);
			k++;

		}

		HimportStream.close();

		if(i==0) continue;


	}

	return;

}

double importH::taylorArgument(Eigen::VectorXd& testNoiseVector,Eigen::VectorXi partialVector){

	double output = 1.0;
	for(int i=0;i<partialVector.size();i++){
		output *= std::pow(testNoiseVector(i),partialVector(i));
	}
	return output;

}

Eigen::MatrixXcd importH::at(Eigen::MatrixXcd& Htot,Eigen::VectorXd& noiseLocation){

	Eigen::MatrixXcd result(subHSDimension,subHSDimension);

	result = Htot.block(0,0,subHSDimension,subHSDimension);

	int k = 1;
	for(int i=0;i<taylorExpansionOrder;i++){
		for(int j=0;j<numbPartialDerivMatricesAt(i);j++){
			result += (1.0/doublefactorial(i+1)) * Htot.block(0,k*subHSDimension,subHSDimension,subHSDimension) * taylorArgument(noiseLocation,partialDerivativeIndexMatrix.row(k-1));
			k++;
		}
	}

	return result;

}

int importH::genNumbPartialDerivMatrices(int& taylorExpansionOrder,int& numbFluxParameters){

	int output = 0;
	for(int i=0;i<=taylorExpansionOrder;i++){

		output += doublefactorial(i+numbFluxParameters-1)/(doublefactorial(i)*doublefactorial(numbFluxParameters-1));

	}

	return output;

}

inline int importH::g(int n,int m){
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


Eigen::MatrixXi importH::generateSubBasisVector(int subPhotons, int subModes){
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

double importH::doublefactorial(int x){

    double total=1.0;
    assert(x>=0 && "Invalid Factorial");

    for(int i=x;i>0;i--){
        total=i*total;
    }

    return total;

}
