/*! \file finiteDifferenceFuncs.cpp
 *  \brief functions used to generate the coefficients in finite difference formulas.
 *
 *
 *     This is an implementation of Fornberg's Algorithm. See:
 *      Fornberg, Bengt (1988), "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids",
 *      Mathematics of Computation, 51 (184): 699–706, doi:10.1090/S0025-5718-1988-0935077-0, ISSN 0025-5718
 *
 * \author Jake Smith <jsmith74@tulane.edu>
 */



#include "reducedTaylorCoefficientGen.h"
#include <iostream>
#include <cmath>



void finiteDifferenceFormula::setfiniteDifferenceFormula(int NderivPoints,int derivOrder,double pointSpacing){

	weights.resize(NderivPoints);
	M = derivOrder;
	N = NderivPoints - 1;
	alpha.resize(NderivPoints);
	int midpoint = (NderivPoints - 1)/2;
	for(int i=0;i<NderivPoints;i++){

		alpha(i) = -1*midpoint*pointSpacing + i*pointSpacing;

	}

	Eigen::VectorXd delta((N+1)*(N+1)*(M+1));

	delta(0) = 1.0;
	double c1 = 1.0;
	for(int n=1;n<=N;n++){

		double c2 = 1.0;
		for(int v=0;v<=n-1;v++){
			double c3 = alpha(n) - alpha(v);
			c2 *= c3;
			if(n<=M) delta(dIn(n-1,v,n)) = 0.0;
			for(int m=0;m<=std::min(n,M);m++){
				if(m>0) delta(dIn(n,v,m)) = (alpha(n) * delta(dIn(n-1,v,m)) - m * delta(dIn(n-1,v,m-1)))/c3;
				else delta(dIn(n,v,m)) = (alpha(n) * delta(dIn(n-1,v,m)))/c3;
			}
		}
		for(int m=0;m<=std::min(n,M);m++){
			if(m>0) delta(dIn(n,n,m)) = (c1/c2) * (m * delta(dIn(n-1,n-1,m-1)) - alpha(n-1) * delta(dIn(n-1,n-1,m)));
			else delta(dIn(n,n,m)) = -1.0* (c1/c2) * ( alpha(n-1) * delta(dIn(n-1,n-1,m)));
		}

		c1 = c2;
	}

	for(int i=0;i<NderivPoints;i++){

		weights(i) = delta(dIn(N,i,M));

	}

	return;
}

finiteDifferenceFormula::finiteDifferenceFormula(){

}


int finiteDifferenceFormula::dIn(int n,int v,int m){
	return n + v*(N+1)+m*(N+1)*(N+1);
}
