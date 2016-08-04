/*! \file finiteDifferenceFuncs.h
 *  \brief header file for finiteDifferenceFuncs.cpp
 *
 *     This file contains function and class definitions for finiteDifferenceFuncs.cpp
 *
 *     This is an implementation of Fornberg's Algorithm. See:
 *      Fornberg, Bengt (1988), "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids",
 *      Mathematics of Computation, 51 (184): 699–706, doi:10.1090/S0025-5718-1988-0935077-0, ISSN 0025-5718
 *
 * \author Jake Smith <jsmith74@tulane.edu>
 */



#include "wstp.h"
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <Eigen/Eigenvalues>
#include <complex>
#include <fstream>


class finiteDifferenceFormula{

    public:

        void setfiniteDifferenceFormula(int NderivPoints,int derivOrder,double pointSpacing);
        finiteDifferenceFormula();
        Eigen::VectorXd weights;

    private:

        int M,N;
        int dIn(int n,int v,int m);
        Eigen::VectorXd alpha;

};
