//
// Created by Lakshit Dabas on 12/16/19.
//

#include "gates.h"
#include <complex>
using namespace std;
#include "Eigen/Dense"
#include <iostream>


using namespace Eigen;

MatrixXcd makeHadamard() {
    complex<double> c(1, 0);
    MatrixXd _2x2(2,2);
    _2x2(0, 0) = 1/sqrt(2);
    _2x2(0, 1) = 1/sqrt(2);
    _2x2(1, 0) = 1/sqrt(2);
    _2x2(1, 1) = -1/sqrt(2);
    //MatrixXcd C1 = c * _2x2; // multiply complex*real
    MatrixXcd hadamard = c * _2x2; // complex scalar times complex matrix
    return hadamard;
};

MatrixXcd _H = makeHadamard();

Gate H(){ return _H;}