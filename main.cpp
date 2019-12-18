#include <iostream>
#include "Eigen/Dense"
#include <complex>
using namespace std;

using namespace Eigen;


int main() {
    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);

    complex<double> c(1, 0);
    MatrixXd hadamard(2,2);
    hadamard(0, 0) = 1/sqrt(2);
    hadamard(0, 1) = 1/sqrt(2);
    hadamard(1, 0) = 1/sqrt(2);
    hadamard(1, 1) = -1/sqrt(2);
    //MatrixXcd C1 = c * hadamard; // multiply complex*real
    MatrixXcd C2 = c * hadamard; // complex scalar times complex matrix


    std::cout << C2 << std::endl;
    return 0;
}
