#ifndef MATRIXALGEBRA_H
#define MATRIXALGEBRA_H

#include "vectoralgebra.h"
#include <cmath>
using namespace std;

class MatrixAlgebra
{
public:
    MatrixAlgebra();

    static vector<vector<double> > computeIdentityMatrix(int dimension);

    static vector<vector<double> > addMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
    static vector<vector<double> > subtractMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);

    static vector<vector<double> > multiplyMatrix(double scalar1, vector<vector<double> > matrix1);
    static vector<vector<double> > multiplyMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
    static vector<vector<double> > multiplyMatricesElementWise(vector<vector<double> > matrix1, vector<vector<double> > matrix2);

    static vector<double> multiplyMatrixByVector(vector<vector<double> > matrix1, vector<double> vector1);

    static vector<vector<double> > computeTranspose(vector<vector<double> > matrix1);

    static double computeTrace(vector<vector<double> > matrix1);
    static double computeDeterminant(vector<vector<double> > matrix1);

    static double computeComponentSum(vector<vector<double> > matrix1);

    static vector<vector<double> > computeInverseMatrix(vector<vector<double> > matrix1);

    static double computeLargestEigenvalue(vector<vector<double> > matrix1);
};

#endif // MATRIXALGEBRA_H
