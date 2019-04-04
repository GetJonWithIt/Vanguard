#ifndef TENSORALGEBRA_H
#define TENSORALGEBRA_H

#include "matrixalgebra.h"
using namespace std;

class TensorAlgebra
{
public:
    TensorAlgebra();

    static double computeFirstInvariant(vector<vector<double> > tensor1);
    static double computeSecondInvariant(vector<vector<double> > tensor1);
    static double computeThirdInvariant(vector<vector<double> > tensor1);

    static vector<vector<double> > computeGramianMatrix(vector<vector<double> > tensor1);
};

#endif // TENSORALGEBRA_H
