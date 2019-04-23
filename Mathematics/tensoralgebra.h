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
    static vector<vector<double> > computeDeviator(vector<vector<double> > tensor1);

    static vector<vector<vector<double> > > reshapeMatrixAtLevelOne(vector<vector<double> > matrix1);
    static vector<vector<vector<double> > > reshapeMatrixAtLevelTwo(vector<vector<double> > matrix1);

    static vector<vector<vector<vector<double> > > > multiplyThreeTensors(vector<vector<vector<double> > > threeTensor1, vector<vector<vector<double> > > threeTensor2);

    static vector<vector<vector<vector<double> > > > swapFourTensorAxes(vector<vector<vector<vector<double> > > > fourTensor1, int axis1, int axis2);

    static vector<vector<vector<vector<double> > > > addFourTensors(vector<vector<vector<vector<double> > > > fourTensor1, vector<vector<vector<vector<double> > > > fourTensor2);
    static vector<vector<vector<vector<double> > > > subtractFourTensors(vector<vector<vector<vector<double> > > > fourTensor1, vector<vector<vector<vector<double> > > > fourTensor2);

    static vector<vector<vector<vector<double> > > > multiplyFourTensor(double scalar, vector<vector<vector<vector<double> > > > fourTensor1);

    static double computeSigmaNorm(vector<vector<double> > tensor1);
};

#endif // TENSORALGEBRA_H
