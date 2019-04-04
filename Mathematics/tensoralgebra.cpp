#include "tensoralgebra.h"

TensorAlgebra::TensorAlgebra()
{
}

double TensorAlgebra::computeFirstInvariant(vector<vector<double> > tensor1)
{
    return MatrixAlgebra::computeTrace(tensor1);
}

double TensorAlgebra::computeSecondInvariant(vector<vector<double> > tensor1)
{
    return 0.5 * ((MatrixAlgebra::computeTrace(tensor1) * MatrixAlgebra::computeTrace(tensor1)) - MatrixAlgebra::computeComponentSum(MatrixAlgebra::multiplyMatricesElementWise(tensor1, tensor1)));
}

double TensorAlgebra::computeThirdInvariant(vector<vector<double> > tensor1)
{
    return MatrixAlgebra::computeDeterminant(tensor1);
}

vector<vector<double> > TensorAlgebra::computeGramianMatrix(vector<vector<double> > tensor1)
{
    return MatrixAlgebra::multiplyMatrices(MatrixAlgebra::computeTranspose(tensor1), tensor1);
}
