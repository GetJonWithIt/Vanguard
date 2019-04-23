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

vector<vector<double> > TensorAlgebra::computeDeviator(vector<vector<double> > tensor1)
{
    return MatrixAlgebra::subtractMatrices(tensor1, MatrixAlgebra::multiplyMatrix((MatrixAlgebra::computeTrace(tensor1) / 3.0), MatrixAlgebra::computeIdentityMatrix(3)));
}

vector<vector<vector<double> > > TensorAlgebra::reshapeMatrixAtLevelOne(vector<vector<double> > matrix1)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<vector<double> > > reshapedTensor(rowCount, vector<vector<double> >(columnCount, vector<double>(1)));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            reshapedTensor[i][j][0] = matrix1[i][j];
        }
    }

    return reshapedTensor;
}


vector<vector<vector<double> > > TensorAlgebra::reshapeMatrixAtLevelTwo(vector<vector<double> > matrix1)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<vector<double> > > reshapedTensor(rowCount, vector<vector<double> >(1, vector<double>(columnCount)));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            reshapedTensor[i][0][j] = matrix1[i][j];
        }
    }

    return reshapedTensor;
}

vector<vector<vector<vector<double> > > > TensorAlgebra::multiplyThreeTensors(vector<vector<vector<double> > > threeTensor1, vector<vector<vector<double> > > threeTensor2)
{
    int threeTensor1RowCount = threeTensor1.size();
    int threeTensor1ColumnCount = threeTensor1[0].size();
    int threeTensor1LayerCount = threeTensor1[0][0].size();;

    int threeTensor2RowCount = threeTensor2.size();
    int threeTensor2LayerCount = threeTensor2[0][0].size();

    vector<vector<vector<vector<double> > > > productTensor(threeTensor1RowCount, vector<vector<vector<double> > >(threeTensor1ColumnCount, vector<vector<double> >(threeTensor2RowCount, vector<double>(threeTensor2LayerCount))));

    for (int i = 0; i < threeTensor1RowCount; i++)
    {
        for (int j = 0; j < threeTensor1ColumnCount; j++)
        {
            for (int k = 0; k < threeTensor2RowCount; k++)
            {
                for (int l = 0; l < threeTensor2LayerCount; l++)
                {
                    for (int m = 0; m < threeTensor1LayerCount; m++)
                    {
                        productTensor[i][j][k][l] += threeTensor1[i][j][m] * threeTensor2[k][m][l];
                    }
                }
            }
        }
    }

    return productTensor;
}

vector<vector<vector<vector<double> > > > TensorAlgebra::swapFourTensorAxes(vector<vector<vector<vector<double> > > > fourTensor1, int axis1, int axis2)
{
    int rowCount = fourTensor1.size();
    int columnCount = fourTensor1[0].size();
    int layer1Count = fourTensor1[0][0].size();
    int layer2Count = fourTensor1[0][0][0].size();

    if ((axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 ==0))
    {
        vector<vector<vector<vector<double> > > > fourTensorSwapped(columnCount, vector<vector<vector<double> > >(rowCount, vector<vector<double> >(layer1Count, vector<double>(layer2Count))));

        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                for (int k = 0; k < layer1Count; k++)
                {
                    for (int l = 0; l < layer2Count; l++)
                    {
                        fourTensorSwapped[j][i][k][l] = fourTensor1[i][j][k][l];
                    }
                }
            }
        }

        return fourTensorSwapped;
    }
    else if ((axis1 == 0 && axis2 == 2) || (axis1 == 2 && axis2 == 0))
    {
        vector<vector<vector<vector<double> > > > fourTensorSwapped(layer1Count, vector<vector<vector<double> > >(columnCount, vector<vector<double> >(rowCount, vector<double>(layer2Count))));

        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                for (int k = 0; k < layer1Count; k++)
                {
                    for (int l = 0; l < layer2Count; l++)
                    {
                        fourTensorSwapped[k][j][i][l] = fourTensor1[i][j][k][l];
                    }
                }
            }
        }

        return fourTensorSwapped;
    }
    else if ((axis1 == 0 && axis2 == 3) || (axis1 == 3 && axis2 == 0))
    {
        vector<vector<vector<vector<double> > > > fourTensorSwapped(layer2Count, vector<vector<vector<double> > >(columnCount, vector<vector<double> >(layer1Count, vector<double>(rowCount))));

        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                for (int k = 0; k < layer1Count; k++)
                {
                    for (int l = 0; l < layer2Count; l++)
                    {
                        fourTensorSwapped[l][j][k][i] = fourTensor1[i][j][k][l];
                    }
                }
            }
        }

        return fourTensorSwapped;
    }
    else if ((axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1))
    {
        vector<vector<vector<vector<double> > > > fourTensorSwapped(rowCount, vector<vector<vector<double> > >(layer1Count, vector<vector<double> >(columnCount, vector<double>(layer2Count))));

        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                for (int k = 0; k < layer1Count; k++)
                {
                    for (int l = 0; l < layer2Count; l++)
                    {
                        fourTensorSwapped[i][k][j][l] = fourTensor1[i][j][k][l];
                    }
                }
            }
        }

        return fourTensorSwapped;
    }
    else if ((axis1 == 1 && axis2 == 3) || (axis1 == 3 && axis2 == 1))
    {
        vector<vector<vector<vector<double> > > > fourTensorSwapped(rowCount, vector<vector<vector<double> > >(layer2Count, vector<vector<double> >(layer1Count, vector<double>(columnCount))));

        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                for (int k = 0; k < layer1Count; k++)
                {
                    for (int l = 0; l < layer2Count; l++)
                    {
                        fourTensorSwapped[i][l][k][j] = fourTensor1[i][j][k][l];
                    }
                }
            }
        }

        return fourTensorSwapped;
    }
    else if ((axis1 == 2 && axis2 == 3) || (axis1 == 3 && axis2 == 2))
    {
        vector<vector<vector<vector<double> > > > fourTensorSwapped(rowCount, vector<vector<vector<double> > >(columnCount, vector<vector<double> >(layer2Count, vector<double>(layer1Count))));

        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
            {
                for (int k = 0; k < layer1Count; k++)
                {
                    for (int l = 0; l < layer2Count; l++)
                    {
                        fourTensorSwapped[i][j][l][k] = fourTensor1[i][j][k][l];
                    }
                }
            }
        }

        return fourTensorSwapped;
    }
    else
    {
        return fourTensor1;
    }
}

vector<vector<vector<vector<double> > > > TensorAlgebra::addFourTensors(vector<vector<vector<vector<double> > > > fourTensor1, vector<vector<vector<vector<double> > > > fourTensor2)
{
    int rowCount = fourTensor1.size();
    int columnCount = fourTensor1[0].size();
    int layer1Count = fourTensor1[0][0].size();
    int layer2Count = fourTensor1[0][0][0].size();

    vector<vector<vector<vector<double> > > > sumTensor(rowCount, vector<vector<vector<double> > >(columnCount, vector<vector<double> >(layer1Count, vector<double>(layer2Count))));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            sumTensor[i][j] = MatrixAlgebra::addMatrices(fourTensor1[i][j], fourTensor2[i][j]);
        }
    }

    return sumTensor;
}

vector<vector<vector<vector<double> > > > TensorAlgebra::subtractFourTensors(vector<vector<vector<vector<double> > > > fourTensor1, vector<vector<vector<vector<double> > > > fourTensor2)
{
    return TensorAlgebra::addFourTensors(fourTensor1, TensorAlgebra::multiplyFourTensor(-1.0, fourTensor2));
}

vector<vector<vector<vector<double> > > > TensorAlgebra::multiplyFourTensor(double scalar, vector<vector<vector<vector<double> > > > fourTensor1)
{
    int rowCount = fourTensor1.size();
    int columnCount = fourTensor1[0].size();
    int layer1Count = fourTensor1[0][0].size();
    int layer2Count = fourTensor1[0][0][0].size();

    vector<vector<vector<vector<double> > > > productTensor(rowCount, vector<vector<vector<double> > >(columnCount, vector<vector<double> >(layer1Count, vector<double>(layer2Count))));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            productTensor[i][j] = MatrixAlgebra::multiplyMatrix(scalar, fourTensor1[i][j]);
        }
    }

    return productTensor;
}

double TensorAlgebra::computeSigmaNorm(vector<vector<double> > tensor1)
{
    double coefficient1 = ((tensor1[0][0] - tensor1[1][1]) * (tensor1[0][0] - tensor1[1][1])) + ((tensor1[1][1] - tensor1[2][2]) * (tensor1[1][1] - tensor1[2][2])) +
            ((tensor1[2][2] - tensor1[0][0]) * (tensor1[2][2] - tensor1[0][0]));
    double coefficient2 = (tensor1[0][1] * tensor1[0][1]) + (tensor1[1][2] * tensor1[1][2]) + (tensor1[2][0] * tensor1[2][0]);

    return sqrt((0.5 * coefficient1) + (3.0 * coefficient2));
}
