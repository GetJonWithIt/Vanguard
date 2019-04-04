#include "matrixalgebra.h"

MatrixAlgebra::MatrixAlgebra()
{
}

vector<vector<double> > MatrixAlgebra::computeIdentityMatrix(int dimension)
{
    vector<vector<double> > identityMatrix(dimension, vector<double>(dimension));

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            if (i == j)
            {
                identityMatrix[i][j] = 1.0;
            }
            else
            {
                identityMatrix[i][j] = 0.0;
            }
        }
    }

    return identityMatrix;
}

vector<vector<double> > MatrixAlgebra::addMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<double> > sumMatrix(rowCount, vector<double>(columnCount));

    for (int i = 0; i < rowCount; i++)
    {
        sumMatrix[i] = VectorAlgebra::addVectors(matrix1[i], matrix2[i]);
    }

    return sumMatrix;
}

vector<vector<double> > MatrixAlgebra::subtractMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    return addMatrices(matrix1, multiplyMatrix(-1.0, matrix2));
}

vector<vector<double> > MatrixAlgebra::multiplyMatrix(double scalar1, vector<vector<double> > matrix1)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<double> > productMatrix(rowCount, vector<double>(columnCount));

    for (int i = 0; i < rowCount; i++)
    {
        productMatrix[i] = VectorAlgebra::multiplyVector(scalar1, matrix1[i]);
    }

    return productMatrix;
}

vector<vector<double> > MatrixAlgebra::multiplyMatrices(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    int matrix1RowCount = matrix1.size();
    int matrix1ColumnCount = matrix1[0].size();
    int matrix2ColumnCount = matrix2[0].size();
    vector<vector<double> > productMatrix(matrix1RowCount, vector<double>(matrix2ColumnCount));

    for (int i = 0; i < matrix1RowCount; i++)
    {
        for (int j = 0; j < matrix2ColumnCount; j++)
        {
            for (int k = 0; k < matrix1ColumnCount; k++)
            {
                productMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return productMatrix;
}

vector<vector<double> > MatrixAlgebra::multiplyMatricesElementWise(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<double> > productMatrix(rowCount, vector<double>(columnCount));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            productMatrix[i][j] = matrix1[i][j] * matrix2[i][j];
        }
    }

    return productMatrix;
}

vector<double> MatrixAlgebra::multiplyMatrixByVector(vector<vector<double> > matrix1, vector<double> vector1)
{
    int componentCount = vector1.size();

    vector<vector<double> > matrix2Transpose(1, vector<double>(componentCount));
    matrix2Transpose[0] = vector1;
    vector<vector<double> > matrix2 = computeTranspose(matrix2Transpose);

    vector<vector<double> > productMatrix = multiplyMatrices(matrix1, matrix2);
    return computeTranspose(productMatrix)[0];
}

vector<vector<double> > MatrixAlgebra::computeTranspose(vector<vector<double> > matrix1)
{
    int rowCount = matrix1.size();
    int columnCount = matrix1[0].size();
    vector<vector<double> > transposeMatrix(columnCount, vector<double>(rowCount));

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            transposeMatrix[j][i] = matrix1[i][j];
        }
    }

    return transposeMatrix;
}

double MatrixAlgebra::computeTrace(vector<vector<double> > matrix1)
{
    double trace = 0.0;
    int rowCount = matrix1.size();

    for (int i = 0; i < rowCount; i++)
    {
        trace += matrix1[i][i];
    }

    return trace;
}

double MatrixAlgebra::computeDeterminant(vector<vector<double> > matrix1)
{
    return (matrix1[0][0] * ((matrix1[1][1] * matrix1[2][2]) - (matrix1[2][1] * matrix1[1][2]))) - (matrix1[0][1] * ((matrix1[1][0] * matrix1[2][2]) - (matrix1[1][2] * matrix1[2][0]))) +
            (matrix1[0][2] * ((matrix1[1][0] * matrix1[2][1]) - (matrix1[1][1] * matrix1[2][0])));
}

double MatrixAlgebra::computeComponentSum(vector<vector<double> > matrix1)
{
    double componentSum = 0.0;
    int rowCount = matrix1.size();

    for (int i = 0; i < rowCount; i++)
    {
        componentSum += VectorAlgebra::computeComponentSum(matrix1[i]);
    }

    return componentSum;
}

vector<vector<double> > MatrixAlgebra::computeInverseMatrix(vector<vector<double> > matrix1)
{
    int rowCount = matrix1.size();
    vector<vector<double> > inverseMatrix(rowCount, vector<double>(rowCount));
    double determinant = computeDeterminant(matrix1);

    inverseMatrix[0][0] = ((matrix1[1][1] * matrix1[2][2]) - (matrix1[2][1] * matrix1[1][2])) / determinant;
    inverseMatrix[0][1] = ((matrix1[0][2] * matrix1[2][1]) - (matrix1[0][1] * matrix1[2][2])) / determinant;
    inverseMatrix[0][2] = ((matrix1[0][1] * matrix1[1][2]) - (matrix1[0][2] * matrix1[1][1])) / determinant;

    inverseMatrix[1][0] = ((matrix1[1][2] * matrix1[2][0]) - (matrix1[1][0] * matrix1[2][2])) / determinant;
    inverseMatrix[1][1] = ((matrix1[0][0] * matrix1[2][2]) - (matrix1[0][2] * matrix1[2][0])) / determinant;
    inverseMatrix[1][2] = ((matrix1[1][0] * matrix1[0][2]) - (matrix1[0][0] * matrix1[1][2])) / determinant;

    inverseMatrix[2][0] = ((matrix1[1][0] * matrix1[2][1]) - (matrix1[2][0] * matrix1[1][1])) / determinant;
    inverseMatrix[2][1] = ((matrix1[2][0] * matrix1[0][1]) - (matrix1[0][0] * matrix1[2][1])) / determinant;
    inverseMatrix[2][2] = ((matrix1[0][0] * matrix1[1][1]) - (matrix1[1][0] * matrix1[0][1])) / determinant;

    return inverseMatrix;
}

double MatrixAlgebra::computeLargestEigenvalue(vector<vector<double> > matrix1)
{
    int columnCount = matrix1[0].size();

    vector<double> eigenvectorGuess(columnCount);
    for (int i = 0; i < columnCount; i++)
    {
        eigenvectorGuess[i] = 1.0;
    }

    for (int i = 0; i < 3; i++)
    {
        eigenvectorGuess = multiplyMatrixByVector(matrix1, eigenvectorGuess);
        eigenvectorGuess = VectorAlgebra::multiplyVector((1.0 / VectorAlgebra::computeNorm(eigenvectorGuess)), eigenvectorGuess);
    }

    return VectorAlgebra::computeNorm(MatrixAlgebra::multiplyMatrixByVector(matrix1, eigenvectorGuess));
}
