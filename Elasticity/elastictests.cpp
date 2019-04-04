#include "elastictests.h"

ElasticTests::ElasticTests()
{
}

void ElasticTests::solveBartonTest1(int cellCount)
{
    double cellSpacing = 1.0 / cellCount;

    vector<ElasticStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters materialParameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationTensor(3, vector<double>(3));
    leftDeformationTensor[0][0] = 0.98;
    leftDeformationTensor[0][1] = 0.0;
    leftDeformationTensor[0][2] = 0.0;
    leftDeformationTensor[1][0] = 0.02;
    leftDeformationTensor[1][1] = 1.0;
    leftDeformationTensor[1][2] = 0.1;
    leftDeformationTensor[2][0] = 0.0;
    leftDeformationTensor[2][1] = 0.0;
    leftDeformationTensor[2][2] = 1.0;

    vector<vector<double> > rightDeformationTensor(3, vector<double>(3));
    rightDeformationTensor[0][0] = 1.0;
    rightDeformationTensor[0][1] = 0.0;
    rightDeformationTensor[0][2] = 0.0;
    rightDeformationTensor[1][0] = 0.0;
    rightDeformationTensor[1][1] = 1.0;
    rightDeformationTensor[1][2] = 0.1;
    rightDeformationTensor[2][0] = 0.0;
    rightDeformationTensor[2][1] = 0.0;
    rightDeformationTensor[2][2] = 1.0;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeInverseMatrix(leftDeformationTensor);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeInverseMatrix(rightDeformationTensor);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = ElasticStateVector(0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), materialParameters);
        }
        else
        {
            initialCells[i] = ElasticStateVector(0.0, 0.0, 0.0, rightDistortionTensor, 0.0, materialParameters);
        }
    }

    outputSolution(SecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, 0, materialParameters));
}

void ElasticTests::solveBartonTest2(int cellCount)
{
    double cellSpacing = 1.0 / cellCount;

    vector<ElasticStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters materialParameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationTensor(3, vector<double>(3));
    leftDeformationTensor[0][0] = 1.0;
    leftDeformationTensor[0][1] = 0.0;
    leftDeformationTensor[0][2] = 0.0;
    leftDeformationTensor[1][0] = -0.01;
    leftDeformationTensor[1][1] = 0.95;
    leftDeformationTensor[1][2] = 0.02;
    leftDeformationTensor[2][0] = -0.015;
    leftDeformationTensor[2][1] = 0.0;
    leftDeformationTensor[2][2] = 0.9;

    vector<vector<double> > rightDeformationTensor(3, vector<double>(3));
    rightDeformationTensor[0][0] = 1.0;
    rightDeformationTensor[0][1] = 0.0;
    rightDeformationTensor[0][2] = 0.0;
    rightDeformationTensor[1][0] = 0.015;
    rightDeformationTensor[1][1] = 0.95;
    rightDeformationTensor[1][2] = 0.0;
    rightDeformationTensor[2][0] = -0.01;
    rightDeformationTensor[2][1] = 0.0;
    rightDeformationTensor[2][2] = 0.9;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeInverseMatrix(leftDeformationTensor);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeInverseMatrix(rightDeformationTensor);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = ElasticStateVector(2.0, 0.0, 0.1, leftDistortionTensor, 0.0, materialParameters);
        }
        else
        {
            initialCells[i] = ElasticStateVector(0.0, 0.0, 0.0, rightDistortionTensor, 0.0, materialParameters);
        }
    }

    outputSolution(SecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, 0, materialParameters));
}

void ElasticTests::solve2DBartonTest1(int cellCount)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<ElasticStateVector> > initialCells(cellCount, vector<ElasticStateVector>(cellCount));
    HyperelasticMaterialParameters materialParameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationTensor(3, vector<double>(3));
    leftDeformationTensor[0][0] = 0.98;
    leftDeformationTensor[0][1] = 0.0;
    leftDeformationTensor[0][2] = 0.0;
    leftDeformationTensor[1][0] = 0.02;
    leftDeformationTensor[1][1] = 1.0;
    leftDeformationTensor[1][2] = 0.1;
    leftDeformationTensor[2][0] = 0.0;
    leftDeformationTensor[2][1] = 0.0;
    leftDeformationTensor[2][2] = 1.0;

    vector<vector<double> > rightDeformationTensor(3, vector<double>(3));
    rightDeformationTensor[0][0] = 1.0;
    rightDeformationTensor[0][1] = 0.0;
    rightDeformationTensor[0][2] = 0.0;
    rightDeformationTensor[1][0] = 0.0;
    rightDeformationTensor[1][1] = 1.0;
    rightDeformationTensor[1][2] = 0.1;
    rightDeformationTensor[2][0] = 0.0;
    rightDeformationTensor[2][1] = 0.0;
    rightDeformationTensor[2][2] = 1.0;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeInverseMatrix(leftDeformationTensor);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeInverseMatrix(rightDeformationTensor);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = ElasticStateVector(0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), materialParameters);
            }
            else
            {
                initialCells[i][j] = ElasticStateVector(0.0, 0.0, 0.0, rightDistortionTensor, 0.0, materialParameters);
            }
        }
    }

    outputSolution2D(SecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.03, 0.0, 0, 0, materialParameters));
}

void ElasticTests::solve2DBartonTest2(int cellCount)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<ElasticStateVector> > initialCells(cellCount, vector<ElasticStateVector>(cellCount));
    HyperelasticMaterialParameters materialParameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationTensor(3, vector<double>(3));
    leftDeformationTensor[0][0] = 1.0;
    leftDeformationTensor[0][1] = 0.0;
    leftDeformationTensor[0][2] = 0.0;
    leftDeformationTensor[1][0] = -0.01;
    leftDeformationTensor[1][1] = 0.95;
    leftDeformationTensor[1][2] = 0.02;
    leftDeformationTensor[2][0] = -0.015;
    leftDeformationTensor[2][1] = 0.0;
    leftDeformationTensor[2][2] = 0.9;

    vector<vector<double> > rightDeformationTensor(3, vector<double>(3));
    rightDeformationTensor[0][0] = 1.0;
    rightDeformationTensor[0][1] = 0.0;
    rightDeformationTensor[0][2] = 0.0;
    rightDeformationTensor[1][0] = 0.015;
    rightDeformationTensor[1][1] = 0.95;
    rightDeformationTensor[1][2] = 0.0;
    rightDeformationTensor[2][0] = -0.01;
    rightDeformationTensor[2][1] = 0.0;
    rightDeformationTensor[2][2] = 0.9;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeInverseMatrix(leftDeformationTensor);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeInverseMatrix(rightDeformationTensor);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = ElasticStateVector(2.0, 0.0, 0.1, leftDistortionTensor, 0.0, materialParameters);
            }
            else
            {
                initialCells[i][j] = ElasticStateVector(0.0, 0.0, 0.0, rightDistortionTensor, 0.0, materialParameters);
            }
        }
    }

    outputSolution2D(SecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.03, 0.0, 0, 0, materialParameters));
}

void ElasticTests::outputSolution(vector<ElasticStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream densityFile("density.dat");

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << solution[i].getDensity() << endl;
    }

    densityFile.close();
}

void ElasticTests::outputSolution2D(vector<vector<ElasticStateVector> > solution)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream densityFile("density.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getDensity() << endl;
        }
    }

    densityFile.close();
}
