#include "elasticreducedtests.h"

ElasticReducedTests::ElasticReducedTests()
{
}

void ElasticReducedTests::solveBartonTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<ElasticReducedStateVector> initialCells(cellCount);
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
            initialCells[i] = ElasticReducedStateVector(0.999, 0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, materialParameters, materialParameters);
        }
        else
        {
            initialCells[i] = ElasticReducedStateVector(0.001, 0.0, 0.0, 0.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, materialParameters, materialParameters);
        }
    }

    outputSolution(ElasticSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.03, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void ElasticReducedTests::solveZhangTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<ElasticReducedStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters material1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HyperelasticMaterialParameters material2Parameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

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
            initialCells[i] = ElasticReducedStateVector(0.999, 2.0, 0.0, 0.1, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, material1Parameters, material2Parameters);
        }
        else
        {
            initialCells[i] = ElasticReducedStateVector(0.001, 0.0, -0.03, -0.01, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, material1Parameters, material2Parameters);
        }
    }

    outputSolution(ElasticSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.025, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void ElasticReducedTests::solve2DBartonTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<ElasticReducedStateVector> > initialCells(cellCount, vector<ElasticReducedStateVector>(cellCount));
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
                initialCells[i][j] = ElasticReducedStateVector(0.999, 0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, materialParameters,
                                                               materialParameters);
            }
            else
            {
                initialCells[i][j] = ElasticReducedStateVector(0.001, 0.0, 0.0, 0.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, materialParameters,
                                                               materialParameters);
            }
        }
    }

    outputSolution2D(ElasticSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.015, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void ElasticReducedTests::solve2DZhangTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<ElasticReducedStateVector> > initialCells(cellCount, vector<ElasticReducedStateVector>(cellCount));
    HyperelasticMaterialParameters material1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HyperelasticMaterialParameters material2Parameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

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
                initialCells[i][j] = ElasticReducedStateVector(0.999, 2.0, 0.0, 0.1, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, material1Parameters, material2Parameters);
            }
            else
            {
                initialCells[i][j] = ElasticReducedStateVector(0.001, 0.0, -0.03, -0.01, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, material1Parameters, material2Parameters);
            }
        }
    }

    outputSolution2D(ElasticSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.0125, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void ElasticReducedTests::solve2DZhangTest2(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<ElasticReducedStateVector> > initialCells(cellCount, vector<ElasticReducedStateVector>(cellCount));
    HyperelasticMaterialParameters material1Parameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HyperelasticMaterialParameters material2Parameters(0.01, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > distortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            initialCells[i][j] = ElasticReducedStateVector(0.001, 0.0, 0.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, material1Parameters, material2Parameters);
        }
    }

    for (int i = (0.4 * cellCount); i < (0.5 * cellCount); i++)
    {
        for (int j = (0.45 * cellCount); j < (0.55 * cellCount); j++)
        {
            initialCells[i][j] = ElasticReducedStateVector(0.999, 0.0, 800.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, material1Parameters, material2Parameters);
        }
    }

    for (int i = (0.5 * cellCount); i < (0.6 * cellCount); i++)
    {
        for (int j = (0.25 * cellCount); j < (0.75 * cellCount); j++)
        {
            initialCells[i][j] = ElasticReducedStateVector(0.999, 0.0, 0.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, material1Parameters, material2Parameters);
        }
    }

    outputSolution2D(ElasticSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 2.3 * pow(10.0, -4.0), 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void ElasticReducedTests::outputSolution(vector<ElasticReducedStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream volumeFractionFile("volumeFraction.dat");
    ofstream densityFile("density.dat");

    for (int i = 0; i < cellCount; i++)
    {
        volumeFractionFile << (cellSpacing * i) << " " << solution[i].getMaterial1VolumeFraction() << endl;
        densityFile << (cellSpacing * i) << " " << solution[i].computeTotalDensity() << endl;
    }

    volumeFractionFile.close();
    densityFile.close();
}

void ElasticReducedTests::outputSolution2D(vector<vector<ElasticReducedStateVector> > solution)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream volumeFractionFile("volumeFraction.dat");
    ofstream densityFile("density.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            volumeFractionFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getMaterial1VolumeFraction() << endl;
            densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].computeTotalDensity() << endl;
        }
    }

    volumeFractionFile.close();
    densityFile.close();
}
