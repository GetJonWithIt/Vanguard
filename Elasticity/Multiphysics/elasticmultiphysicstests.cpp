#include "elasticmultiphysicstests.h"

ElasticMultiphysicsTests::ElasticMultiphysicsTests()
{
}

void ElasticMultiphysicsTests::solveBartonTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<ElasticMultiphysicsStateVector> initialCells(cellCount);
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
            initialCells[i] = ElasticMultiphysicsStateVector(0.999, 0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, materialParameters, materialParameters);
        }
        else
        {
            initialCells[i] = ElasticMultiphysicsStateVector(0.001, 0.0, 0.0, 0.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, materialParameters, materialParameters);
        }
    }

    outputSolution(ElasticSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void ElasticMultiphysicsTests::solveBartonTest2(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<ElasticMultiphysicsStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters materialParameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 10.0, 3.0, 2.0);

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
            initialCells[i] = ElasticMultiphysicsStateVector(0.999, 2.0, 0.0, 0.1, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, materialParameters, materialParameters);
        }
        else
        {
            initialCells[i] = ElasticMultiphysicsStateVector(0.001, 0.0, -0.03, -0.01, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, materialParameters, materialParameters);
        }
    }

    outputSolution(ElasticSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void ElasticMultiphysicsTests::solveZhangTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<ElasticMultiphysicsStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters material1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HyperelasticMaterialParameters material2Parameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationTensor(3, vector<double>(3));
    leftDeformationTensor[0][0] = 1.0;
    leftDeformationTensor[0][1] = 0.0;
    leftDeformationTensor[0][2] = 0.0;
    leftDeformationTensor[1][0] = -0.01;
    leftDeformationTensor[1][1] = 0.95;
    leftDeformationTensor[1][2] = 0.02;
    leftDeformationTensor[2][0]= -0.015;
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
            initialCells[i] = ElasticMultiphysicsStateVector(0.999, 2.0, 0.0, 0.1, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, material1Parameters, material2Parameters);
        }
        else
        {
            initialCells[i] = ElasticMultiphysicsStateVector(0.001, 0.0, -0.03, -0.01, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, material1Parameters, material2Parameters);
        }
    }

    outputSolution(ElasticSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.05, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void ElasticMultiphysicsTests::outputSolution(vector<ElasticMultiphysicsStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream volumeFractionFile("volumeFraction.dat");
    ofstream densityFile("density.dat");
    ofstream xVelocityFile("xVelocity.dat");

    for (int i = 0; i < cellCount; i++)
    {
        volumeFractionFile << (cellSpacing * i) << " " << solution[i].getMaterial1VolumeFraction() << endl;
        densityFile << (cellSpacing * i) << " " << solution[i].computeTotalDensity() << endl;
        xVelocityFile << (cellSpacing * i) << " " << solution[i].getMaterial2Entropy() << endl;
    }

    volumeFractionFile.close();
    densityFile.close();
    xVelocityFile.close();
}
