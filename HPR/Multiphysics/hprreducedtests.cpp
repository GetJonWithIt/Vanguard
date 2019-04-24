#include "hprreducedtests.h"

HPRReducedTests::HPRReducedTests()
{
}

void HPRReducedTests::solveBartonTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<HPRReducedStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters hyperelasticMaterialParameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1, pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);
    materialParameters.configureThermalImpulseRelaxationTime(4.0 * pow(10.0, -8.0));

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
            initialCells[i] = HPRReducedStateVector(0.999, 0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, hyperelasticMaterialParameters, hyperelasticMaterialParameters, materialParameters, materialParameters);
        }
        else
        {
            initialCells[i] = HPRReducedStateVector(0.001, 0.0, 0.0, 0.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, hyperelasticMaterialParameters, hyperelasticMaterialParameters, materialParameters, materialParameters);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void HPRReducedTests::solveZhangTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<HPRReducedStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters hyperelasticMaterial1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters material1Parameters("GodunovRomenski", false, false, 8.31445985, 2.71, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16, pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.338, 2.0, 1.0, 3.577, 2.088);
    material1Parameters.configureThermalImpulseRelaxationTime(4.0 * pow(10.0, -8.0));

    HyperelasticMaterialParameters hyperelasticMaterial2Parameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters material2Parameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1, pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);

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
            initialCells[i] = HPRReducedStateVector(0.999, 2.0, 0.0, 0.1, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterial1Parameters, hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
        else
        {
            initialCells[i] = HPRReducedStateVector(0.001, 0.0, -0.03, -0.01, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterial1Parameters, hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.05, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void HPRReducedTests::outputSolution(vector<HPRReducedStateVector> solution)
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
