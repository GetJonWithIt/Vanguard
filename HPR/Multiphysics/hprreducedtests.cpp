#include "hprreducedtests.h"

HPRReducedTests::HPRReducedTests()
{
}

void HPRReducedTests::solveBartonTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<HPRReducedStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters hyperelasticMaterialParameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1, pow(10.0, 8.0),
                                             0.0, 0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);
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
            initialCells[i] = HPRReducedStateVector(0.999, 0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, hyperelasticMaterialParameters,
                                                    hyperelasticMaterialParameters, materialParameters, materialParameters);
        }
        else
        {
            initialCells[i] = HPRReducedStateVector(0.001, 0.0, 0.0, 0.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, hyperelasticMaterialParameters,
                                                    hyperelasticMaterialParameters, materialParameters, materialParameters);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void HPRReducedTests::solveZhangTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<HPRReducedStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters hyperelasticMaterial1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters material1Parameters("GodunovRomenski", false, false, 8.31445985, 2.71, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16,
                                              pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.338, 2.0, 1.0, 3.577, 2.088);
    material1Parameters.configureThermalImpulseRelaxationTime(4.0 * pow(10.0, -8.0));

    HyperelasticMaterialParameters hyperelasticMaterial2Parameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters material2Parameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1,
                                              pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);

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
            initialCells[i] = HPRReducedStateVector(0.999, 2.0, 0.0, 0.1, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                    hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
        else
        {
            initialCells[i] = HPRReducedStateVector(0.001, 0.0, -0.03, -0.01, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                    hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.05, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void HPRReducedTests::solveZhangElasticPlasticTest(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<HPRReducedStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters hyperelasticMaterialParameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, true, 8.31445985, 2.71, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16,
                                             pow(10.0, -9.0), 0.0, 24.8 * pow(10.0, 9.0), 0.75, 0.2976 * pow(10.0, 9.0), 1.338, 2.0, 1.0, 3.577, 2.088);

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = HPRReducedStateVector(0.999, 0.75, 0.0, 0.0, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterialParameters,
                                                    hyperelasticMaterialParameters, materialParameters, materialParameters);
        }
        else
        {
            initialCells[i] = HPRReducedStateVector(0.001, -0.75, 0.0, 0.0, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterialParameters,
                                                    hyperelasticMaterialParameters, materialParameters, materialParameters);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, subcyclingIterations, reinitialisationFrequency, materialParameters, materialParameters));
}

void HPRReducedTests::solve2DBartonTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<HPRReducedStateVector> > initialCells(cellCount, vector<HPRReducedStateVector>(cellCount));
    HyperelasticMaterialParameters hyperelasticMaterialParameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1, pow(10.0, 8.0),
                                             0.0, 0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);
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
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = HPRReducedStateVector(0.999, 0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, hyperelasticMaterialParameters,
                                                           hyperelasticMaterialParameters, materialParameters, materialParameters);
            }
            else
            {
                initialCells[i][j] = HPRReducedStateVector(0.001, 0.0, 0.0, 0.0, leftDistortionTensor, pow(10.0, -3.0), rightDistortionTensor, 0.0, hyperelasticMaterialParameters,
                                                           hyperelasticMaterialParameters, materialParameters, materialParameters);
            }
        }
    }

    outputSolution2D(HPRSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.03, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void HPRReducedTests::solve2DZhangTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<HPRReducedStateVector> > initialCells(cellCount, vector<HPRReducedStateVector>(cellCount));
    HyperelasticMaterialParameters hyperelasticMaterial1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters material1Parameters("GodunovRomenski", false, false, 8.31445985, 2.71, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16,
                                              pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.338, 2.0, 1.0, 3.577, 2.088);
    material1Parameters.configureThermalImpulseRelaxationTime(4.0 * pow(10.0, -8.0));

    HyperelasticMaterialParameters hyperelasticMaterial2Parameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters material2Parameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1,
                                              pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);

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
                initialCells[i][j] = HPRReducedStateVector(0.999, 2.0, 0.0, 0.1, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                           hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
            }
            else
            {
                initialCells[i][j] = HPRReducedStateVector(0.001, 0.0, -0.03, -0.01, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                           hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
            }
        }
    }

    outputSolution2D(HPRSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.025, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void HPRReducedTests::solve2DZhangTest2(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<HPRReducedStateVector> > initialCells(cellCount, vector<HPRReducedStateVector>(cellCount));
    HyperelasticMaterialParameters hyperelasticMaterial1Parameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters material1Parameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1,
                                              pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);

    HyperelasticMaterialParameters hyperelasticMaterial2Parameters(0.01, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters material2Parameters("GodunovRomenski", false, false, 8.31445985, 0.01, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1,
                                              pow(10.0, 8.0), 0.0, 0.0, 0.75, 0.0, 1.0, 0.01 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);

    vector<vector<double> > distortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            initialCells[i][j] = HPRReducedStateVector(0.001, 0.0, 0.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                       hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    for (int i = (0.4 * cellCount); i < (0.5 * cellCount); i++)
    {
        for (int j = (0.45 * cellCount); j < (0.55 * cellCount); j++)
        {
            initialCells[i][j] = HPRReducedStateVector(0.999, 0.0, 800.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                       hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    for (int i = (0.5 * cellCount); i < (0.6 * cellCount); i++)
    {
        for (int j = (0.25 * cellCount); j < (0.75 * cellCount); j++)
        {
            initialCells[i][j] = HPRReducedStateVector(0.999, 0.0, 0.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                       hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    outputSolution2D(HPRSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 2.0 * pow(10.0, -4.0), 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void HPRReducedTests::solve2DZhangElasticPlasticTest(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<HPRReducedStateVector> > initialCells(cellCount, vector<HPRReducedStateVector>(cellCount));
    HyperelasticMaterialParameters hyperelasticMaterialParameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, true, 8.31445985, 2.71, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16,
                                              pow(10.0, -9.0), 0.0, 24.8 * pow(10.0, 9.0), 0.75, 0.2976 * pow(10.0, 9.0), 1.338, 2.0, 1.0, 3.577, 2.088);

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 *  cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = HPRReducedStateVector(0.999, 0.75, 0.0, 0.0, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterialParameters,
                                                           hyperelasticMaterialParameters, materialParameters, materialParameters);
            }
            else
            {
                initialCells[i][j] = HPRReducedStateVector(0.001, 0.0, -0.03, -0.01, leftDistortionTensor, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterialParameters,
                                                           hyperelasticMaterialParameters, materialParameters, materialParameters);
            }
        }
    }

    outputSolution2D(HPRSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.03, 0.0, 0, subcyclingIterations, reinitialisationFrequency, materialParameters, materialParameters));
}

void HPRReducedTests::solve2DZhangElasticPlasticTest2(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<HPRReducedStateVector> > initialCells(cellCount, vector<HPRReducedStateVector>(cellCount));
    HyperelasticMaterialParameters hyperelasticMaterial1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters material1Parameters("GodunovRomenski", false, true, 8.31445985, 2.71, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16,
                                              pow(10.0, -9.0), 0.0, 24.8 * pow(10.0, 9.0), 0.75, 0.2976 * pow(10.0, 9.0), 1.338, 2.0, 1.0, 3.577, 2.088);

    HyperelasticMaterialParameters hyperelasticMaterial2Parameters(0.01, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters material2Parameters("GodunovRomenski", false, true, 8.31445985, 0.01, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16,
                                              pow(10.0, -9.0), 0.0, 24.8 * pow(10.0, 9.0), 0.75, 0.2976 * pow(10.0, 9.0), 1.338, 2.0, 1.0, 3.577, 2.088);

    vector<vector<double> > distortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            initialCells[i][j] = HPRReducedStateVector(0.001, 0.0, 0.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                       hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    for (int i = (0.4 * cellCount); i < (0.5 * cellCount); i++)
    {
        for (int j = (0.45 * cellCount); j < (0.55 * cellCount); j++)
        {
            initialCells[i][j] = HPRReducedStateVector(0.999, 0.0, 800.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                       hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    for (int i = (0.5 * cellCount); i < (0.6 * cellCount); i++)
    {
        for (int j = (0.25 * cellCount); j < (0.75 * cellCount); j++)
        {
            initialCells[i][j] = HPRReducedStateVector(0.999, 0.0, 0.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                       hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    outputSolution2D(HPRSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 3.0 * pow(10.0, -5.0), 0.0, 0, subcyclingIterations, reinitialisationFrequency, material1Parameters,
                                                   material2Parameters));
}

void HPRReducedTests::solve2DHowellBallTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    /*
    double cellSpacing = 1.0 / cellCount;

    vector<vector<HPRReducedStateVector> > initialCells(cellCount, vector<HPRReducedStateVector>(cellCount));

    vector<vector<double> > distortionTensor = MatrixAlgebra::computeIDentityMatrix(3);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            initialCells[i][j] = HPRReducedStateVector(0.001, 0.0, 0.0, 0.0, distortionTensor, 0.0, distortionTensor, 0.0, hyperelasticMaterial1Parameters,
                                                       hyperelasticMaterial2Parameters, material1Parameters, material2Parameters);
        }
    }

    for (int i = (0.1 * cellCount); i < (0.15 * cellCount); i++)
    {
        for (int j = (0.44 * cellCount); j < (0.56 * cellCount); j++)
        {
            initialCells[i][j] = HPRReducedStateVector(0.999, 0.0, )
        }
    }
    */
}

void HPRReducedTests::outputSolution(vector<HPRReducedStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream volumeFractionFile("volumeFraction.dat");
    ofstream densityFile("density.dat");
    ofstream yVelocityFile("yVelocity.dat");

    for (int i = 0; i < cellCount; i++)
    {
        volumeFractionFile << (cellSpacing * i) << " " << solution[i].getMaterial1VolumeFraction() << endl;
        densityFile << (cellSpacing * i) << " " << solution[i].computeTotalDensity() << endl;
        yVelocityFile << (cellSpacing * i) << " " << solution[i].getInterfaceYVelocity() << endl;
    }

    volumeFractionFile.close();
    densityFile.close();
    yVelocityFile.close();
}

void HPRReducedTests::outputSolution2D(vector<vector<HPRReducedStateVector> > solution)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream volumeFractionFile("volumeFraction.dat");
    ofstream densityFile("density.dat");
    ofstream yVelocityFile("yVelocity.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            volumeFractionFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getMaterial1VolumeFraction() << endl;
            densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].computeTotalDensity() << endl;
            yVelocityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getInterfaceYVelocity() << endl;
        }
    }

    volumeFractionFile.close();
    densityFile.close();
    yVelocityFile.close();
}
