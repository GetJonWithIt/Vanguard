#include "hprtests.h"

HPRTests::HPRTests()
{
}

void HPRTests::solveStokesFirstProblem(int cellCount, int subcyclingIterations)
{
    vector<HPRStateVector> initialCells(cellCount);
    HPRMaterialParameters materialParameters("StiffenedGas", false, false, 8.31445985, 1.0, 1.0 / 1.4, 0.0, 0.0, 1.0, 1.4, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, pow(10.0, -4.0), 0.75, 0.0, 1.0,
                                             pow(10.0, -16.0), 0.0, 0.0, 0.0);

    materialParameters.configureStrainDissipationTime();
    materialParameters.configureThermalImpulseRelaxationTime();

    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = HPRStateVector(1.0, 1.0 / 1.4, MatrixAlgebra::computeIdentityMatrix(3), 0.0, -0.1, 0.0, 0.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = HPRStateVector(1.0, 1.0 / 1.4, MatrixAlgebra::computeIdentityMatrix(3), 0.0, 0.1, 0.0, 0.0, 0.0, 0.0);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 1.0, 0.0, 0, subcyclingIterations, materialParameters), materialParameters);
}

void HPRTests::solveHeatConductionProblem(int cellCount, int subcyclingIterations)
{
    vector<HPRStateVector> initialCells(cellCount);
    HPRMaterialParameters materialParameters("StiffenedGas", true, false, 8.31445985, 1.0, 1.0, 0.0, 1.0, 2.5, 1.4, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, pow(10.0, -2.0), 0.75, 0.0, 1.0, 2.0,
                                             0.0, 0.0, 0.0);

    materialParameters.configureStrainDissipationTime();
    materialParameters.configureThermalImpulseRelaxationTime(pow(10.0, -2.0));

    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = HPRStateVector(2.0, 1.0, MatrixAlgebra::multiplyMatrix(pow(2.0, (1.0 / 3.0)), MatrixAlgebra::computeIdentityMatrix(3)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = HPRStateVector(0.5, 1.0, MatrixAlgebra::multiplyMatrix(pow(0.5, (1.0 / 3.0)), MatrixAlgebra::computeIdentityMatrix(3)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 1.0, 0.0, 0, subcyclingIterations, materialParameters), materialParameters);
}

void HPRTests::solveBartonTest1(int cellCount, int subcyclingIterations)
{
    vector<HPRStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters hyperelasticMaterialParameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1, pow(10.0, 8.0), 0.0,
                                             0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);
    materialParameters.configureThermalImpulseRelaxationTime(4.0 * pow(10.0, -8.0));

    double cellSpacing = 1.0 / cellCount;

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
            initialCells[i] = HPRStateVector(0.0, 0.5, 1.0, leftDistortionTensor, 0.001, hyperelasticMaterialParameters, materialParameters);
        }
        else
        {
            initialCells[i] = HPRStateVector(0.0, 0.0, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterialParameters, materialParameters);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, subcyclingIterations, materialParameters), materialParameters);
}

void HPRTests::solveZhangElasticPlasticTest(int cellCount, int subcyclingIterations)
{
    vector<HPRStateVector> initialCells(cellCount);
    HyperelasticMaterialParameters hyperelasticMaterialParameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, true, 8.31445985, 2.71, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16, pow(10.0, -9.0),
                                             0.0, 24.8 * pow(10.0, 9.0), 0.75, 0.2976 * pow(10.0, 8.0), 1.338, 2.0, 1.0, 3.577, 2.088);

    double cellSpacing = 1.0 / cellCount;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = HPRStateVector(0.75, 0.0, 0.0, leftDistortionTensor, 0.0, hyperelasticMaterialParameters, materialParameters);
        }
        else
        {
            initialCells[i] = HPRStateVector(-0.75, 0.0, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterialParameters, materialParameters);
        }
    }

    outputSolution(HPRSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.06, 0.0, 0, subcyclingIterations, materialParameters), materialParameters);
}

void HPRTests::solve2DBartonTest1(int cellCount, int subcyclingIterations)
{
    vector<vector<HPRStateVector> > initialCells(cellCount, vector<HPRStateVector>(cellCount));
    HyperelasticMaterialParameters hyperelasticMaterialParameters(8.9, 4.6, 2.1, 3.94 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, false, 8.31445985, 8.9, 0.0, 300.0, 300.0, 3.94 * pow(10.0, -4.0), 1.0, 0.0, 3.909, 0.0, 0.0, 0.0, 2.1, pow(10.0, 8.0),
                                             0.0, 0.0, 0.75, 0.0, 1.0, 8.9 * 3.909 * sqrt(4.0 * pow(10.0, -4.0) / 300.0), 1.0, 3.0, 2.0);
    materialParameters.configureThermalImpulseRelaxationTime(4.0 * pow(10.0, -8.0));

    double cellSpacing = 1.0 / cellCount;

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
                initialCells[i][j] = HPRStateVector(0.0, 0.5, 1.0, leftDistortionTensor, 0.001, hyperelasticMaterialParameters, materialParameters);
            }
            else
            {
                initialCells[i][j] = HPRStateVector(0.0, 0.0, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterialParameters, materialParameters);
            }
        }
    }

    outputSolution2D(HPRSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.03, 0.0, 0, subcyclingIterations, materialParameters), materialParameters);
}

void HPRTests::solve2DZhangElasticPlasticTest(int cellCount, int subcyclingIterations)
{
    vector<vector<HPRStateVector> > initialCells(cellCount, vector<HPRStateVector>(cellCount));
    HyperelasticMaterialParameters hyperelasticMaterialParameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HPRMaterialParameters materialParameters("GodunovRomenski", false, true, 8.31445985, 2.71, 73.0, 300.0, 300.0, 9.0 * pow(10.0, -4.0), 1.0, 0.0, 5.037, 1.0, 1.338, 0.0, 3.16,
                                             pow(10.0, -9.0), 0.0, 24.8 * pow(10.0, 9.0), 0.75, 0.2976 * pow(10.0, 9.0), 1.338, 2.0, 1.0, 3.577, 2.088);

    double cellSpacing = 1.0 / cellCount;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = HPRStateVector(0.75, 0.0, 0.0, leftDistortionTensor, 0.0, hyperelasticMaterialParameters, materialParameters);
            }
            else
            {
                initialCells[i][j] = HPRStateVector(-0.75, 0.0, 0.0, rightDistortionTensor, 0.0, hyperelasticMaterialParameters, materialParameters);
            }
        }
    }

    outputSolution2D(HPRSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.03, 0.0, 0, subcyclingIterations, materialParameters), materialParameters);
}

void HPRTests::outputSolution(vector<HPRStateVector> solution, HPRMaterialParameters materialParameters)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream densityFile("density.dat");
    ofstream xVelocityFile("xVelocity.dat");
    ofstream yVelocityFile("yVelocity.dat");

    ofstream temperatureFile("temperature.dat");
    ofstream heatFluxFile("heatFlux.dat");

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << solution[i].getDensity() << endl;
        xVelocityFile << (cellSpacing * i) << " " << solution[i].getXVelocity() << endl;
        yVelocityFile << (cellSpacing * i) << " " << solution[i].getYVelocity() << endl;

        temperatureFile << (cellSpacing * i) << " " << solution[i].computeTemperature(materialParameters) << endl;
        heatFluxFile << (cellSpacing * i) << " " << solution[i].computeHeatFluxVector(materialParameters)[0] << endl;
    }

    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();

    temperatureFile.close();
    heatFluxFile.close();
}

void HPRTests::outputSolution2D(vector<vector<HPRStateVector> > solution, HPRMaterialParameters materialParameters)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream densityFile("density.dat");
    ofstream xVelocityFile("xVelocity.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getDensity() << endl;
            xVelocityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getXVelocity() << endl;
        }
    }

    densityFile.close();
    xVelocityFile.close();
}
