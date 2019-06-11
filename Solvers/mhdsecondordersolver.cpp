#include "mhdsecondordersolver.h"

MHDSecondOrderSolver::MHDSecondOrderSolver()
{
}

vector<double> MHDSecondOrderSolver::computeXSLICFlux(MHDStateVector leftLeftStateVector, MHDStateVector leftStateVector, MHDStateVector rightStateVector,
                                                      MHDStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      MHDMaterialParameters materialParameters)
{
    MHDStateVector evolvedRightStateVector = MHDSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                    1, materialParameters);
    MHDStateVector evolvedLeftStateVector = MHDSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                   0, materialParameters);

    return MHDFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> MHDSecondOrderSolver::computeXSLICFlux(MHDReducedStateVector leftLeftStateVector, MHDReducedStateVector leftStateVector, MHDReducedStateVector rightStateVector,
                                                      MHDReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    MHDReducedStateVector evolvedRightStateVector = MHDSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                           1, material1Parameters, material2Parameters);
    MHDReducedStateVector evolvedLeftStateVector = MHDSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                          0, material1Parameters, material2Parameters);

    return MHDFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> MHDSecondOrderSolver::computeYSLICFlux(MHDStateVector topTopStateVector, MHDStateVector topStateVector, MHDStateVector bottomStateVector,
                                                      MHDStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      MHDMaterialParameters materialParameters)
{
    MHDStateVector evolvedBottomStateVector = MHDSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                     1, materialParameters);
    MHDStateVector evolvedTopStateVector = MHDSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                  0, materialParameters);

    return MHDFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> MHDSecondOrderSolver::computeYSLICFlux(MHDReducedStateVector topTopStateVector, MHDReducedStateVector topStateVector, MHDReducedStateVector bottomStateVector,
                                                      MHDReducedStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    MHDReducedStateVector evolvedBottomStateVector = MHDSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                            1, material1Parameters, material2Parameters);
    MHDReducedStateVector evolvedTopStateVector = MHDSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                         0, material1Parameters, material2Parameters);

    return MHDFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

void MHDSecondOrderSolver::computeSLICTimeStep(vector<MHDStateVector> & currentCells, vector<MHDStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                               int slopeLimiter, MHDMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);

        vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
        vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3], currentCellsWithBoundary[i + 4],
                cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
    }
}

void MHDSecondOrderSolver::computeSLICTimeStep(vector<MHDReducedStateVector> & currentCells, vector<MHDReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                               double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);

        vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3], currentCellsWithBoundary[i + 4],
                cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), material1Parameters,
                                                   material2Parameters);
    }
}

void MHDSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  double bias, int slopeLimiter, MHDMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);

            vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 2][j], currentCellsWithBoundary[i + 2][j + 1], currentCellsWithBoundary[i + 2][j + 2],
                    currentCellsWithBoundary[i + 2][j + 3], cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
            vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 2][j + 1], currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 2][j + 3],
                    currentCellsWithBoundary[i + 2][j + 4], cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                          materialParameters);
        }
    }
}

void MHDSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary,
                                                  double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 2][j], currentCellsWithBoundary[i + 2][j + 1], currentCellsWithBoundary[i + 2][j + 2],
                    currentCellsWithBoundary[i + 2][j + 3], cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
            vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 2][j + 1], currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 2][j + 3],
                    currentCellsWithBoundary[i + 2][j + 4], cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                          material1Parameters, material2Parameters);
        }
    }
}

void MHDSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  double bias, int slopeLimiter, MHDMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);

            vector<double> topFluxVector = computeYSLICFlux(currentCellsWithBoundary[i][j + 2], currentCellsWithBoundary[i + 1][j + 2], currentCellsWithBoundary[i + 2][j + 2],
                    currentCellsWithBoundary[i + 3][j + 2], cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
            vector<double> bottomFluxVector = computeYSLICFlux(currentCellsWithBoundary[i + 1][j + 2], currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 3][j + 2],
                    currentCellsWithBoundary[i + 4][j + 2], cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep),
                                                          materialParameters);
        }
    }
}

void MHDSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary,
                                                  double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> topFluxVector = computeYSLICFlux(currentCellsWithBoundary[i][j + 2], currentCellsWithBoundary[i + 1][j + 2], currentCellsWithBoundary[i + 2][j + 2],
                    currentCellsWithBoundary[i + 3][j + 2], cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
            vector<double> bottomFluxVector = computeYSLICFlux(currentCellsWithBoundary[i + 1][j + 2], currentCellsWithBoundary[i + 2][j + 2], currentCellsWithBoundary[i + 3][j + 2],
                    currentCellsWithBoundary[i + 4][j + 2], cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep),
                                                          material1Parameters, material2Parameters);
        }
    }
}

vector<MHDStateVector> MHDSecondOrderSolver::solve(vector<MHDStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                   int subcyclingIterations, MHDMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<MHDStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<MHDStateVector> currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = MHDSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        double maximumWaveSpeed = MHDSolvers::computeMaximumWaveSpeed(currentCellsWithBoundary, materialParameters);
        materialParameters.setHyperbolicWaveSpeed(maximumWaveSpeed);
        materialParameters.configureParabolicDamping();

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                        materialParameters);
        }

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 2);
        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                        materialParameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<MHDReducedStateVector> MHDSecondOrderSolver::solve(vector<MHDReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                          int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, MHDMaterialParameters material1Parameters,
                                                          MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<MHDReducedStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<MHDReducedStateVector> currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = MHDSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                            material2Parameters);

        double maximumWaveSpeed = MHDSolvers::computeMaximumWaveSpeed(currentCellsWithBoundary, material1Parameters, material2Parameters);
        material1Parameters.setHyperbolicWaveSpeed(maximumWaveSpeed);
        material2Parameters.setHyperbolicWaveSpeed(maximumWaveSpeed);

        material1Parameters.configureParabolicDamping();
        material2Parameters.configureParabolicDamping();

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                        material1Parameters, material2Parameters);
        }

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 2);
        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                        material1Parameters, material2Parameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        if (reinitialisationFrequency != 0 && currentIteration != 0)
        {
            if ((currentIteration % reinitialisationFrequency) == 0)
            {
                MultiphysicsSolvers::reinitialiseVolumeFraction(currentCells, material1Parameters, material2Parameters);
            }
        }

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<MHDStateVector> > MHDSecondOrderSolver::solve2D(vector<vector<MHDStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                              int slopeLimiter, int subcyclingIterations, MHDMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<MHDStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<MHDStateVector> > currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = MHDSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        double maximumWaveSpeed = MHDSolvers::computeMaximumWaveSpeed2D(currentCellsWithBoundary, materialParameters);
        materialParameters.setHyperbolicWaveSpeed(maximumWaveSpeed);
        materialParameters.configureParabolicDamping();

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          materialParameters);
        }


        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, materialParameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          materialParameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<MHDReducedStateVector> > MHDSecondOrderSolver::solve2D(vector<vector<MHDReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                     double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                     MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<MHDReducedStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<MHDReducedStateVector> > currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = MHDSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                              material2Parameters);

        double maximumWaveSpeed = MHDSolvers::computeMaximumWaveSpeed2D(currentCellsWithBoundary, material1Parameters, material2Parameters);
        material1Parameters.setHyperbolicWaveSpeed(maximumWaveSpeed);
        material2Parameters.setHyperbolicWaveSpeed(maximumWaveSpeed);

        material1Parameters.configureParabolicDamping();
        material2Parameters.configureParabolicDamping();

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          material1Parameters, material2Parameters);
        }

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          material1Parameters, material2Parameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        if (reinitialisationFrequency != 0 && currentIteration != 0)
        {
            if ((currentIteration % reinitialisationFrequency) == 0)
            {
                MultiphysicsSolvers::reinitialiseVolumeFraction(currentCells, material1Parameters, material2Parameters);
            }
        }

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}
