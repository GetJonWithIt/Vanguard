#include "hprsecondordersolver.h"

HPRSecondOrderSolver::HPRSecondOrderSolver()
{
}

vector<double> HPRSecondOrderSolver::computeXSLICFlux(HPRStateVector leftLeftStateVector, HPRStateVector leftStateVector, HPRStateVector rightStateVector, HPRStateVector rightRightStateVector,
                                                      double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    HPRStateVector evolvedRightStateVector = HPRSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                    materialParameters);
    HPRStateVector evolvedLeftStateVector = HPRSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                   materialParameters);

    return HPRFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> HPRSecondOrderSolver::computeYSLICFlux(HPRStateVector topTopStateVector, HPRStateVector topStateVector, HPRStateVector bottomStateVector, HPRStateVector bottomBottomStateVector,
                                                      double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    HPRStateVector evolvedBottomStateVector = HPRSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                     materialParameters);
    HPRStateVector evolvedTopStateVector = HPRSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                  materialParameters);

    return HPRFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, materialParameters);
}

void HPRSecondOrderSolver::computeSLICTimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                               int slopeLimiter, HPRMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

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

void HPRSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

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

void HPRSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

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

vector<HPRStateVector> HPRSecondOrderSolver::solve(vector<HPRStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                   int subcyclingIterations, HPRMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<HPRStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<HPRStateVector> currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = HPRSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, (timeStep / subcyclingIterations), bias, slopeLimiter, materialParameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<HPRStateVector> > HPRSecondOrderSolver::solve2D(vector<vector<HPRStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                              int slopeLimiter, int subcyclingIterations, HPRMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<HPRStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<HPRStateVector> > currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = HPRSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, (timeStep / subcyclingIterations), bias, slopeLimiter, materialParameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}
