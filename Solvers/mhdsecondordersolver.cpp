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

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}
