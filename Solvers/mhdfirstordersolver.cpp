#include "mhdfirstordersolver.h"

MHDFirstOrderSolver::MHDFirstOrderSolver()
{
}

vector<double> MHDFirstOrderSolver::computeXLaxFriedrichsFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep,
                                                              MHDMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> MHDFirstOrderSolver::computeXRichtmyerFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep,
                                                          MHDMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing,
                                                                                    timeStep);
    return MHDStateVector::computeXFluxVector(intermediateStateVector, materialParameters);
}

vector<double> MHDFirstOrderSolver::computeXFORCEFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep,
                                                      MHDMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

void MHDFirstOrderSolver::computeFORCETimeStep(vector<MHDStateVector> & currentCells, vector<MHDStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                               MHDMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, materialParameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, materialParameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
    }
}

vector<MHDStateVector> MHDFirstOrderSolver::solve(vector<MHDStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                  MHDMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<MHDStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<MHDStateVector> currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = MHDSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}
