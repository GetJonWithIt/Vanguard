#include "hprfirstordersolver.h"

HPRFirstOrderSolver::HPRFirstOrderSolver()
{
}

vector<double> HPRFirstOrderSolver::computeXLaxFriedrichsFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep,
                                                              HPRMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> HPRFirstOrderSolver::computeYLaxFriedrichsFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                              HPRMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> HPRFirstOrderSolver::computeXRichtmyerFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep,
                                                          HPRMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return HPRStateVector::computeXFluxVector(intermediateStateVector, materialParameters);
}

vector<double> HPRFirstOrderSolver::computeYRichtmyerFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                          HPRMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return HPRStateVector::computeYFluxVector(intermediateStateVector, materialParameters);
}

vector<double> HPRFirstOrderSolver::computeXFORCEFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> HPRFirstOrderSolver::computeYFORCEFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

void HPRFirstOrderSolver::computeFORCETimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                               HPRMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, materialParameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, materialParameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
    }
}

void HPRFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  HPRMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);

            vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep, materialParameters);
            vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2], cellSpacing, timeStep, materialParameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
        }
    }
}

void HPRFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  HPRMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);

            vector<double> topFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep, materialParameters);
            vector<double> bottomFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1], cellSpacing, timeStep, materialParameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep), materialParameters);
        }
    }
}

vector<HPRStateVector> HPRFirstOrderSolver::solve(vector<HPRStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                  HPRMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<HPRStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<HPRStateVector> currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = HPRSolvers::computeStableTimeStep(currentCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<HPRStateVector> > HPRFirstOrderSolver::solve2D(vector<vector<HPRStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                             HPRMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<HPRStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<HPRStateVector> > currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = HPRSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);
        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}