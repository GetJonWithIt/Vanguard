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

vector<double> HPRSecondOrderSolver::computeXSLICFlux(HPRMultiphysicsStateVector leftLeftStateVector, HPRMultiphysicsStateVector leftStateVector, HPRMultiphysicsStateVector rightStateVector,
                                                      HPRMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRMultiphysicsStateVector evolvedRightStateVector = HPRSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                1, material1Parameters, material2Parameters);
    HPRMultiphysicsStateVector evolvedLeftStateVector = HPRSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                               0, material1Parameters, material2Parameters);

    return HPRFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> HPRSecondOrderSolver::computeXSLICFlux(HPRIntermediateStateVector leftLeftStateVector, HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector rightStateVector,
                                                      HPRIntermediateStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRIntermediateStateVector evolvedRightStateVector = HPRSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                1, material1Parameters, material2Parameters);
    HPRIntermediateStateVector evolvedLeftStateVector = HPRSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                               0, material1Parameters, material2Parameters);

    return HPRFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> HPRSecondOrderSolver::computeXSLICFlux(HPRReducedStateVector leftLeftStateVector, HPRReducedStateVector leftStateVector, HPRReducedStateVector rightStateVector,
                                                      HPRReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRReducedStateVector evolvedRightStateVector = HPRSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                           material1Parameters, material2Parameters);
    HPRReducedStateVector evolvedLeftStateVector = HPRSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                          material1Parameters, material2Parameters);

    return HPRFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
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

vector<double> HPRSecondOrderSolver::computeYSLICFlux(HPRIntermediateStateVector topTopStateVector, HPRIntermediateStateVector topStateVector, HPRIntermediateStateVector bottomStateVector,
                                                      HPRIntermediateStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRIntermediateStateVector evolvedBottomStateVector = HPRSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias,
                                                                                                 slopeLimiter, 1, material1Parameters, material2Parameters);
    HPRIntermediateStateVector evolvedTopStateVector = HPRSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias,
                                                                                              slopeLimiter, 0, material1Parameters, material2Parameters);

    return HPRFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> HPRSecondOrderSolver::computeYSLICFlux(HPRReducedStateVector topTopStateVector, HPRReducedStateVector topStateVector, HPRReducedStateVector bottomStateVector,
                                                      HPRReducedStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                      HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRReducedStateVector evolvedBottomStateVector = HPRSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias,
                                                                                            slopeLimiter, 1, material1Parameters, material2Parameters);
    HPRReducedStateVector evolvedTopStateVector = HPRSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias,
                                                                                         slopeLimiter, 0, material1Parameters, material2Parameters);

    return HPRFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

void HPRSecondOrderSolver::computeSLICTimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                               int slopeLimiter, HPRMaterialParameters materialParameters)
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

void HPRSecondOrderSolver::computeSLICTimeStep(vector<HPRMultiphysicsStateVector> & currentCells, vector<HPRMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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

void HPRSecondOrderSolver::computeSLICTimeStep(vector<HPRIntermediateStateVector> & currentCells, vector<HPRIntermediateStateVector> & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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

void HPRSecondOrderSolver::computeSLICTimeStep(vector<HPRReducedStateVector> & currentCells, vector<HPRReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                               double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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

void HPRSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
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

void HPRSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<HPRIntermediateStateVector> > & currentCells, vector<vector<HPRIntermediateStateVector> > & currentCellsWithBoundary,
                                                  double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters,
                                                  HPRMaterialParameters material2Parameters)
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

void HPRSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, vector<vector<HPRReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                                  double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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

void HPRSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
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

void HPRSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<HPRIntermediateStateVector> > & currentCells, vector<vector<HPRIntermediateStateVector> > & currentCellsWithBoundary,
                                                  double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters,
                                                  HPRMaterialParameters material2Parameters)
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

void HPRSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, vector<vector<HPRReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                                  double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                        materialParameters);
        }

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 2);
        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                        materialParameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<HPRMultiphysicsStateVector> HPRSecondOrderSolver::solve(vector<HPRMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                               int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HPRMaterialParameters material1Parameters,
                                                               HPRMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<HPRMultiphysicsStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<HPRMultiphysicsStateVector> currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = HPRSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                            material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            // Runge-Kutta goes here.
        }

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 2);
        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            // Runge-Kutta goes here.
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

vector<HPRIntermediateStateVector> HPRSecondOrderSolver::solve(vector<HPRIntermediateStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                               int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HPRMaterialParameters material1Parameters,
                                                               HPRMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<HPRIntermediateStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<HPRIntermediateStateVector> currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = HPRSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                            material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter, material1Parameters,
                                                        material2Parameters);
        }

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 2);
        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter, material1Parameters,
                                                        material2Parameters);
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

vector<HPRReducedStateVector> HPRSecondOrderSolver::solve(vector<HPRReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                          int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HPRMaterialParameters material1Parameters,
                                                          HPRMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<HPRReducedStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<HPRReducedStateVector> currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = HPRSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                            material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter, material1Parameters,
                                                        material2Parameters);
        }

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 2);
        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter, material1Parameters,
                                                        material2Parameters);
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

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          materialParameters);
        }

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, materialParameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          materialParameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<HPRIntermediateStateVector> > HPRSecondOrderSolver::solve2D(vector<vector<HPRIntermediateStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                          double finalTime, double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                          HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<HPRIntermediateStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<HPRIntermediateStateVector> > currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = HPRSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                              material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);

            // Runge-Kutta goes here.
        }

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);
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

vector<vector<HPRReducedStateVector> > HPRSecondOrderSolver::solve2D(vector<vector<HPRReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                     double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                     HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<HPRReducedStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<HPRReducedStateVector> > currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = HPRSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                              material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          material1Parameters, material2Parameters);
        }

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentCellsWithBoundary = HPRSolvers::insertBoundaryCells2D(currentCells, 1);

            HPRForcingSolver::computeRungeKuttaTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
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
