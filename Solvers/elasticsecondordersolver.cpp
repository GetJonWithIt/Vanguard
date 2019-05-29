#include "elasticsecondordersolver.h"

ElasticSecondOrderSolver::ElasticSecondOrderSolver()
{
}

vector<double> ElasticSecondOrderSolver::computeXSLICFlux(ElasticStateVector leftLeftStateVector, ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                          ElasticStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                          HyperelasticMaterialParameters materialParameters)
{
    ElasticStateVector evolvedRightStateVector = ElasticSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                            slopeLimiter, 1, materialParameters);
    ElasticStateVector evolvedLeftStateVector = ElasticSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias,
                                                                                           slopeLimiter, 0, materialParameters);

    return ElasticFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> ElasticSecondOrderSolver::computeXSLICFlux(ElasticMultiphysicsStateVector leftLeftStateVector, ElasticMultiphysicsStateVector leftStateVector,
                                                          ElasticMultiphysicsStateVector rightStateVector, ElasticMultiphysicsStateVector rightRightStateVector, double cellSpacing,
                                                          double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                                          HyperelasticMaterialParameters material2Parameters)
{
    ElasticMultiphysicsStateVector evolvedRightStateVector = ElasticSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                        bias, slopeLimiter, 1, material1Parameters, material2Parameters);
    ElasticMultiphysicsStateVector evolvedLeftStateVector = ElasticSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep,
                                                                                                       bias, slopeLimiter, 0, material1Parameters, material2Parameters);

    return ElasticFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> ElasticSecondOrderSolver::computeXSLICFlux(ElasticReducedStateVector leftLeftStateVector, ElasticReducedStateVector leftStateVector,
                                                          ElasticReducedStateVector rightStateVector, ElasticReducedStateVector rightRightStateVector, double cellSpacing,
                                                          double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                                          HyperelasticMaterialParameters material2Parameters)
{
    ElasticReducedStateVector evolvedRightStateVector = ElasticSolvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                   bias, slopeLimiter, 1, material1Parameters, material2Parameters);
    ElasticReducedStateVector evolvedLeftStateVector = ElasticSolvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep,
                                                                                                  bias, slopeLimiter, 0, material1Parameters, material2Parameters);

    return ElasticFirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> ElasticSecondOrderSolver::computeYSLICFlux(ElasticStateVector topTopStateVector, ElasticStateVector topStateVector, ElasticStateVector bottomStateVector,
                                                          ElasticStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                          HyperelasticMaterialParameters materialParameters)
{
    ElasticStateVector evolvedBottomStateVector = ElasticSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias,
                                                                                             slopeLimiter, 1, materialParameters);
    ElasticStateVector evolvedTopStateVector = ElasticSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias,
                                                                                          slopeLimiter, 0, materialParameters);

    return ElasticFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> ElasticSecondOrderSolver::computeYSLICFlux(ElasticReducedStateVector topTopStateVector, ElasticReducedStateVector topStateVector,
                                                          ElasticReducedStateVector bottomStateVector, ElasticReducedStateVector bottomBottomStateVector, double cellSpacing,
                                                          double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                                          HyperelasticMaterialParameters material2Parameters)
{
    ElasticReducedStateVector evolvedBottomStateVector = ElasticSolvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep,
                                                                                                    bias, slopeLimiter, 1, material1Parameters, material2Parameters);
    ElasticReducedStateVector evolvedTopStateVector = ElasticSolvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep,
                                                                                                 bias, slopeLimiter, 0, material1Parameters, material2Parameters);

    return ElasticFirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

void ElasticSecondOrderSolver::computeSLICTimeStep(vector<ElasticStateVector> & currentCells, vector<ElasticStateVector> & currentCellsWithBoundary, double cellSpacing,
                                                   double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);

        vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                cellSpacing, timeStep, bias, slopeLimiter, materialParameters);
        vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                currentCellsWithBoundary[i + 4], cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                   materialParameters);
    }
}

void ElasticSecondOrderSolver::computeSLICTimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, vector<ElasticMultiphysicsStateVector> & currentCellsWithBoundary,
                                                   double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                                   HyperelasticMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);

        vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                currentCellsWithBoundary[i + 4], cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                   material1Parameters, material2Parameters);
    }
}

void ElasticSecondOrderSolver::computeSLICTimeStep(vector<ElasticReducedStateVector> & currentCells, vector<ElasticReducedStateVector> & currentCellsWithBoundary, double cellSpacing,
                                                   double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                                   HyperelasticMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);

        vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                currentCellsWithBoundary[i + 4], cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                   material1Parameters, material2Parameters);
    }
}

void ElasticSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary,
                                                      double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters)
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

void ElasticSecondOrderSolver::computeXSLICTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                                      double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                                      HyperelasticMaterialParameters material2Parameters)
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

void ElasticSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary,
                                                      double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters)
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

void ElasticSecondOrderSolver::computeYSLICTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                                      double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                                      HyperelasticMaterialParameters material2Parameters)
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

vector<ElasticStateVector> ElasticSecondOrderSolver::solve(vector<ElasticStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                           int slopeLimiter, int subcyclingIterations, HyperelasticMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticStateVector> currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = ElasticSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<ElasticMultiphysicsStateVector> ElasticSecondOrderSolver::solve(vector<ElasticMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient,
                                                                       double finalTime, double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                       HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticMultiphysicsStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticMultiphysicsStateVector> currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = ElasticSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                         material2Parameters);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

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

vector<ElasticReducedStateVector> ElasticSecondOrderSolver::solve(vector<ElasticReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                  double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                  HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticReducedStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticReducedStateVector> currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentCells, 2);
        double timeStep = ElasticSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                         material2Parameters);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        if (reinitialisationFrequency != 0 && currentIteration != 0)
        {
            if ((currentIteration % reinitialisationFrequency) == 0)
            {
                MultiphysicsSolvers::reinitialiseVolumeFraction(currentCells);
            }
        }

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<ElasticStateVector> > ElasticSecondOrderSolver::solve2D(vector<vector<ElasticStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                      double finalTime, double bias, int slopeLimiter, int subcyclingIterations,
                                                                      HyperelasticMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<ElasticStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<ElasticStateVector> > currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = ElasticSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<ElasticReducedStateVector> > ElasticSecondOrderSolver::solve2D(vector<vector<ElasticReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                             double finalTime, double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                             HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<ElasticReducedStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<ElasticReducedStateVector> > currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = ElasticSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                                  material2Parameters);

        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 2);
        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        if (reinitialisationFrequency != 0 && currentIteration != 0)
        {
            if ((currentIteration % reinitialisationFrequency) == 0)
            {
                MultiphysicsSolvers::reinitialiseVolumeFraction(currentCells);
            }
        }

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}
