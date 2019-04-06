#include "secondordersolver.h"

SecondOrderSolver::SecondOrderSolver()
{
}

vector<double> SecondOrderSolver::computeXSLICFlux(EulerStateVector leftLeftStateVector, EulerStateVector leftStateVector, EulerStateVector rightStateVector,
                                                   EulerStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                   EulerMaterialParameters materialParameters)
{
    EulerStateVector evolvedRightStateVector = Solvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                   materialParameters);
    EulerStateVector evolvedLeftStateVector = Solvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                  materialParameters);

    return FirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> SecondOrderSolver::computeXSLICFlux(EulerMultiphysicsStateVector leftLeftStateVector, EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector,
                                                   EulerMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                   EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    EulerMultiphysicsStateVector evolvedRightStateVector = Solvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                               1, material1Parameters, material2Parameters);
    EulerMultiphysicsStateVector evolvedLeftStateVector = Solvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                              0, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> SecondOrderSolver::computeXSLICFlux(ElasticStateVector leftLeftStateVector, ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                   ElasticStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                   HyperelasticMaterialParameters materialParameters)
{
    ElasticStateVector evolvedRightStateVector = Solvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                     materialParameters);
    ElasticStateVector evolvedLeftStateVector = Solvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                    materialParameters);

    return FirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> SecondOrderSolver::computeXSLICFlux(ElasticMultiphysicsStateVector leftLeftStateVector, ElasticMultiphysicsStateVector leftStateVector,
                                                   ElasticMultiphysicsStateVector rightStateVector, ElasticMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep,
                                                   double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    ElasticMultiphysicsStateVector evolvedRightStateVector = Solvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                 1, material1Parameters, material2Parameters);
    ElasticMultiphysicsStateVector evolvedLeftStateVector = Solvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                0, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> SecondOrderSolver::computeXSLICFlux(ElasticReducedStateVector leftLeftStateVector, ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector,
                                                   ElasticReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                   HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    ElasticReducedStateVector evolvedRightStateVector = Solvers::evolveStateByHalfXTimeStep(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                            1, material1Parameters, material2Parameters);
    ElasticReducedStateVector evolvedLeftStateVector = Solvers::evolveStateByHalfXTimeStep(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                           0, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeXFORCEFlux(evolvedRightStateVector, evolvedLeftStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> SecondOrderSolver::computeYSLICFlux(EulerStateVector topTopStateVector, EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerStateVector bottomBottomStateVector,
                                                   double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
{
    EulerStateVector evolvedBottomStateVector = Solvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                    materialParameters);
    EulerStateVector evolvedTopStateVector = Solvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                 materialParameters);

    return FirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> SecondOrderSolver::computeYSLICFlux(EulerMultiphysicsStateVector topTopStateVector, EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector,
                                                   EulerMultiphysicsStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                   EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    EulerMultiphysicsStateVector evolvedBottomStateVector = Solvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                1, material1Parameters, material2Parameters);
    EulerMultiphysicsStateVector evolvedTopStateVector = Solvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                             0, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

vector<double> SecondOrderSolver::computeYSLICFlux(ElasticStateVector topTopStateVector, ElasticStateVector topStateVector, ElasticStateVector bottomStateVector,
                                                   ElasticStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                   HyperelasticMaterialParameters materialParameters)
{
    ElasticStateVector evolvedBottomStateVector = Solvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter, 1,
                                                                                      materialParameters);
    ElasticStateVector evolvedTopStateVector = Solvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias, slopeLimiter, 0,
                                                                                   materialParameters);

    return FirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, materialParameters);
}

vector<double> SecondOrderSolver::computeYSLICFlux(ElasticReducedStateVector topTopStateVector, ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector,
                                                   ElasticReducedStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                   HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    ElasticReducedStateVector evolvedBottomStateVector = Solvers::evolveStateByHalfYTimeStep(topTopStateVector, topStateVector, bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                             1, material1Parameters, material2Parameters);
    ElasticReducedStateVector evolvedTopStateVector = Solvers::evolveStateByHalfYTimeStep(topStateVector, bottomStateVector, bottomBottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                          0, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeYFORCEFlux(evolvedBottomStateVector, evolvedTopStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
}

void SecondOrderSolver::computeSLICTimeStep(vector<EulerStateVector> & currentCells, vector<EulerStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                            int slopeLimiter, EulerMaterialParameters materialParameters)
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

void SecondOrderSolver::computeSLICTimeStep(vector<EulerMultiphysicsStateVector> & currentCells, vector<EulerMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                            double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);

        vector<double> leftFluxVector = computeXSLICFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3],
                cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXSLICFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], currentCellsWithBoundary[i + 3], currentCellsWithBoundary[i + 4],
                cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                    material1Parameters, material2Parameters);
    }
}

void SecondOrderSolver::computeSLICTimeStep(vector<ElasticStateVector> & currentCells, vector<ElasticStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                            double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters)
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

void SecondOrderSolver::computeSLICTimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, vector<ElasticMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                            double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

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

void SecondOrderSolver::computeSLICTimeStep(vector<ElasticReducedStateVector> & currentCells, vector<ElasticReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                            double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

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

void SecondOrderSolver::computeXSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
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

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
        }
    }
}

void SecondOrderSolver::computeXSLICTimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                               double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters,
                                               EulerMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

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

void SecondOrderSolver::computeXSLICTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters)
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

void SecondOrderSolver::computeXSLICTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                               double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                               HyperelasticMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

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

void SecondOrderSolver::computeYSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters)
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

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep), materialParameters);
        }
    }
}

void SecondOrderSolver::computeYSLICTimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                               double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters,
                                               EulerMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

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

void SecondOrderSolver::computeYSLICTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters)
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

void SecondOrderSolver::computeYSLICTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                               double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                               HyperelasticMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

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

vector<EulerStateVector> SecondOrderSolver::solve(vector<EulerStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                  int subcyclingIterations, EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<EulerStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<EulerStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            // Runge-Kutta goes here.
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<EulerMultiphysicsStateVector> SecondOrderSolver::solve(vector<EulerMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                              int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters,
                                                              EulerMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<EulerMultiphysicsStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<EulerMultiphysicsStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters, material2Parameters);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
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

vector<ElasticStateVector> SecondOrderSolver::solve(vector<ElasticStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                    int subcyclingIterations, HyperelasticMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeSLICTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<ElasticMultiphysicsStateVector> SecondOrderSolver::solve(vector<ElasticMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                                int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                                HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticMultiphysicsStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticMultiphysicsStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
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

vector<ElasticReducedStateVector> SecondOrderSolver::solve(vector<ElasticReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                           int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                           HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticReducedStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticReducedStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
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

vector<vector<EulerStateVector> > SecondOrderSolver::solve2D(vector<vector<EulerStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                             int slopeLimiter, int subcyclingIterations, EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<EulerStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<EulerStateVector> > currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            // Runge-Kutta goes here.
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<EulerMultiphysicsStateVector> > SecondOrderSolver::solve2D(vector<vector<EulerMultiphysicsStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                         double finalTime, double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                         EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<EulerMultiphysicsStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<EulerMultiphysicsStateVector> > currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                           material2Parameters);

        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
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

vector<vector<ElasticStateVector> > SecondOrderSolver::solve2D(vector<vector<ElasticStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                               double bias, int slopeLimiter, int subcyclingIterations, HyperelasticMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<ElasticStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<ElasticStateVector> > currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<ElasticReducedStateVector> > SecondOrderSolver::solve2D(vector<vector<ElasticReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                      double finalTime, double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                      HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<ElasticReducedStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<ElasticReducedStateVector> > currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                           material2Parameters);

        computeXSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

        currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 2);
        computeYSLICTimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

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
