#include "elasticfirstordersolver.h"

ElasticFirstOrderSolver::ElasticFirstOrderSolver()
{
}

vector<double> ElasticFirstOrderSolver::computeXLaxFriedrichsFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                                  HyperelasticMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> ElasticFirstOrderSolver::computeXLaxFriedrichsFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector,
                                                                  double cellSpacing, double timeStep, HyperelasticMaterialParameters material1Parameters,
                                                                  HyperelasticMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> ElasticFirstOrderSolver::computeXLaxFriedrichsFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing,
                                                                  double timeStep, HyperelasticMaterialParameters material1Parameters,
                                                                  HyperelasticMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> ElasticFirstOrderSolver::computeYLaxFriedrichsFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                                  HyperelasticMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> ElasticFirstOrderSolver::computeYLaxFriedrichsFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing,
                                                                  double timeStep, HyperelasticMaterialParameters material1Parameters,
                                                                  HyperelasticMaterialParameters material2Parameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(material1Parameters, material2Parameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> ElasticFirstOrderSolver::computeXRichtmyerFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                              HyperelasticMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                    cellSpacing, timeStep);
    return ElasticStateVector::computeXFluxVector(intermediateStateVector, materialParameters);
}

vector<double> ElasticFirstOrderSolver::computeXRichtmyerFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing,
                                                              double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                    cellSpacing, timeStep);
    return ElasticMultiphysicsStateVector::computeXFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> ElasticFirstOrderSolver::computeXRichtmyerFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing,
                                                              double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector,
                                                                                    cellSpacing, timeStep);
    return ElasticReducedStateVector::computeXFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> ElasticFirstOrderSolver::computeYRichtmyerFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                              HyperelasticMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);

    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                    cellSpacing, timeStep);
    return ElasticStateVector::computeYFluxVector(intermediateStateVector, materialParameters);
}

vector<double> ElasticFirstOrderSolver::computeYRichtmyerFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing,
                                                              double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(material1Parameters, material2Parameters);

    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector,
                                                                                    cellSpacing, timeStep);
    return ElasticReducedStateVector::computeYFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> ElasticFirstOrderSolver::computeXFORCEFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                          HyperelasticMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> ElasticFirstOrderSolver::computeXFORCEFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing,
                                                          double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> ElasticFirstOrderSolver::computeXFORCEFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing,
                                                          double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> ElasticFirstOrderSolver::computeYFORCEFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                          HyperelasticMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> ElasticFirstOrderSolver::computeYFORCEFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing,
                                                          double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

void ElasticFirstOrderSolver::computeFORCETimeStep(vector<ElasticStateVector> & currentCells, vector<ElasticStateVector> & currentCellsWithBoundary, double cellSpacing,
                                                   double timeStep, HyperelasticMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, materialParameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, materialParameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                   materialParameters);
    }
}

void ElasticFirstOrderSolver::computeFORCETimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, vector<ElasticMultiphysicsStateVector> & currentCellsWithBoundary,
                                                   double cellSpacing, double timeStep, HyperelasticMaterialParameters material1Parameters,
                                                   HyperelasticMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, material1Parameters,
                material2Parameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, material1Parameters,
                material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                   material1Parameters, material2Parameters);
    }
}

void ElasticFirstOrderSolver::computeFORCETimeStep(vector<ElasticReducedStateVector> & currentCells, vector<ElasticReducedStateVector> & currentCellsWithBoundary,
                                                   double cellSpacing, double timeStep, HyperelasticMaterialParameters material1Parameters,
                                                   HyperelasticMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, material1Parameters,
                material2Parameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, material1Parameters,
                material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                   material1Parameters, material2Parameters);
    }
}

void ElasticFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary,
                                                      double cellSpacing, double timeStep, HyperelasticMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);

            vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep,
                    materialParameters);
            vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2], cellSpacing, timeStep,
                   materialParameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                          materialParameters);
        }
    }
}

void ElasticFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                               double cellSpacing, double timeStep, HyperelasticMaterialParameters material1Parameters,
                                                      HyperelasticMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);
            vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2], cellSpacing, timeStep,
                    material1Parameters, material2Parameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                          material1Parameters, material2Parameters);
        }
    }
}

void ElasticFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary,
                                                      double cellSpacing, double timeStep, HyperelasticMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(materialParameters);

            vector<double> topFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep,
                    materialParameters);
            vector<double> bottomFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1], cellSpacing, timeStep,
                    materialParameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep),
                                                          materialParameters);
        }
    }
}

void ElasticFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                                      double cellSpacing, double timeStep, HyperelasticMaterialParameters material1Parameters,
                                                      HyperelasticMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> topFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep,
                    material1Parameters, material2Parameters);
            vector<double> bottomFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1], cellSpacing, timeStep,
                    material1Parameters, material2Parameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep),
                                                          material1Parameters, material2Parameters);
        }
    }
}

vector<ElasticStateVector> ElasticFirstOrderSolver::solve(vector<ElasticStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                          int subcyclingIterations, HyperelasticMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticStateVector> currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = ElasticSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<ElasticMultiphysicsStateVector> ElasticFirstOrderSolver::solve(vector<ElasticMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient,
                                                                      double finalTime, int subcyclingIterations, HyperelasticMaterialParameters material1Parameters,
                                                                      HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticMultiphysicsStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticMultiphysicsStateVector> currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = ElasticSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                                material2Parameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<ElasticReducedStateVector> ElasticFirstOrderSolver::solve(vector<ElasticReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                 int subcyclingIterations, HyperelasticMaterialParameters material1Parameters,
                                                                 HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticReducedStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticReducedStateVector> currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = ElasticSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                                material2Parameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<ElasticStateVector> > ElasticFirstOrderSolver::solve2D(vector<vector<ElasticStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                     double finalTime, int subcyclingIterations, HyperelasticMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<ElasticStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<ElasticStateVector> > currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = ElasticSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);
        currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<ElasticReducedStateVector> > ElasticFirstOrderSolver::solve2D(vector<vector<ElasticReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                            double finalTime, int subcyclingIterations, HyperelasticMaterialParameters material1Parameters,
                                                                            HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<ElasticReducedStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<ElasticReducedStateVector> > currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = ElasticSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                           material2Parameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);
        currentCellsWithBoundary = ElasticSolvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}
