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

vector<double> MHDFirstOrderSolver::computeXLaxFriedrichsFlux(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                              MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> MHDFirstOrderSolver::computeXLaxFriedrichsFlux(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                                              MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> MHDFirstOrderSolver::computeXLaxFriedrichsFlux(MHDReducedStateVector leftStateVector, MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                              MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> MHDFirstOrderSolver::computeYLaxFriedrichsFlux(MHDStateVector topStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                              MHDMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> MHDFirstOrderSolver::computeYLaxFriedrichsFlux(MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                              MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(material1Parameters, material2Parameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> MHDFirstOrderSolver::computeYLaxFriedrichsFlux(MHDReducedStateVector topStateVector, MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                              MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(material1Parameters, material2Parameters);

    return FirstOrderSolver::computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
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

vector<double> MHDFirstOrderSolver::computeXRichtmyerFlux(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                          MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing,
                                                                                    timeStep);

    return MHDMultiphysicsStateVector::computeXFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> MHDFirstOrderSolver::computeXRichtmyerFlux(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                                          MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing,
                                                                                    timeStep);

    return MHDIntermediateStateVector::computeXFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> MHDFirstOrderSolver::computeXRichtmyerFlux(MHDReducedStateVector leftStateVector, MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                          MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing,
                                                                                    timeStep);

    return MHDReducedStateVector::computeXFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> MHDFirstOrderSolver::computeYRichtmyerFlux(MHDStateVector topStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                          MHDMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing,
                                                                                    timeStep);

    return MHDStateVector::computeYFluxVector(intermediateStateVector, materialParameters);
}

vector<double> MHDFirstOrderSolver::computeYRichtmyerFlux(MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                          MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing,
                                                                                    timeStep);

    return MHDIntermediateStateVector::computeYFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> MHDFirstOrderSolver::computeYRichtmyerFlux(MHDReducedStateVector topStateVector, MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                          MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> intermediateStateVector = FirstOrderSolver::computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing,
                                                                                    timeStep);

    return MHDReducedStateVector::computeYFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> MHDFirstOrderSolver::computeXFORCEFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep,
                                                      MHDMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> MHDFirstOrderSolver::computeXFORCEFlux(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                      MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> MHDFirstOrderSolver::computeXFORCEFlux(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                                      MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> MHDFirstOrderSolver::computeXFORCEFlux(MHDReducedStateVector leftStateVector, MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                      MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> MHDFirstOrderSolver::computeYFORCEFlux(MHDStateVector topStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                      MHDMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> MHDFirstOrderSolver::computeYFORCEFlux(MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                      MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return FirstOrderSolver::computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> MHDFirstOrderSolver::computeYFORCEFlux(MHDReducedStateVector topStateVector, MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                      MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

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

void MHDFirstOrderSolver::computeFORCETimeStep(vector<MHDMultiphysicsStateVector> & currentCells, vector<MHDMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), material1Parameters,
                                                   material2Parameters);
    }
}

void MHDFirstOrderSolver::computeFORCETimeStep(vector<MHDIntermediateStateVector> & currentCells, vector<MHDIntermediateStateVector> & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), material1Parameters,
                                                   material2Parameters);
    }
}

void MHDFirstOrderSolver::computeFORCETimeStep(vector<MHDReducedStateVector> & currentCells, vector<MHDReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                               MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), material1Parameters,
                                                   material2Parameters);
    }
}

void MHDFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  MHDMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
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

void MHDFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<MHDIntermediateStateVector> > & currentCells, vector<vector<MHDIntermediateStateVector> > & currentCellsWithBoundary,
                                                  double cellSpacing, double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);
            vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                          material1Parameters, material2Parameters);
        }
    }
}

void MHDFirstOrderSolver::computeXFORCETimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                                  double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);
            vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep),
                                                          material1Parameters, material2Parameters);
        }
    }
}

void MHDFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                  MHDMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
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

void MHDFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<MHDIntermediateStateVector> > & currentCells, vector<vector<MHDIntermediateStateVector> > & currentCellsWithBoundary,
                                                  double cellSpacing, double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> topFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);
            vector<double> bottomFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep),
                                                          material1Parameters, material2Parameters);
        }
    }
}

void MHDFirstOrderSolver::computeYFORCETimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                                  double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> topFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);
            vector<double> bottomFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);

            currentCells[i][j].setConservedVariableVector(FirstOrderSolver::computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep),
                                                          material1Parameters, material2Parameters);
        }
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

vector<MHDMultiphysicsStateVector> MHDFirstOrderSolver::solve(vector<MHDMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                              int subcyclingIterations, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<MHDMultiphysicsStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<MHDMultiphysicsStateVector> currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = MHDSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                            material2Parameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<MHDIntermediateStateVector> MHDFirstOrderSolver::solve(vector<MHDIntermediateStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                              int subcyclingIterations, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<MHDIntermediateStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<MHDIntermediateStateVector> currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = MHDSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                            material2Parameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<MHDReducedStateVector> MHDFirstOrderSolver::solve(vector<MHDReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                         MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<MHDReducedStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<MHDReducedStateVector> currentCellsWithBoundary = MHDSolvers::insertBoundaryCells(currentCells, 1);
        double timeStep = MHDSolvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                            material2Parameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<MHDStateVector> > MHDFirstOrderSolver::solve2D(vector<vector<MHDStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                             MHDMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<MHDStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<MHDStateVector> > currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = MHDSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, materialParameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<MHDIntermediateStateVector> > MHDFirstOrderSolver::solve2D(vector<vector<MHDIntermediateStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                         int subcyclingIterations, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<MHDIntermediateStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<MHDIntermediateStateVector> > currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = MHDSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                              material2Parameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, material1Parameters, material2Parameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<MHDReducedStateVector> > MHDFirstOrderSolver::solve2D(vector<vector<MHDReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                    int subcyclingIterations, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<MHDReducedStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<MHDReducedStateVector> > currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = MHDSolvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                              material2Parameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, material1Parameters, material2Parameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentCellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentCells, 1);
        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, 0.5 * timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}
