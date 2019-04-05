#include "firstordersolver.h"

FirstOrderSolver::FirstOrderSolver()
{
}

vector<double> FirstOrderSolver::computeXLaxFriedrichsFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                          EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);

    return computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> FirstOrderSolver::computeXLaxFriedrichsFlux(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                           EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);

    return computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> FirstOrderSolver::computeXLaxFriedrichsFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                           HyperelasticMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);

    return computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> FirstOrderSolver::computeXLaxFriedrichsFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing,
                                                           double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);

    return computeLaxFriedrichsFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> FirstOrderSolver::computeYLaxFriedrichsFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                           EulerMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);

    return computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> FirstOrderSolver::computeYLaxFriedrichsFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                           EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(material1Parameters, material2Parameters);

    return computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> FirstOrderSolver::computeYLaxFriedrichsFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                           HyperelasticMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);

    return computeLaxFriedrichsFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
}

vector<double> FirstOrderSolver::computeLaxFriedrichsFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector, vector<double> leftFluxVector,
                                                          vector<double> rightFluxVector, double cellSpacing, double timeStep)
{
    return VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(VectorAlgebra::addVectors(leftFluxVector, rightFluxVector), VectorAlgebra::multiplyVector(
                                                                            (cellSpacing / timeStep), VectorAlgebra::subtractVectors(leftConservedVariableVector, rightConservedVariableVector))));
}

vector<double> FirstOrderSolver::computeXRichtmyerFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                      EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);

    vector<double> intermediateStateVector = computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
    return EulerStateVector::computeXFluxVector(intermediateStateVector, materialParameters);
}

vector<double> FirstOrderSolver::computeXRichtmyerFlux(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                       EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);

    vector<double> intermediateStateVector = computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
    return EulerMultiphysicsStateVector::computeXFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> FirstOrderSolver::computeXRichtmyerFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                       HyperelasticMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(materialParameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(materialParameters);
    vector<double> intermediateStateVector = computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
    return ElasticStateVector::computeXFluxVector(intermediateStateVector, materialParameters);
}

vector<double> FirstOrderSolver::computeXRichtmyerFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                       HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> leftFluxVector = leftStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> rightFluxVector = rightStateVector.computeXFluxVector(material1Parameters, material2Parameters);
    vector<double> intermediateStateVector = computeRichtmyerFlux(leftConservedVariableVector, rightConservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
    return ElasticMultiphysicsStateVector::computeXFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> FirstOrderSolver::computeYRichtmyerFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                       EulerMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);

    vector<double> intermediateStateVector = computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
    return EulerStateVector::computeYFluxVector(intermediateStateVector, materialParameters);
}

vector<double> FirstOrderSolver::computeYRichtmyerFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                       EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(material1Parameters, material2Parameters);

    vector<double> intermediateStateVector = computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
    return EulerMultiphysicsStateVector::computeYFluxVector(intermediateStateVector, material1Parameters, material2Parameters);
}

vector<double> FirstOrderSolver::computeYRichtmyerFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                       HyperelasticMaterialParameters materialParameters)
{
    vector<double> topConservedVariableVector = topStateVector.computeConservedVariableVector(materialParameters);
    vector<double> bottomConservedVariableVector = bottomStateVector.computeConservedVariableVector(materialParameters);

    vector<double> topFluxVector = topStateVector.computeYFluxVector(materialParameters);
    vector<double> bottomFluxVector = bottomStateVector.computeYFluxVector(materialParameters);

    vector<double> intermediateStateVector = computeRichtmyerFlux(topConservedVariableVector, bottomConservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep);
    return ElasticStateVector::computeYFluxVector(intermediateStateVector, materialParameters);
}

vector<double> FirstOrderSolver::computeRichtmyerFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector, vector<double> leftFluxVector,
                                                      vector<double> rightFluxVector, double cellSpacing, double timeStep)
{
    return VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(VectorAlgebra::addVectors(leftConservedVariableVector, rightConservedVariableVector),
                                                                        VectorAlgebra::multiplyVector((timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector))));
}

vector<double> FirstOrderSolver::computeXFORCEFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                   EulerMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);

    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> FirstOrderSolver::computeXFORCEFlux(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                   EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> FirstOrderSolver::computeXFORCEFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                   HyperelasticMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, materialParameters);

    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> FirstOrderSolver::computeXFORCEFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                   HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeXLaxFriedrichsFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeXRichtmyerFlux(leftStateVector, rightStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> FirstOrderSolver::computeYFORCEFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                   EulerMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);

    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> FirstOrderSolver::computeYFORCEFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                   EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, material1Parameters, material2Parameters);

    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> FirstOrderSolver::computeYFORCEFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                   HyperelasticMaterialParameters materialParameters)
{
    vector<double> laxFriedrichsFlux = computeYLaxFriedrichsFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);
    vector<double> richtmyerFlux = computeYRichtmyerFlux(topStateVector, bottomStateVector, cellSpacing, timeStep, materialParameters);

    return computeFORCEFlux(laxFriedrichsFlux, richtmyerFlux);
}

vector<double> FirstOrderSolver::computeFORCEFlux(vector<double> laxFriedrichsFlux, vector<double> richtmyerFlux)
{
    return VectorAlgebra::multiplyVector(0.5, VectorAlgebra::addVectors(laxFriedrichsFlux, richtmyerFlux));
}

void FirstOrderSolver::computeFORCETimeStep(vector<EulerStateVector> & currentCells, vector<EulerStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                            EulerMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, materialParameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, materialParameters);

        currentCells[i].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
    }
}

void FirstOrderSolver::computeFORCETimeStep(vector<EulerMultiphysicsStateVector> & currentCells, vector<EulerMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                            double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), material1Parameters, material2Parameters);
    }
}

void FirstOrderSolver::computeFORCETimeStep(vector<ElasticStateVector> & currentCells, vector<ElasticStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(materialParameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, materialParameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, materialParameters);

        currentCells[i].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
    }
}

void FirstOrderSolver::computeFORCETimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, vector<ElasticMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                            double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = currentCells[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> leftFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i], currentCellsWithBoundary[i + 1], cellSpacing, timeStep, material1Parameters, material2Parameters);
        vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1], currentCellsWithBoundary[i + 2], cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), material1Parameters, material2Parameters);
    }
}

void FirstOrderSolver::computeXFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                              EulerMaterialParameters materialParameters)
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

            currentCells[i][j].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
        }
    }
}

void FirstOrderSolver::computeXFORCETimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                               double cellSpacing, double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
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
            vector<double> rightFluxVector = computeXFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 1][j + 2], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);

            currentCells[i][j].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), material1Parameters,
                                                          material2Parameters);
        }
    }
}

void FirstOrderSolver::computeXFORCETimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, HyperelasticMaterialParameters materialParameters)
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

            currentCells[i][j].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, leftFluxVector, rightFluxVector, cellSpacing, timeStep), materialParameters);
        }
    }
}

void FirstOrderSolver::computeYFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                               EulerMaterialParameters materialParameters)
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

            currentCells[i][j].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep), materialParameters);
        }
    }
}

void FirstOrderSolver::computeYFORCETimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                               double cellSpacing, double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> conservedVariableVector = currentCells[i][j].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> topFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i][j + 1], currentCellsWithBoundary[i + 1][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);
            vector<double> bottomFluxVector = computeYFORCEFlux(currentCellsWithBoundary[i + 1][j + 1], currentCellsWithBoundary[i + 2][j + 1], cellSpacing, timeStep, material1Parameters,
                    material2Parameters);

            currentCells[i][j].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep), material1Parameters,
                                                          material2Parameters);
        }
    }
}

void FirstOrderSolver::computeYFORCETimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                               double timeStep, HyperelasticMaterialParameters materialParameters)
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

            currentCells[i][j].setConservedVariableVector(computeFORCEUpdate(conservedVariableVector, topFluxVector, bottomFluxVector, cellSpacing, timeStep), materialParameters);
        }
    }
}

vector<double> FirstOrderSolver::computeFORCEUpdate(vector<double> conservedVariableVector, vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep)
{
    return VectorAlgebra::addVectors(conservedVariableVector, VectorAlgebra::multiplyVector((timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector)));
}

vector<EulerStateVector> FirstOrderSolver::solve(vector<EulerStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                 EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<EulerStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<EulerStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

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

vector<EulerMultiphysicsStateVector> FirstOrderSolver::solve(vector<EulerMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                             int subcyclingIterations, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<EulerMultiphysicsStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<EulerMultiphysicsStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters, material2Parameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

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

vector<ElasticStateVector> FirstOrderSolver::solve(vector<ElasticStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                   HyperelasticMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<ElasticMultiphysicsStateVector> FirstOrderSolver::solve(vector<ElasticMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                               int subcyclingIterations, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<ElasticMultiphysicsStateVector> currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<ElasticMultiphysicsStateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                         material2Parameters);

        computeFORCETimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}

vector<vector<EulerStateVector> > FirstOrderSolver::solve2D(vector<vector<EulerStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                            EulerMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<EulerStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<EulerStateVector> > currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

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

vector<vector<EulerMultiphysicsStateVector> > FirstOrderSolver::solve2D(vector<vector<EulerMultiphysicsStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                        int subcyclingIterations, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<EulerMultiphysicsStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<EulerMultiphysicsStateVector> > currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters,
                                                           material2Parameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);
        currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, material1Parameters, material2Parameters);

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

vector<vector<ElasticStateVector> > FirstOrderSolver::solve2D(vector<vector<ElasticStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                              int subcyclingIterations, HyperelasticMaterialParameters materialParameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<vector<ElasticStateVector> > currentCells = initialCells;

    while (currentTime < finalTime)
    {
        vector<vector<ElasticStateVector> > currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep2D(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, materialParameters);

        computeXFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);
        currentCellsWithBoundary = Solvers::insertBoundaryCells2D(currentCells, 1);
        computeYFORCETimeStep2D(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, materialParameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return currentCells;
}
