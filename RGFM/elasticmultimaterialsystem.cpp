#include "elasticmultimaterialsystem.h"

ElasticMultimaterialSystem::ElasticMultimaterialSystem()
{
}

ElasticMultimaterialSystem::ElasticMultimaterialSystem(vector<ElasticStateVector> newMaterial1Cells, vector<ElasticStateVector> newMaterial2Cells, vector<double> newLevelSetFunction)
{
    material1Cells = newMaterial1Cells;
    material2Cells = newMaterial2Cells;
    levelSetFunction = newLevelSetFunction;
}

void ElasticMultimaterialSystem::setMaterial1Cells(vector<ElasticStateVector> newMaterial1Cells)
{
    material1Cells = newMaterial1Cells;
}

void ElasticMultimaterialSystem::setMaterial2Cells(vector<ElasticStateVector> newMaterial2Cells)
{
    material2Cells = newMaterial2Cells;
}

void ElasticMultimaterialSystem::setLevelSetFunction(vector<double> newLevelSetFunction)
{
    levelSetFunction = newLevelSetFunction;
}

vector<ElasticStateVector> ElasticMultimaterialSystem::getMaterial1Cells()
{
    return material1Cells;
}

vector<ElasticStateVector> ElasticMultimaterialSystem::getMaterial2Cells()
{
    return material2Cells;
}

vector<double> ElasticMultimaterialSystem::getLevelSetFunction()
{
    return levelSetFunction;
}
