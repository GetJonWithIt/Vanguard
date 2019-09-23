#include "mhdmultimaterialsystem.h"

MHDMultimaterialSystem::MHDMultimaterialSystem()
{
}

MHDMultimaterialSystem::MHDMultimaterialSystem(vector<MHDStateVector> newMaterial1Cells, vector<MHDStateVector> newMaterial2Cells, vector<double> newLevelSetFunction)
{
    material1Cells = newMaterial1Cells;
    material2Cells = newMaterial2Cells;
    levelSetFunction = newLevelSetFunction;
}

void MHDMultimaterialSystem::setMaterial1Cells(vector<MHDStateVector> newMaterial1Cells)
{
    material1Cells = newMaterial1Cells;
}

void MHDMultimaterialSystem::setMaterial2Cells(vector<MHDStateVector> newMaterial2Cells)
{
    material2Cells = newMaterial2Cells;
}

void MHDMultimaterialSystem::setLevelSetFunction(vector<double> newLevelSetFunction)
{
    levelSetFunction = newLevelSetFunction;
}

vector<MHDStateVector> MHDMultimaterialSystem::getMaterial1Cells()
{
    return material1Cells;
}

vector<MHDStateVector> MHDMultimaterialSystem::getMaterial2Cells()
{
    return material2Cells;
}

vector<double> MHDMultimaterialSystem::getLevelSetFunction()
{
    return levelSetFunction;
}
