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

MHDMultimaterialSystem::MHDMultimaterialSystem(vector<vector<MHDStateVector> > newMaterial1Cells2D, vector<vector<MHDStateVector> > newMaterial2Cells2D,
                                               vector<vector<double> > newLevelSetFunction2D)
{
    material1Cells2D = newMaterial1Cells2D;
    material2Cells2D = newMaterial2Cells2D;
    levelSetFunction2D = newLevelSetFunction2D;
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

void MHDMultimaterialSystem::setMaterial1Cells2D(vector<vector<MHDStateVector> > newMaterial1Cells2D)
{
    material1Cells2D = newMaterial1Cells2D;
}

void MHDMultimaterialSystem::setMaterial2Cells2D(vector<vector<MHDStateVector> > newMaterial2Cells2D)
{
    material2Cells2D = newMaterial2Cells2D;
}

void MHDMultimaterialSystem::setLevelSetFunction2D(vector<vector<double> > newLevelSetFunction2D)
{
    levelSetFunction2D = newLevelSetFunction2D;
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

vector<vector<MHDStateVector> > MHDMultimaterialSystem::getMaterial1Cells2D()
{
    return material1Cells2D;
}

vector<vector<MHDStateVector> > MHDMultimaterialSystem::getMaterial2Cells2D()
{
    return material2Cells2D;
}

vector<vector<double> > MHDMultimaterialSystem::getLevelSetFunction2D()
{
    return levelSetFunction2D;
}
