#include "multimaterialsystem.h"

MultimaterialSystem::MultimaterialSystem()
{
}

MultimaterialSystem::MultimaterialSystem(vector<EulerStateVector> newMaterial1Cells, vector<EulerStateVector> newMaterial2Cells, vector<double> newLevelSetFunction)
{
    material1Cells = newMaterial1Cells;
    material2Cells = newMaterial2Cells;
    levelSetFunction = newLevelSetFunction;
}

MultimaterialSystem::MultimaterialSystem(vector<vector<EulerStateVector> > newMaterial1Cells2D, vector<vector<EulerStateVector> > newMaterial2Cells2D,
                                         vector<vector<double> > newLevelSetFunction2D)
{
    material1Cells2D = newMaterial1Cells2D;
    material2Cells2D = newMaterial2Cells2D;
    levelSetFunction2D = newLevelSetFunction2D;
}

void MultimaterialSystem::setMaterial1Cells(vector<EulerStateVector> newMaterial1Cells)
{
    material1Cells = newMaterial1Cells;
}

void MultimaterialSystem::setMaterial2Cells(vector<EulerStateVector> newMaterial2Cells)
{
    material2Cells = newMaterial2Cells;
}

void MultimaterialSystem::setLevelSetFunction(vector<double> newLevelSetFunction)
{
    levelSetFunction = newLevelSetFunction;
}

void MultimaterialSystem::setMaterial1Cells2D(vector<vector<EulerStateVector> > newMaterial1Cells2D)
{
    material1Cells2D = newMaterial1Cells2D;
}

void MultimaterialSystem::setMaterial2Cells2D(vector<vector<EulerStateVector> > newMaterial2Cells2D)
{
    material2Cells2D = newMaterial2Cells2D;
}

void MultimaterialSystem::setLevelSetFunction2D(vector<vector<double> > newLevelSetFunction2D)
{
    levelSetFunction2D = newLevelSetFunction2D;
}

vector<EulerStateVector> MultimaterialSystem::getMaterial1Cells()
{
    return material1Cells;
}

vector<EulerStateVector> MultimaterialSystem::getMaterial2Cells()
{
    return material2Cells;
}

vector<double> MultimaterialSystem::getLevelSetFunction()
{
    return levelSetFunction;
}

vector<vector<EulerStateVector> > MultimaterialSystem::getMaterial1Cells2D()
{
    return material1Cells2D;
}

vector<vector<EulerStateVector> > MultimaterialSystem::getMaterial2Cells2D()
{
    return material2Cells2D;
}

vector<vector<double> > MultimaterialSystem::getLevelSetFunction2D()
{
    return levelSetFunction2D;
}
