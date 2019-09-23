#ifndef MHDMULTIMATERIALSYSTEM_H
#define MHDMULTIMATERIALSYSTEM_H

#include "mhdhllcsolver.h"
using namespace std;

class MHDMultimaterialSystem
{
public:
    MHDMultimaterialSystem();
    MHDMultimaterialSystem(vector<MHDStateVector> newMaterial1Cells, vector<MHDStateVector> newMaterial2Cells, vector<double> newLevelSetFunction);

    void setMaterial1Cells(vector<MHDStateVector> newMaterial1Cells);
    void setMaterial2Cells(vector<MHDStateVector> newMaterial2Cells);
    void setLevelSetFunction(vector<double> newLevelSetFunction);

    vector<MHDStateVector> getMaterial1Cells();
    vector<MHDStateVector> getMaterial2Cells();
    vector<double> getLevelSetFunction();

private:
    vector<MHDStateVector> material1Cells;
    vector<MHDStateVector> material2Cells;
    vector<double> levelSetFunction;
};

#endif // MHDMULTIMATERIALSYSTEM_H
