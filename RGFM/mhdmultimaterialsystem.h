#ifndef MHDMULTIMATERIALSYSTEM_H
#define MHDMULTIMATERIALSYSTEM_H

#include "mhdhllcsolver.h"
using namespace std;

class MHDMultimaterialSystem
{
public:
    MHDMultimaterialSystem();
    MHDMultimaterialSystem(vector<MHDStateVector> newMaterial1Cells, vector<MHDStateVector> newMaterial2Cells, vector<double> newLevelSetFunction);
    MHDMultimaterialSystem(vector<vector<MHDStateVector> > newMaterial1Cells2D, vector<vector<MHDStateVector> > newMaterial2Cells2D, vector<vector<double> > newLevelSetFunction2D);

    void setMaterial1Cells(vector<MHDStateVector> newMaterial1Cells);
    void setMaterial2Cells(vector<MHDStateVector> newMaterial2Cells);
    void setLevelSetFunction(vector<double> newLevelSetFunction);

    void setMaterial1Cells2D(vector<vector<MHDStateVector> > newMaterial1Cells2D);
    void setMaterial2Cells2D(vector<vector<MHDStateVector> > newMaterial2Cells2D);
    void setLevelSetFunction2D(vector<vector<double> > newLevelSetFunction2D);

    vector<MHDStateVector> getMaterial1Cells();
    vector<MHDStateVector> getMaterial2Cells();
    vector<double> getLevelSetFunction();

    vector<vector<MHDStateVector> > getMaterial1Cells2D();
    vector<vector<MHDStateVector> > getMaterial2Cells2D();
    vector<vector<double> > getLevelSetFunction2D();

private:
    vector<MHDStateVector> material1Cells;
    vector<MHDStateVector> material2Cells;
    vector<double> levelSetFunction;

    vector<vector<MHDStateVector> > material1Cells2D;
    vector<vector<MHDStateVector> > material2Cells2D;
    vector<vector<double> > levelSetFunction2D;
};

#endif // MHDMULTIMATERIALSYSTEM_H
