#ifndef MULTIMATERIALSYSTEM_H
#define MULTIMATERIALSYSTEM_H

#include "exactsolver.h"
#include "hllcsolver.h"
using namespace std;

class MultimaterialSystem
{
public:
    MultimaterialSystem();
    MultimaterialSystem(vector<EulerStateVector> newMaterial1Cells, vector<EulerStateVector> newMaterial2Cells, vector<double> newLevelSetFunction);
    MultimaterialSystem(vector<vector<EulerStateVector> > newMaterial1Cells2D, vector<vector<EulerStateVector> > newMaterial2Cells2D, vector<vector<double> > newLevelSetFunction2D);

    void setMaterial1Cells(vector<EulerStateVector> newMaterial1Cells);
    void setMaterial2Cells(vector<EulerStateVector> newMaterial2Cells);
    void setLevelSetFunction(vector<double> newLevelSetFunction);

    void setMaterial1Cells2D(vector<vector<EulerStateVector> > newMaterial1Cells2D);
    void setMaterial2Cells2D(vector<vector<EulerStateVector> > newMaterial2Cells2D);
    void setLevelSetFunction2D(vector<vector<double> > newLevelSetFunction2D);

    vector<EulerStateVector> getMaterial1Cells();
    vector<EulerStateVector> getMaterial2Cells();
    vector<double> getLevelSetFunction();

    vector<vector<EulerStateVector> > getMaterial1Cells2D();
    vector<vector<EulerStateVector> > getMaterial2Cells2D();
    vector<vector<double> > getLevelSetFunction2D();

private:
    vector<EulerStateVector> material1Cells;
    vector<EulerStateVector> material2Cells;
    vector<double> levelSetFunction;

    vector<vector<EulerStateVector> > material1Cells2D;
    vector<vector<EulerStateVector> > material2Cells2D;
    vector<vector<double> > levelSetFunction2D;
};

#endif // MULTIMATERIALSYSTEM_H
