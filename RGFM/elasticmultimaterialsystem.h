#ifndef ELASTICMULTIMATERIALSYSTEM_H
#define ELASTICMULTIMATERIALSYSTEM_H

#include "elastichllcsolver.h"
using namespace std;

class ElasticMultimaterialSystem
{
public:
    ElasticMultimaterialSystem();
    ElasticMultimaterialSystem(vector<ElasticStateVector> newMaterial1Cells, vector<ElasticStateVector> newMaterial2Cells, vector<double> newLevelSetFunction);

    void setMaterial1Cells(vector<ElasticStateVector> newMaterial1Cells);
    void setMaterial2Cells(vector<ElasticStateVector> newMaterial2Cells);
    void setLevelSetFunction(vector<double> newLevelSetFunction);

    vector<ElasticStateVector> getMaterial1Cells();
    vector<ElasticStateVector> getMaterial2Cells();
    vector<double> getLevelSetFunction();

private:
    vector<ElasticStateVector> material1Cells;
    vector<ElasticStateVector> material2Cells;
    vector<double> levelSetFunction;
};

#endif // ELASTICMULTIMATERIALSYSTEM_H
