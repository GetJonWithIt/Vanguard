#ifndef MULTIPHYSICSSOLVERS_H
#define MULTIPHYSICSSOLVERS_H

#include "solvers.h"
using namespace std;

class MultiphysicsSolvers
{
public:
    MultiphysicsSolvers();

    static void reinitialiseVolumeFraction(vector<EulerMultiphysicsStateVector> & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static void reinitialiseVolumeFraction(vector<vector<EulerMultiphysicsStateVector> > & currentCells, EulerMaterialParameters material1Parameters,
                                           EulerMaterialParameters material2Parameters);

    static void reinitialiseVolumeFraction(vector<EulerReducedStateVector> & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static void reinitialiseVolumeFraction(vector<vector<EulerReducedStateVector> > & currentCells, EulerMaterialParameters material1Parameters,
                                           EulerMaterialParameters material2Parameters);

    static void reinitialiseVolumeFraction(vector<ElasticMultiphysicsStateVector> & currentCells, HyperelasticMaterialParameters material1Parameters,
                                           HyperelasticMaterialParameters material2Parameters);

    static void reinitialiseVolumeFraction(vector<ElasticReducedStateVector> & currentCells);
    static void reinitialiseVolumeFraction(vector<vector<ElasticReducedStateVector> > & currentCells);

    static void reinitialiseVolumeFraction(vector<HPRReducedStateVector> & currentCells, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void reinitialiseVolumeFraction(vector<vector<HPRReducedStateVector> > & currentCells, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
};

#endif // MULTIPHYSICSSOLVERS_H
