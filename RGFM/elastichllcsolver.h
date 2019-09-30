#ifndef ELASTICHLLCSOLVER_H
#define ELASTICHLLCSOLVER_H

#include "Solvers/elasticsecondordersolver.h"
using namespace std;

class ElasticHLLCSolver
{
public:
    ElasticHLLCSolver();

    static double computeTildeRegionDensity(double density, double waveSpeed, double velocity, double tildeRegionWaveSpeed);

    static double computeLeftTildeRegionDensity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                HyperelasticMaterialParameters material2Parameters);
    static double computeRightTildeRegionDensity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                 HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeLeftTildeRegionConservedVariableVector(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                        HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeRightTildeRegionConservedVariableVector(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                         HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeLeftStarRegionConservedVariableVector(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                       HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeRightStarRegionConservedVariableVector(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                        HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static ElasticStateVector solveTildeX(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                          HyperelasticMaterialParameters material2Parameters);

    static ElasticStateVector solveStarX(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                         HyperelasticMaterialParameters material2Parameters);

    static double computeTildeRegionXVelocity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                              HyperelasticMaterialParameters material2Parameters);

    static vector<vector<double> > computeLeftTildeRegionTotalStressTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                           HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<vector<double> > computeRightTildeRegionTotalStressTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                            HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<vector<double> > computeLeftTildeRegionDistortionTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                          HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<vector<double> > computeRightTildeRegionDistortionTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                           HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static double computeLeftTildeRegionTotalEnergy(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                    HyperelasticMaterialParameters material2Parameters);
    static double computeRightTildeRegionTotalEnergy(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                     HyperelasticMaterialParameters material2Parameters);

    static double computeXStarRegionYVelocity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                              HyperelasticMaterialParameters material2Parameters);
    static double computeXStarRegionZVelocity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                              HyperelasticMaterialParameters material2Parameters);

    static vector<vector<double> > computeLeftStarRegionTotalStressTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                          HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<vector<double> > computeRightStarRegionTotalStressTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                           HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<vector<double> > computeLeftStarRegionDistortionTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                         HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<vector<double> > computeRightStarRegionDistortionTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                          HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static double computeLeftStarRegionTotalEnergy(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                   HyperelasticMaterialParameters material2Parameters);
    static double computeRightStarRegionTotalEnergy(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                    HyperelasticMaterialParameters material2Parameters);

    static double computeLeftWaveSpeed(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                       HyperelasticMaterialParameters material2Parameters);
    static double computeRightWaveSpeed(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                        HyperelasticMaterialParameters material2Parameters);
};

#endif // ELASTICHLLCSOLVER_H
