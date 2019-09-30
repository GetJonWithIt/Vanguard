#ifndef ELASTICRGFMSOLVER_H
#define ELASTICRGFMSOLVER_H

#include "elasticmultimaterialsystem.h"
#include "rgfmsolver.h"
using namespace std;

class ElasticRGFMSolver
{
public:
    ElasticRGFMSolver();

    static ElasticMultimaterialSystem applyTildeRGFMBoundaryConditions(ElasticMultimaterialSystem multimaterialSystem, HyperelasticMaterialParameters material1Parameters,
                                                                       HyperelasticMaterialParameters material2Parameters);

    static ElasticMultimaterialSystem applyStarRGFMBoundaryConditions(ElasticMultimaterialSystem multimaterialSystem, HyperelasticMaterialParameters material1Parameters,
                                                                      HyperelasticMaterialParameters material2Parameters);

    static vector<double> updateLevelSetFunction(vector<double> levelSetFunction, double cellSpacing, double timeStep, vector<ElasticStateVector> material1Cells,
                                                 vector<ElasticStateVector> material2Cells);

    static ElasticMultimaterialSystem solveTilde(ElasticMultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                 int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                 HyperelasticMaterialParameters material2Parameters);

    static ElasticMultimaterialSystem solveStar(ElasticMultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                HyperelasticMaterialParameters material2Parameters);
};

#endif // ELASTICRGFMSOLVER_H
