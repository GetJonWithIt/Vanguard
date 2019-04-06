#ifndef ELASTICREDUCEDSTATEVECTOR_H
#define ELASTICREDUCEDSTATEVECTOR_H

#include "Elasticity/elasticequationofstate.h"
#include "Elasticity/elasticacoustictensor.h"
using namespace std;

class ElasticReducedStateVector
{
public:
    ElasticReducedStateVector();
    ElasticReducedStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                              vector<vector<double> > newMaterial1DistortionTensor, double newMaterial1Entropy, vector<vector<double> > newMaterial2DistortionTensor, double newMaterial2Entropy,
                              HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    vector<double> computeXFluxVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    vector<double> computeYFluxVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    double computeMaterial1SoundSpeed(HyperelasticMaterialParameters material1Parameters, int direction);
    double computeMaterial2SoundSpeed(HyperelasticMaterialParameters material2Parameters, int direction);

    double computeTotalDensity();
    void relaxTotalDensity();

    void setMaterial1VolumeFraction(double newMaterial1VolumeFraction);
    void setInterfaceXVelocity(double newInterfaceXVelocity);
    void setInterfaceYVelocity(double newInterfaceYVelocity);
    void setInterfaceZVelocity(double newInterfaceZVelocity);

    void setInterfaceDistortionTensor(vector<vector<double> > newInterfaceDistortionTensor);
    void setInterfaceEntropy(double newInterfaceEntropy);

    void setMaterial1Density(double newMaterial1Density);
    void setMaterial2Density(double newMaterial2Density);

    double getMaterial1VolumeFraction();
    double getInterfaceXVelocity();
    double getInterfaceYVelocity();
    double getInterfaceZVelocity();

    vector<vector<double> > getInterfaceDistortionTensor();
    double getInterfaceEntropy();

    double getMaterial1Density();
    double getMaterial2Density();

private:
    double material1VolumeFraction;
    double interfaceXVelocity;
    double interfaceYVelocity;
    double interfaceZVelocity;

    vector<vector<double> > interfaceDistortionTensor;
    double interfaceEntropy;

    double material1Density;
    double material2Density;
};

#endif // ELASTICREDUCEDSTATEVECTOR_H
