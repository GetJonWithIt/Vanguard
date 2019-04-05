#ifndef ELASTICMULTIPHYSICSSTATEVECTOR_H
#define ELASTICMULTIPHYSICSSTATEVECTOR_H

#include "Elasticity/elasticequationofstate.h"
#include "Elasticity/elasticacoustictensor.h"
using namespace std;

class ElasticMultiphysicsStateVector
{
public:
    ElasticMultiphysicsStateVector();
    ElasticMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                   vector<vector<double> > newMaterial1DistortionTensor, double newMaterial1Entropy, vector<vector<double> > newMaterial2DistortionTensor,
                                   double newMaterial2Entropy, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    vector<double> computeXFluxVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    double computeMaterial1TotalEnergy(HyperelasticMaterialParameters material1Parameters);
    double computeMaterial1SoundSpeed(HyperelasticMaterialParameters material1Parameters, int direction);

    double computeMaterial2TotalEnergy(HyperelasticMaterialParameters material2Parameters);
    double computeMaterial2SoundSpeed(HyperelasticMaterialParameters material2Parameters, int direction);

    double computeTotalDensity();
    vector<vector<double> > computeTotalDistortionTensor();
    double computeTotalEntropy();

    void relaxTotalDensity();
    void relaxTotalDistortionTensor();
    void relaxTotalEntropy(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    void setMaterial1VolumeFraction(double newMaterial1VolumeFraction);
    void setInterfaceXVelocity(double newInterfaceXVelocity);
    void setInterfaceYVelocity(double newInterfaceYVelocity);
    void setInterfaceZVelocity(double newInterfaceZVelocity);

    void setMaterial1Density(double newMaterial1Density);
    void setMaterial1DistortionTensor(vector<vector<double> > newMaterial1DistortionTensor);
    void setMaterial1Entropy(double newMaterial1Entropy);

    void setMaterial2Density(double newMaterial2Density);
    void setMaterial2DistortionTensor(vector<vector<double> > newMaterial2DistortionTensor);
    void setMaterial2Entropy(double newMaterial2Entropy);

    double getMaterial1VolumeFraction();
    double getInterfaceXVelocity();
    double getInterfaceYVelocity();
    double getInterfaceZVelocity();

    double getMaterial1Density();
    vector<vector<double> > getMaterial1DistortionTensor();
    double getMaterial1Entropy();

    double getMaterial2Density();
    vector<vector<double> > getMaterial2DistortionTensor();
    double getMaterial2Entropy();

private:
    double material1VolumeFraction;
    double interfaceXVelocity;
    double interfaceYVelocity;
    double interfaceZVelocity;

    double material1Density;
    vector<vector<double> > material1DistortionTensor;
    double material1Entropy;

    double material2Density;
    vector<vector<double> > material2DistortionTensor;
    double material2Entropy;
};

#endif // ELASTICMULTIPHYSICSSTATEVECTOR_H
