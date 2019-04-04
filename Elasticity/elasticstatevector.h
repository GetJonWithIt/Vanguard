#ifndef ELASTICSTATEVECTOR_H
#define ELASTICSTATEVECTOR_H

#include "Mathematics/matrixalgebra.h"
#include "elasticequationofstate.h"
#include "elasticacoustictensor.h"
using namespace std;

class ElasticStateVector
{
public:
    ElasticStateVector();
    ElasticStateVector(double newXVelocity, double newYVelocity, double newZVelocity, vector<vector<double> > newDistortionTensor, double newEntropy,
                       HyperelasticMaterialParameters materialParameters);
    ElasticStateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, vector<vector<double> > newDistortionTensor, double newEntropy);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, HyperelasticMaterialParameters materialParameters);

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(HyperelasticMaterialParameters materialParameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters materialParameters);
    vector<double> computeXFluxVector(HyperelasticMaterialParameters materialParameters);

    static vector<double> computeYFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters materialParameters);
    vector<double> computeYFluxVector(HyperelasticMaterialParameters materialParameters);

    double computeSoundSpeed(HyperelasticMaterialParameters materialParameters, int direction);

    void setDensity(double newDensity);
    void setXVelocity(double newXVelocity);
    void setYVelocity(double newYVelocity);
    void setZVelocity(double newZVelocity);
    void setDistortionTensor(vector<vector<double> > newDistortionTensor);
    void setEntropy(double newEntropy);

    double getDensity();
    double getXVelocity();
    double getYVelocity();
    double getZVelocity();
    vector<vector<double> > getDistortionTensor();
    double getEntropy();

private:
    double density;
    double xVelocity;
    double yVelocity;
    double zVelocity;
    vector<vector<double> > distortionTensor;
    double entropy;
};

#endif // ELASTICSTATEVECTOR_H
