#ifndef HYPERELASTICMATERIALPARAMETERS_H
#define HYPERELASTICMATERIALPARAMETERS_H

#include <cmath>
using namespace std;

class HyperelasticMaterialParameters
{
public:
    HyperelasticMaterialParameters();
    HyperelasticMaterialParameters(double newReferenceMassDensity, double newSoundSpeed, double newShearWaveSpeed, double newSpecificHeatCapacity, double newInitialTemperature,
                                   double newAlphaParameter, double newBetaParameter, double newGammaParameter);

    double computeBulkSoundSpeedSquared();
    double computeShearWaveSpeedSquared();

    void setReferenceMassDensity(double newReferenceMassDensity);
    void setSoundSpeed(double newSoundSpeed);
    void setShearWaveSpeed(double newShearWaveSpeed);
    void setSpecificHeatCapacity(double newSpecificHeatCapacity);
    void setInitialTemperature(double newInitialTemperature);

    void setAlphaParameter(double newAlphaParameter);
    void setBetaParameter(double newBetaParameter);
    void setGammaParameter(double newGammaParameter);

    double getReferenceMassDensity();
    double getSoundSpeed();
    double getShearWaveSpeed();
    double getSpecificHeatCapacity();
    double getInitialTemperature();

    double getAlphaParameter();
    double getBetaParameter();
    double getGammaParameter();

private:
    double referenceMassDensity;
    double soundSpeed;
    double shearWaveSpeed;
    double specificHeatCapacity;
    double initialTemperature;

    double alphaParameter;
    double betaParameter;
    double gammaParameter;
};

#endif // HYPERELASTICMATERIALPARAMETERS_H
