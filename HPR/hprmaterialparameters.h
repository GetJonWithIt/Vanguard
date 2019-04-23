#ifndef HPRMATERIALPARAMETERS_H
#define HPRMATERIALPARAMETERS_H

#include <cmath>
#include <string>
using namespace std;

class HPRMaterialParameters
{
public:
    HPRMaterialParameters();
    HPRMaterialParameters(string newEquationOfState, bool newIsThermal, bool newIsPlastic, double newIdealGasConstant, double newReferenceDensity, double newReferencePressure,
                          double newReferenceTemperature, double newInitialTemperature, double newSpecificHeatCapacity, double newAdiabaticIndex, double newStiffeningParameter,
                          double newReferenceSoundSpeed, double newReferenceGruneisenCoefficient, double newHugoniotSlopeCoefficient, double newReferenceInternalEnergy, double newTransverseWaveSpeed,
                          double newStrainDissipationTime, double newThermalImpulseRelaxationTime, double newDynamicViscosityCoefficient, double newPrandtlNumber,
                          double newElasticPlasticTransitionParameter, double newPowerLawIndex, double newHeatWaveSpeed, double newAlphaParameter, double newBetaParameter, double newGammaParameter);

    void configureStrainDissipationTime();
    void configureThermalImpulseRelaxationTime();
    void configureThermalImpulseRelaxationTime(double heatConductionCoefficient);

    double computeReferenceSoundSpeedSquared();
    double computeTransverseWaveSpeedSquared();
    double computeHeatWaveSpeedSquared();

    double computeStrainDissipationTime();
    double computeHeatConductionCoefficient();
    double computeThermalImpulseRelaxationTime();
    double computeThermalImpulseRelaxationTime(double heatConductionCoefficient);

    void setEquationOfState(string newEquationOfState);
    void setIsThermal(bool newIsThermal);
    void setIsPlastic(bool newIsPlastic);

    void setIdealGasConstant(double newIdealGasConstant);
    void setReferenceDensity(double newReferenceDensity);
    void setReferencePressure(double newReferencePressure);
    void setReferenceTemperature(double newReferenceTemperature);
    void setInitialTemperature(double newInitialTemperature);
    void setSpecificHeatCapacity(double newSpecificHeatCapacity);
    void setAdiabaticIndex(double newAdiabaticIndex);
    void setStiffeningParameter(double newStiffeningParameter);
    void setReferenceSoundSpeed(double newReferenceSoundSpeed);
    void setReferenceGruneisenCoefficient(double newReferenceGruneisenCoefficient);
    void setHugoniotSlopeCoefficient(double newHugoniotSlopeCoefficient);
    void setReferenceInternalEnergy(double newReferenceInternalEnergy);
    void setTransverseWaveSpeed(double newTransverseWaveSpeed);

    void setStrainDissipationTime(double newStrainDissipationTime);
    void setThermalImpulseRelaxationTime(double newThermalImpulseRelaxationTime);
    void setDynamicViscosityCoefficient(double newDynamicViscosityCoefficient);
    void setPrandtlNumber(double newPrandtlNumber);
    void setElasticPlasticTransitionParameter(double newElasticPlasticTransitionParameter);
    void setPowerLawIndex(double newPowerLawIndex);
    void setHeatWaveSpeed(double newHeatWaveSpeed);

    void setAlphaParameter(double newAlphaParameter);
    void setBetaParameter(double newBetaParameter);
    void setGammaParameter(double newGammaParameter);

    string getEquationOfState();
    bool getIsThermal();
    bool getIsPlastic();

    double getIdealGasConstant();
    double getReferenceDensity();
    double getReferencePressure();
    double getReferenceTemperature();
    double getInitialTemperature();
    double getSpecificHeatCapacity();
    double getAdiabaticIndex();
    double getStiffeningParameter();
    double getReferenceSoundSpeed();
    double getReferenceGruneisenCoefficient();
    double getHugoniotSlopeCoefficient();
    double getReferenceInternalEnergy();
    double getTransverseWaveSpeed();

    double getStrainDissipationTime();
    double getThermalImpulseRelaxationTime();
    double getDynamicViscosityCoefficient();
    double getPrandtlNumber();
    double getElasticPlasticTransitionParameter();
    double getPowerLawIndex();
    double getHeatWaveSpeed();

    double getAlphaParameter();
    double getBetaParameter();
    double getGammaParameter();

private:
    string equationOfState;
    bool isThermal;
    bool isPlastic;

    double idealGasConstant;
    double referenceDensity;
    double referencePressure;
    double referenceTemperature;
    double initialTemperature;
    double specificHeatCapacity;
    double adiabaticIndex;
    double stiffeningParameter;
    double referenceSoundSpeed;
    double referenceGruneisenCoefficient;
    double hugoniotSlopeCoefficient;
    double referenceInternalEnergy;
    double transverseWaveSpeed;

    double strainDissipationTime;
    double thermalImpulseRelaxationTime;
    double dynamicViscosityCoefficient;
    double prandtlNumber;
    double elasticPlasticTransitionParameter;
    double powerLawIndex;
    double heatWaveSpeed;

    double alphaParameter;
    double betaParameter;
    double gammaParameter;
};

#endif // HPRMATERIALPARAMETERS_H
