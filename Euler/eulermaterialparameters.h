#ifndef EULERMATERIALPARAMETERS_H
#define EULERMATERIALPARAMETERS_H


class EulerMaterialParameters
{
public:
    EulerMaterialParameters();
    EulerMaterialParameters(double newAdiabaticIndex);
    EulerMaterialParameters(double newAdiabaticIndex, double newStiffeningParameter);
    EulerMaterialParameters(double newAdiabaticIndex, double newStiffeningParameter, double newSpecificHeatCapacity, double newEnergyOfFormation, double newIgnitionTemperature, double newReactionRate);

    void setAdiabaticIndex(double newAdiabaticIndex);
    void setStiffeningParameter(double newStiffeningParameter);

    void setSpecificHeatCapacity(double newSpecificHeatCapacity);
    void setEnergyOfFormation(double newEnergyOfFormation);
    void setIgnitionTemperature(double newIgnitionTemperature);
    void setReactionRate(double newReactionRate);

    double getAdiabaticIndex();
    double getStiffeningParameter();

    double getSpecificHeatCapacity();
    double getEnergyOfFormation();
    double getIgnitionTemperature();
    double getReactionRate();

private:
    double adiabaticIndex;
    double stiffeningParameter;

    double specificHeatCapacity;
    double energyOfFormation;
    double ignitionTemperature;
    double reactionRate;
};

#endif // EULERMATERIALPARAMETERS_H
