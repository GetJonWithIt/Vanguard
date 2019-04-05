#ifndef MHDMATERIALPARAMETERS_H
#define MHDMATERIALPARAMETERS_H

#include <cmath>
using namespace std;

class MHDMaterialParameters
{
public:
    MHDMaterialParameters();

    void setAdiabaticIndex(double newAdiabaticIndex);
    void setStiffeningParameter(double newStiffeningParameter);

    void setHyperbolicWaveSpeed(double newHyperbolicWaveSpeed);
    void setParabolicDamping(double newParabolicDamping);

    double getAdiabaticIndex();
    double getStiffeningParameter();

    double getHyperbolicWaveSpeed();
    double getParabolicDamping();

private:
    double adiabaticIndex;
    double stiffeningParameter;

    double hyperbolicWaveSpeed;
    double parabolicDamping;
};

#endif // MHDMATERIALPARAMETERS_H
