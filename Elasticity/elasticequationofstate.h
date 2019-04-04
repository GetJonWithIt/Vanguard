#ifndef ELASTICEQUATIONOFSTATE_H
#define ELASTICEQUATIONOFSTATE_H

#include "Mathematics/tensoralgebra.h"
#include "hyperelasticmaterialparameters.h"
using namespace std;

class ElasticEquationOfState
{
public:
    ElasticEquationOfState();

    static double computeHydrodynamicInternalEnergyFirstComponent(double fingerTensorThirdInvariant, double bulkSoundSpeedSquared, double alphaParameter);
    static double computeHydrodynamicInternalEnergySecondComponent(double fingerTensorThirdInvariant, double entropy, double specificHeatCapacity, double initialTemperature,
                                                                   double gammaParameter);
    static double computeShearDeformationInternalEnergy(double fingerTensorFirstInvariant, double fingerTensorSecondInvariant, double fingerTensorThirdInvariant,
                                                        double shearWaveSpeedSquared, double betaParameter);

    static vector<vector<double> > computeHydrodynamicInternalEnergyFirstComponentDerivative(double fingerTensorThirdInvariant, double bulkSoundSpeedSquared, double alphaParameter);
    static vector<vector<double> > computeHydrodynamicInternalEnergySecondComponentDerivative(double fingerTensorThirdInvariant, double entropy, double specificHeatCapacity,
                                                                                              double initialTemperature, double gammaParameter);
    static vector<vector<double> > computeShearDeformationInternalEnergyDerivative(vector<vector<double> > fingerTensor, double fingerTensorFirstInvariant, double fingerTensorSecondInvariant,
                                                                                   double fingerTensorThirdInvariant, double shearWaveSpeedSquared, double betaParameter);

    static double computeTotalEnergy(vector<vector<double> > distortionTensor, double entropy, double xVelocity, double yVelocity, double zVelocity,
                                     HyperelasticMaterialParameters materialParameters);
    static double computeEntropy(double totalEnergy, vector<vector<double> > distortionTensor, double xVelocity, double yVelocity, double zVelocity,
                                 HyperelasticMaterialParameters materialParameters);
    static double computePressure(double density, vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters);

    static vector<vector<double> > computeTotalStressTensor(double density, vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters);
};

#endif // ELASTICEQUATIONOFSTATE_H
