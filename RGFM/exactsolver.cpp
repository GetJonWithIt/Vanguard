#include "exactsolver.h"

ExactSolver::ExactSolver()
{
}

double ExactSolver::computeACoefficient(double density, double adiabaticIndex)
{
    return 2.0 / (density * (adiabaticIndex + 1.0));
}

double ExactSolver::computeBCoefficient(double pressure, double adiabaticIndex, double stiffeningParameter)
{
    return (pressure * ((adiabaticIndex - 1.0) / (adiabaticIndex + 1.0))) + ((2.0 * adiabaticIndex * stiffeningParameter) / (adiabaticIndex + 1.0));
}

double ExactSolver::computeStarRegionSoundSpeed(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    double pressure = stateVector.getPressure();
    double soundSpeed = stateVector.computeSoundSpeed(materialParameters);

    return soundSpeed * pow((starRegionPressure + stiffeningParameter) / (pressure + stiffeningParameter), (adiabaticIndex - 1.0) / (2.0 * adiabaticIndex));
}

vector<double> ExactSolver::computeRarefactionXPrimitiveVariableVector(double waveSpeed, double soundSpeed, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    vector<double> rarefactionPrimitiveVariableVector(6);

    double density = stateVector.getDensity();
    double xVelocity = stateVector.getXVelocity();
    double pressure = stateVector.getPressure();

    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    double waveSpeedCoefficient = (2.0 / (adiabaticIndex + 1.0)) + (adiabaticIndex - 1.0) * ((xVelocity - waveSpeed) / ((adiabaticIndex + 1.0) * soundSpeed));

    rarefactionPrimitiveVariableVector[0] = density * pow(waveSpeedCoefficient, 2.0 / (adiabaticIndex - 1.0));
    rarefactionPrimitiveVariableVector[1] = 2.0 * (soundSpeed + (adiabaticIndex - 1.0) * (xVelocity / 2.0) + waveSpeed) / (adiabaticIndex + 1.0);
    rarefactionPrimitiveVariableVector[2] = 0.0;
    rarefactionPrimitiveVariableVector[3] = 0.0;
    rarefactionPrimitiveVariableVector[4] = ((pressure + stiffeningParameter) * pow(waveSpeedCoefficient, (2.0 * adiabaticIndex) / (adiabaticIndex - 1.0))) - stiffeningParameter;

    rarefactionPrimitiveVariableVector[5] = density;

    return rarefactionPrimitiveVariableVector;
}

vector<double> ExactSolver::computeRarefactionYPrimitiveVariableVector(double waveSpeed, double soundSpeed, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    vector<double> rarefactionPrimitiveVariableVector(6);

    double density = stateVector.getDensity();
    double yVelocity = stateVector.getYVelocity();
    double pressure = stateVector.getPressure();

    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    double waveSpeedCoefficient = (2.0 / (adiabaticIndex + 1.0)) + (adiabaticIndex - 1.0) * ((yVelocity - waveSpeed) / ((adiabaticIndex + 1.0) * soundSpeed));

    rarefactionPrimitiveVariableVector[0] = density * pow(waveSpeedCoefficient, 2.0 / (adiabaticIndex - 1.0));
    rarefactionPrimitiveVariableVector[1] = 0.0;
    rarefactionPrimitiveVariableVector[2] = 2.0 * (soundSpeed + (adiabaticIndex - 1.0) * (yVelocity / 2.0) + waveSpeed) / (adiabaticIndex + 1.0);
    rarefactionPrimitiveVariableVector[3] = 0.0;
    rarefactionPrimitiveVariableVector[4] = ((pressure + stiffeningParameter) * pow(waveSpeedCoefficient, (2.0 * adiabaticIndex) / (adiabaticIndex - 1.0))) - stiffeningParameter;

    rarefactionPrimitiveVariableVector[5] = density;

    return rarefactionPrimitiveVariableVector;
}

double ExactSolver::computeStarRegionShockDensity(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();

    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    double numerator = ((starRegionPressure + stiffeningParameter) / (pressure + stiffeningParameter)) + ((adiabaticIndex - 1.0) / (adiabaticIndex + 1.0));
    double denominator = ((adiabaticIndex - 1.0) / (adiabaticIndex + 1.0)) * ((starRegionPressure + stiffeningParameter) / (pressure + stiffeningParameter)) + 1.0;

    return density * (numerator / denominator);
}

double ExactSolver::computeStarRegionRarefactionDensity(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();

    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    return density * pow((starRegionPressure + stiffeningParameter) / (pressure + stiffeningParameter), 1.0 / adiabaticIndex);
}

double ExactSolver::computeStarRegionXVelocity(double starRegionPressure, EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                               EulerMaterialParameters material2Parameters)
{
    return 0.5 * (leftStateVector.getXVelocity() + rightStateVector.getXVelocity() + computeWaveJumpFunctionComponent(starRegionPressure, rightStateVector, material2Parameters) -
                  computeWaveJumpFunctionComponent(starRegionPressure, leftStateVector, material1Parameters));
}

double ExactSolver::computeStarRegionYVelocity(double starRegionPressure, EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                               EulerMaterialParameters material2Parameters)
{
    return 0.5 * (topStateVector.getYVelocity() + bottomStateVector.getYVelocity() + computeWaveJumpFunctionComponent(starRegionPressure, bottomStateVector, material2Parameters) -
                  computeWaveJumpFunctionComponent(starRegionPressure, topStateVector, material1Parameters));
}

double ExactSolver::computeStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                               EulerMaterialParameters material2Parameters)
{
    double currentPressure = 0.5 * (leftStateVector.getPressure() + rightStateVector.getPressure());
    double newPressure = 3.0 * currentPressure;

    int iterationCount = 0;

    while (computePressureChange(newPressure, currentPressure) >= pow(10.0, -8.0) && iterationCount < 10000)
    {
        currentPressure = newPressure;
        newPressure = currentPressure - (computeXWaveJumpFunction(currentPressure, leftStateVector, rightStateVector, material1Parameters, material2Parameters) /
                                         computeWaveJumpFunctionDerivative(currentPressure, leftStateVector, rightStateVector, material1Parameters, material2Parameters));

        if (newPressure < pow(10.0, -8.0))
        {
            newPressure = pow(10.0, -8.0);
        }

        iterationCount += 1;
    }

    return newPressure;
}

double ExactSolver::computeStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                               EulerMaterialParameters material2Parameters)
{
    double currentPressure = 0.5 * (topStateVector.getPressure() + bottomStateVector.getPressure());
    double newPressure = 3.0 * currentPressure;

    int iterationCount = 0;

    while (computePressureChange(newPressure, currentPressure) >= pow(10.0, -8.0) && iterationCount < 10000)
    {
        currentPressure = newPressure;
        newPressure = currentPressure - (computeYWaveJumpFunction(currentPressure, topStateVector, bottomStateVector, material1Parameters, material2Parameters) /
                                         computeWaveJumpFunctionDerivative(currentPressure, topStateVector, bottomStateVector, material1Parameters, material2Parameters));

        if (newPressure < pow(10.0, -8.0))
        {
            newPressure = pow(10.0, -8.0);
        }

        iterationCount += 1;
    }

    return newPressure;
}

double ExactSolver::computeWaveJumpFunctionComponent(double newPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();
    double soundSpeed = stateVector.computeSoundSpeed(materialParameters);

    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    if (newPressure > pressure)
    {
        double coefficient = sqrt(computeACoefficient(density, adiabaticIndex) / (newPressure + computeBCoefficient(pressure, adiabaticIndex, stiffeningParameter)));

        return (newPressure - pressure) * coefficient;
    }
    else
    {
        double coefficient = pow((newPressure + stiffeningParameter) / (pressure + stiffeningParameter), (adiabaticIndex - 1.0) / (2.0 * adiabaticIndex));

        return 2.0 * (soundSpeed / (adiabaticIndex - 1.0)) * (coefficient - 1.0);
    }
}

double ExactSolver::computeXWaveJumpFunction(double newPressure, EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters)
{
    return computeWaveJumpFunctionComponent(newPressure, leftStateVector, material1Parameters) + computeWaveJumpFunctionComponent(newPressure, rightStateVector, material2Parameters) +
            rightStateVector.getXVelocity() - leftStateVector.getXVelocity();
}

double ExactSolver::computeYWaveJumpFunction(double newPressure, EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters)
{
    return computeWaveJumpFunctionComponent(newPressure, topStateVector, material1Parameters) + computeWaveJumpFunctionComponent(newPressure, bottomStateVector, material2Parameters) +
            bottomStateVector.getYVelocity() - topStateVector.getYVelocity();
}

double ExactSolver::computeWaveJumpFunctionComponentDerivative(double newPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();
    double soundSpeed = stateVector.computeSoundSpeed(materialParameters);

    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    if (newPressure > pressure)
    {
        double coefficient = sqrt(computeACoefficient(density, adiabaticIndex) / (newPressure + computeBCoefficient(pressure, adiabaticIndex, stiffeningParameter)));

        return (1.0 - ((newPressure - pressure) / (2.0 * (newPressure + computeBCoefficient(pressure, adiabaticIndex, stiffeningParameter))))) * coefficient;
    }
    else
    {
        double coefficient = pow((newPressure + stiffeningParameter) / (pressure + stiffeningParameter), -(adiabaticIndex + 1.0) / (2.0 * adiabaticIndex));

        return (coefficient * soundSpeed) / (adiabaticIndex * (pressure + stiffeningParameter));
    }
}

double ExactSolver::computeWaveJumpFunctionDerivative(double newPressure, EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                      EulerMaterialParameters material2Parameters)
{
    return computeWaveJumpFunctionComponentDerivative(newPressure, leftStateVector, material1Parameters) + computeWaveJumpFunctionComponentDerivative(newPressure, rightStateVector,
                                                                                                                                                      material2Parameters);
}

double ExactSolver::computePressureChange(double oldPressure, double newPressure)
{
    return 2.0 * (abs(oldPressure - newPressure) / abs(oldPressure + newPressure));
}

double ExactSolver::computeShockSpeed(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();

    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    return sqrt((starRegionPressure + computeBCoefficient(pressure, adiabaticIndex, stiffeningParameter)) / computeACoefficient(density, adiabaticIndex));
}

double ExactSolver::computeStarRegionDensity(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    double pressure = stateVector.getPressure();

    if (starRegionPressure < pressure)
    {
        return computeStarRegionRarefactionDensity(starRegionPressure, stateVector, materialParameters);
    }
    else
    {
        return computeStarRegionShockDensity(starRegionPressure, stateVector, materialParameters);
    }
}

EulerStateVector ExactSolver::solveX(double position, double time, double interfaceLocation, EulerStateVector leftStateVector, EulerStateVector rightStateVector,
                                     EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double waveSpeed = (position - interfaceLocation) / time;

    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftPressure = leftStateVector.getPressure();

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightPressure = rightStateVector.getPressure();

    double leftSoundSpeed = leftStateVector.computeSoundSpeed(material1Parameters);
    double rightSoundSpeed = rightStateVector.computeSoundSpeed(material2Parameters);

    double starRegionPressure = computeStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionVelocity = computeStarRegionXVelocity(starRegionPressure, leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (waveSpeed < starRegionVelocity)
    {
        if (starRegionPressure < leftPressure)
        {
            if (waveSpeed < (leftXVelocity - leftSoundSpeed))
            {
                return leftStateVector;
            }
            else
            {
                double leftAcousticWaveSpeed = starRegionVelocity - computeStarRegionSoundSpeed(starRegionPressure, leftStateVector, material1Parameters);

                if (waveSpeed < leftAcousticWaveSpeed)
                {
                    return EulerStateVector(computeRarefactionXPrimitiveVariableVector(waveSpeed, leftSoundSpeed, leftStateVector, material1Parameters));
                }
                else
                {
                    return EulerStateVector(computeStarRegionRarefactionDensity(starRegionPressure, leftStateVector, material1Parameters), starRegionVelocity, 0.0, 0.0, starRegionPressure);
                }
            }
        }
        else
        {
            double leftShockSpeed = leftXVelocity - (computeShockSpeed(starRegionPressure, leftStateVector, material1Parameters) / leftDensity);

            if (waveSpeed < leftShockSpeed)
            {
                return leftStateVector;
            }
            else
            {
                return EulerStateVector(computeStarRegionShockDensity(starRegionPressure, leftStateVector, material1Parameters), starRegionVelocity, 0.0, 0.0, starRegionPressure);
            }
        }
    }
    else
    {
        if (starRegionPressure < rightPressure)
        {
            if (waveSpeed > (rightXVelocity + rightSoundSpeed))
            {
                return rightStateVector;
            }
            else
            {
                double rightAcousticWaveSpeed = starRegionVelocity + computeStarRegionSoundSpeed(starRegionPressure, rightStateVector, material2Parameters);

                if (waveSpeed > rightAcousticWaveSpeed)
                {
                    return EulerStateVector(computeRarefactionXPrimitiveVariableVector(waveSpeed, -rightSoundSpeed, rightStateVector, material2Parameters));
                }
                else
                {
                    return EulerStateVector(computeStarRegionRarefactionDensity(starRegionPressure, rightStateVector, material2Parameters), starRegionVelocity, 0.0, 0.0, starRegionPressure);
                }
            }
        }
        else
        {
            double rightShockSpeed = rightXVelocity + (computeShockSpeed(starRegionPressure, rightStateVector, material2Parameters) / rightDensity);

            if (waveSpeed > rightShockSpeed)
            {
                return rightStateVector;
            }
            else
            {
                return EulerStateVector(computeStarRegionShockDensity(starRegionPressure, rightStateVector, material2Parameters), starRegionVelocity, 0.0, 0.0, starRegionPressure);
            }
        }
    }
}

EulerStateVector ExactSolver::solveY(double position, double time, double interfaceLocation, EulerStateVector topStateVector, EulerStateVector bottomStateVector,
                                     EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double waveSpeed = (position - interfaceLocation) / time;

    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();
    double topPressure = topStateVector.getPressure();

    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomPressure = bottomStateVector.getPressure();

    double topSoundSpeed = topStateVector.computeSoundSpeed(material1Parameters);
    double bottomSoundSpeed = bottomStateVector.computeSoundSpeed(material2Parameters);

    double starRegionPressure = computeStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionVelocity = computeStarRegionYVelocity(starRegionPressure, topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (waveSpeed < starRegionVelocity)
    {
        if (starRegionPressure < topPressure)
        {
            if (waveSpeed < (topYVelocity - topSoundSpeed))
            {
                return topStateVector;
            }
            else
            {
                double topAcousticWaveSpeed = starRegionVelocity - computeStarRegionSoundSpeed(starRegionPressure, topStateVector, material1Parameters);

                if (waveSpeed < topAcousticWaveSpeed)
                {
                    return EulerStateVector(computeRarefactionYPrimitiveVariableVector(waveSpeed, topSoundSpeed, topStateVector, material1Parameters));
                }
                else
                {
                    return EulerStateVector(computeStarRegionRarefactionDensity(starRegionPressure, topStateVector, material1Parameters), 0.0, starRegionVelocity, 0.0, starRegionPressure);
                }
            }
        }
        else
        {
            double topShockSpeed = topYVelocity - (computeShockSpeed(starRegionPressure, topStateVector, material1Parameters) / topDensity);

            if (waveSpeed < topShockSpeed)
            {
                return topStateVector;
            }
            else
            {
                return EulerStateVector(computeStarRegionShockDensity(starRegionPressure, topStateVector, material1Parameters), 0.0, starRegionVelocity, 0.0, starRegionPressure);
            }
        }
    }
    else
    {
        if (starRegionPressure < bottomPressure)
        {
            if (waveSpeed > (bottomYVelocity + bottomSoundSpeed))
            {
                return bottomStateVector;
            }
            else
            {
                double bottomAcousticWaveSpeed = starRegionVelocity + computeStarRegionSoundSpeed(starRegionPressure, bottomStateVector, material2Parameters);

                if (waveSpeed > bottomAcousticWaveSpeed)
                {
                    return EulerStateVector(computeRarefactionYPrimitiveVariableVector(waveSpeed, -bottomSoundSpeed, bottomStateVector, material2Parameters));
                }
                else
                {
                    return EulerStateVector(computeStarRegionRarefactionDensity(starRegionPressure, bottomStateVector, material2Parameters), 0.0, starRegionVelocity, 0.0, starRegionPressure);
                }
            }
        }
        else
        {
            double bottomShockSpeed = bottomYVelocity + (computeShockSpeed(starRegionPressure, bottomStateVector, material2Parameters) / bottomDensity);

            if (waveSpeed > bottomShockSpeed)
            {
                return bottomStateVector;
            }
            else
            {
                return EulerStateVector(computeStarRegionShockDensity(starRegionPressure, bottomStateVector, material2Parameters), 0.0, starRegionVelocity, 0.0, starRegionPressure);
            }
        }
    }
}
