#include "hllcsolver.h"

HLLCSolver::HLLCSolver()
{
}

double HLLCSolver::computeStarRegionDensity(double density, double waveSpeed, double velocity, double starRegionWaveSpeed)
{
    return density * ((waveSpeed - velocity) / (waveSpeed - starRegionWaveSpeed));
}

double HLLCSolver::computeLeftStarRegionDensity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                EulerMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity =  leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();

    return computeStarRegionDensity(leftDensity, leftWaveSpeed, leftXVelocity, starRegionWaveSpeed);
}

double HLLCSolver::computeRightStarRegionDensity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                 EulerMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();

    return computeStarRegionDensity(rightDensity, rightWaveSpeed, rightXVelocity, starRegionWaveSpeed);
}

double HLLCSolver::computeTopStarRegionDensity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                               EulerMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();

    return computeStarRegionDensity(topDensity, topWaveSpeed, topYVelocity, starRegionWaveSpeed);
}

double HLLCSolver::computeBottomStarRegionDensity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                  EulerMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    return computeStarRegionDensity(bottomDensity, bottomWaveSpeed, bottomYVelocity, starRegionWaveSpeed);
}

vector<double> HLLCSolver::computeStarRegionXConservedVariableVector(EulerStateVector stateVector, double waveSpeed, double starRegionWaveSpeed, EulerMaterialParameters materialParameters)
{
    vector<double> starRegionConservedVariableVector(6);

    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();

    double xVelocity = stateVector.getXVelocity();
    double yVelocity = stateVector.getYVelocity();
    double zVelocity = stateVector.getZVelocity();

    double totalEnergy = stateVector.computeTotalEnergy(materialParameters);

    double starRegionDensity = computeStarRegionDensity(density, waveSpeed, xVelocity, starRegionWaveSpeed);
    double energyTerm = (starRegionWaveSpeed - xVelocity) * (starRegionWaveSpeed + (pressure / (density * (waveSpeed - xVelocity))));

    starRegionConservedVariableVector[0] = starRegionDensity;
    starRegionConservedVariableVector[1] = starRegionDensity * starRegionWaveSpeed;
    starRegionConservedVariableVector[2] = starRegionDensity * yVelocity;
    starRegionConservedVariableVector[3] = starRegionDensity * zVelocity;
    starRegionConservedVariableVector[4] = starRegionDensity * ((totalEnergy / density) + energyTerm);

    starRegionConservedVariableVector[5] = starRegionDensity;

    return starRegionConservedVariableVector;
}

vector<double> HLLCSolver::computeStarRegionYConservedVariableVector(EulerStateVector stateVector, double waveSpeed, double starRegionWaveSpeed, EulerMaterialParameters materialParameters)
{
    vector<double> starRegionConservedVariableVector(6);

    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();

    double xVelocity = stateVector.getXVelocity();
    double yVelocity = stateVector.getYVelocity();
    double zVelocity = stateVector.getZVelocity();

    double totalEnergy = stateVector.computeTotalEnergy(materialParameters);

    double starRegionDensity = computeStarRegionDensity(density, waveSpeed, yVelocity, starRegionWaveSpeed);
    double energyTerm = (starRegionWaveSpeed - yVelocity) * (starRegionWaveSpeed + (pressure / (density * (waveSpeed - yVelocity))));

    starRegionConservedVariableVector[0] = starRegionDensity;
    starRegionConservedVariableVector[1] = starRegionDensity * xVelocity;
    starRegionConservedVariableVector[2] = starRegionDensity * starRegionWaveSpeed;
    starRegionConservedVariableVector[3] = starRegionDensity * zVelocity;
    starRegionConservedVariableVector[4] = starRegionDensity * ((totalEnergy / density) + energyTerm);

    starRegionConservedVariableVector[5] = starRegionDensity;

    return starRegionConservedVariableVector;
}

EulerStateVector HLLCSolver::solveX(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0 <= leftWaveSpeed)
    {
        return leftStateVector;
    }
    else if (rightWaveSpeed <= 0)
    {
        return rightStateVector;
    }
    else if (0 <= starRegionWaveSpeed)
    {
        vector<double> starRegionConservedVariableVector = computeStarRegionXConservedVariableVector(leftStateVector, leftWaveSpeed, starRegionWaveSpeed, material1Parameters);

        EulerStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material1Parameters);

        return solution;
    }
    else
    {
        vector<double> starRegionConservedVariableVector = computeStarRegionXConservedVariableVector(rightStateVector, rightWaveSpeed, starRegionWaveSpeed, material2Parameters);

        EulerStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material2Parameters);

        return solution;
    }
}

EulerStateVector HLLCSolver::solveY(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0 <= topWaveSpeed)
    {
        return topStateVector;
    }
    else if (bottomWaveSpeed <= 0)
    {
        return bottomStateVector;
    }
    else if (0 <= starRegionWaveSpeed)
    {
        vector<double> starRegionConservedVariableVector = computeStarRegionYConservedVariableVector(topStateVector, topWaveSpeed, starRegionWaveSpeed, material1Parameters);

        EulerStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material1Parameters);

        return solution;
    }
    else
    {
        vector<double> starRegionConservedVariableVector = computeStarRegionYConservedVariableVector(bottomStateVector, bottomWaveSpeed, starRegionWaveSpeed, material2Parameters);

        EulerStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material2Parameters);

        return solution;
    }
}

double HLLCSolver::computePVRSStarRegionXVelocity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                  EulerMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftPressure = leftStateVector.getPressure();

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightPressure = rightStateVector.getPressure();

    double leftSoundSpeed = leftStateVector.computeSoundSpeed(material1Parameters);
    double rightSoundSpeed = rightStateVector.computeSoundSpeed(material2Parameters);

    double linearisedDensitySoundSpeedProduct = 0.25 * (leftDensity + rightDensity) * (leftSoundSpeed + rightSoundSpeed);

    return (0.5 * (leftXVelocity + rightXVelocity)) + (0.5 * (leftPressure - rightPressure) / linearisedDensitySoundSpeedProduct);
}

double HLLCSolver::computePVRSStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                  EulerMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftPressure = leftStateVector.getPressure();

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightPressure = rightStateVector.getPressure();

    double leftSoundSpeed = leftStateVector.computeSoundSpeed(material1Parameters);
    double rightSoundSpeed = rightStateVector.computeSoundSpeed(material2Parameters);

    double linearisedDensitySoundSpeedProduct = 0.25 * (leftDensity + rightDensity) * (leftSoundSpeed + rightSoundSpeed);

    return (0.5 * (leftPressure + rightPressure)) + (0.5 * (leftXVelocity - rightXVelocity) * linearisedDensitySoundSpeedProduct);
}

double HLLCSolver::computePVRSStarRegionYVelocity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                  EulerMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();
    double topPressure = topStateVector.getPressure();

    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomPressure = bottomStateVector.getPressure();

    double topSoundSpeed = topStateVector.computeSoundSpeed(material1Parameters);
    double bottomSoundSpeed = bottomStateVector.computeSoundSpeed(material2Parameters);

    double linearisedDensitySoundSpeedProduct = 0.25 * (topDensity + bottomDensity) * (topSoundSpeed + bottomSoundSpeed);

    return (0.5 * (topYVelocity + bottomYVelocity)) + (0.5 * (topPressure - bottomPressure) / linearisedDensitySoundSpeedProduct);
}

double HLLCSolver::computePVRSStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                  EulerMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();
    double topPressure = topStateVector.getPressure();

    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomPressure = bottomStateVector.getPressure();

    double topSoundSpeed = topStateVector.computeSoundSpeed(material1Parameters);
    double bottomSoundSpeed = bottomStateVector.computeSoundSpeed(material2Parameters);

    double linearisedDensitySoundSpeedProduct = 0.25 * (topDensity + bottomDensity) * (topSoundSpeed + bottomSoundSpeed);

    return (0.5 * (topPressure + bottomPressure)) + (0.5 * (topYVelocity - bottomYVelocity) * linearisedDensitySoundSpeedProduct);
}

double HLLCSolver::computeTwoRarefactionStarRegionXVelocity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                            EulerMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftPressure = leftStateVector.getPressure();

    double rightXVelocity = rightStateVector.getXVelocity();
    double rightPressure = rightStateVector.getPressure();

    double leftSoundSpeed = leftStateVector.computeSoundSpeed(material1Parameters);
    double rightSoundSpeed = rightStateVector.computeSoundSpeed(material2Parameters);

    double material1AdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double material2AdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double averageAdiabaticIndex = 0.5 * (material1AdiabaticIndex + material2AdiabaticIndex);

    double exponent = (averageAdiabaticIndex - 1.0) / (2.0 * averageAdiabaticIndex);
    double pressureRatio = pow((leftPressure / rightPressure), exponent);
    double coefficient = 2.0 / (averageAdiabaticIndex - 1.0);

    return ((pressureRatio * (leftXVelocity / leftSoundSpeed)) + (rightXVelocity / rightSoundSpeed) + (coefficient * (pressureRatio - 1.0))) / ((pressureRatio / leftSoundSpeed) +
                                                                                                                                                (1.0 / rightSoundSpeed));
}

double HLLCSolver::computeTwoRarefactionStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                            EulerMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftPressure = leftStateVector.getPressure();

    double rightXVelocity = rightStateVector.getXVelocity();
    double rightPressure = rightStateVector.getPressure();

    double leftSoundSpeed = leftStateVector.computeSoundSpeed(material1Parameters);
    double rightSoundSpeed = rightStateVector.computeSoundSpeed(material2Parameters);

    double material1AdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double material2AdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double averageAdiabaticIndex = 0.5 * (material1AdiabaticIndex + material2AdiabaticIndex);

    double coefficient = 0.5 * (averageAdiabaticIndex - 1.0);
    double starRegionVelocity = computeTwoRarefactionStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double exponent = 2.0 * (averageAdiabaticIndex / (averageAdiabaticIndex - 1.0));

    double leftVelocityDifference = 1.0 + (coefficient * ((leftXVelocity - starRegionVelocity) / leftSoundSpeed));
    double rightVelocityDifference = 1.0 + (coefficient * ((starRegionVelocity - rightXVelocity) / rightSoundSpeed));

    return 0.5 * ((leftPressure * pow(leftVelocityDifference, exponent)) + (rightPressure * pow(rightVelocityDifference, exponent)));
}

double HLLCSolver::computeTwoRarefactionStarRegionYVelocity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                            EulerMaterialParameters material2Parameters)
{
    double topYVelocity = topStateVector.getYVelocity();
    double topPressure = topStateVector.getPressure();

    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomPressure = bottomStateVector.getPressure();

    double topSoundSpeed = topStateVector.computeSoundSpeed(material1Parameters);
    double bottomSoundSpeed = bottomStateVector.computeSoundSpeed(material2Parameters);

    double material1AdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double material2AdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double averageAdiabaticIndex = 0.5 * (material1AdiabaticIndex + material2AdiabaticIndex);

    double exponent = (averageAdiabaticIndex - 1.0) / (2.0 * averageAdiabaticIndex);
    double pressureRatio = pow((topPressure / bottomPressure), exponent);
    double coefficient = 2.0 / (averageAdiabaticIndex - 1.0);

    return ((pressureRatio * (topYVelocity / topSoundSpeed)) + (bottomYVelocity / bottomSoundSpeed) + (coefficient * (pressureRatio - 1.0))) / ((pressureRatio / topSoundSpeed) +
                                                                                                                                                (1.0 / bottomSoundSpeed));
}

double HLLCSolver::computeTwoRarefactionStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                            EulerMaterialParameters material2Parameters)
{
    double topYVelocity = topStateVector.getYVelocity();
    double topPressure = topStateVector.getPressure();

    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomPressure = bottomStateVector.getPressure();

    double topSoundSpeed = topStateVector.computeSoundSpeed(material1Parameters);
    double bottomSoundSpeed = bottomStateVector.computeSoundSpeed(material2Parameters);

    double material1AdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double material2AdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double averageAdiabaticIndex = 0.5 * (material1AdiabaticIndex + material2AdiabaticIndex);

    double coefficient = 0.5 * (averageAdiabaticIndex - 1.0);
    double starRegionVelocity = computeTwoRarefactionStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double exponent = 2.0 * (averageAdiabaticIndex / (averageAdiabaticIndex - 1.0));

    double topVelocityDifference = 1.0 + (coefficient * ((topYVelocity - starRegionVelocity) / topSoundSpeed));
    double bottomVelocityDifference = 1.0 + (coefficient * ((starRegionVelocity - bottomYVelocity) / bottomSoundSpeed));

    return 0.5 * ((topPressure * pow(topVelocityDifference, exponent)) + (bottomPressure * pow(bottomVelocityDifference, exponent)));
}

double HLLCSolver::computeTwoShockStarRegionXVelocity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                      EulerMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftPressure = leftStateVector.getPressure();

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightPressure = rightStateVector.getPressure();

    double material1AdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double material2AdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double averageAdiabaticIndex = 0.5 * (material1AdiabaticIndex + material2AdiabaticIndex);

    double PVRSStarRegionPressure = computePVRSStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionPressure = computeTwoShockStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double ACoefficient = 2.0 / (averageAdiabaticIndex + 1.0);
    double BCoefficient = (averageAdiabaticIndex - 1.0) / (averageAdiabaticIndex + 1.0);

    double leftCoefficient = sqrt((ACoefficient / leftDensity) / ((BCoefficient * leftPressure) + PVRSStarRegionPressure));
    double rightCoefficient = sqrt((ACoefficient / rightDensity) / ((BCoefficient * rightPressure) + PVRSStarRegionPressure));

    return (0.5 * (leftXVelocity + rightXVelocity)) + (0.5 * ((rightCoefficient * (starRegionPressure - rightPressure)) - (leftCoefficient * (starRegionPressure - leftPressure))));
}

double HLLCSolver::computeTwoShockStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                      EulerMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftPressure = leftStateVector.getPressure();

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightPressure = rightStateVector.getPressure();

    double material1AdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double material2AdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double averageAdiabaticIndex = 0.5 * (material1AdiabaticIndex + material2AdiabaticIndex);

    double PVRSStarRegionPressure = computePVRSStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double ACoefficient = 2.0 / (averageAdiabaticIndex + 1.0);
    double BCoefficient = (averageAdiabaticIndex - 1.0) / (averageAdiabaticIndex + 1.0);

    double leftCoefficient = sqrt((ACoefficient / leftDensity) / ((BCoefficient * leftPressure) + PVRSStarRegionPressure));
    double rightCoefficient = sqrt((ACoefficient / rightDensity) / ((BCoefficient * rightPressure) + PVRSStarRegionPressure));

    return ((leftCoefficient * leftPressure) + (rightCoefficient * rightPressure) - (rightXVelocity - leftXVelocity)) / (leftCoefficient + rightCoefficient);
}

double HLLCSolver::computeTwoShockStarRegionYVelocity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                      EulerMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();
    double topPressure = topStateVector.getPressure();

    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomPressure = bottomStateVector.getPressure();

    double material1AdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double material2AdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double averageAdiabaticIndex = 0.5 * (material1AdiabaticIndex + material2AdiabaticIndex);

    double PVRSStarRegionPressure = computePVRSStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionPressure = computeTwoShockStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double ACoefficient = 2.0 / (averageAdiabaticIndex + 1.0);
    double BCoefficient = (averageAdiabaticIndex - 1.0) / (averageAdiabaticIndex + 1.0);

    double topCoefficient = sqrt((ACoefficient / topDensity) / ((BCoefficient * topPressure) + PVRSStarRegionPressure));
    double bottomCoefficient = sqrt((ACoefficient / bottomDensity) / ((BCoefficient * bottomPressure) + PVRSStarRegionPressure));

    return (0.5 * (topYVelocity + bottomYVelocity)) + (0.5 * ((bottomCoefficient * (starRegionPressure - bottomPressure)) - (topCoefficient * (starRegionPressure - topPressure))));
}

double HLLCSolver::computeTwoShockStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                      EulerMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();
    double topPressure = topStateVector.getPressure();

    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomPressure = bottomStateVector.getPressure();

    double material1AdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double material2AdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double averageAdiabaticIndex = 0.5 * (material1AdiabaticIndex + material2AdiabaticIndex);

    double PVRSStarRegionPressure = computePVRSStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double ACoefficient = 2.0 / (averageAdiabaticIndex + 1.0);
    double BCoefficient = (averageAdiabaticIndex - 1.0) / (averageAdiabaticIndex + 1.0);

    double topCoefficient = sqrt((ACoefficient / topDensity) / ((BCoefficient * topPressure) + PVRSStarRegionPressure));
    double bottomCoefficient = sqrt((ACoefficient / bottomDensity) / ((BCoefficient * bottomPressure) + PVRSStarRegionPressure));

    return ((topCoefficient * topPressure) + (bottomCoefficient * bottomPressure) - (bottomYVelocity - topYVelocity)) / (topCoefficient + bottomCoefficient);
}

double HLLCSolver::computeStarRegionXVelocity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                              EulerMaterialParameters material2Parameters)
{
    double leftPressure = leftStateVector.getPressure();
    double rightPressure = rightStateVector.getPressure();

    double leftAdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double leftStiffeningParameter = material1Parameters.getStiffeningParameter();

    double rightAdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double rightStiffeningParameter = material2Parameters.getStiffeningParameter();

    double minimumPressure = min(leftPressure, rightPressure);
    double maximumPressure = max(leftPressure, rightPressure);
    double maximumPressureRatio = maximumPressure / minimumPressure;

    double PVRSPressure = computePVRSStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (((maximumPressureRatio <= 2.0) && (minimumPressure <= PVRSPressure) && (PVRSPressure <= maximumPressure)) || abs(leftAdiabaticIndex - rightAdiabaticIndex) > pow(10.0, -8.0) ||
            leftStiffeningParameter > 0.0 || rightStiffeningParameter > 0.0)
    {
        return computePVRSStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    }
    else
    {
        if (PVRSPressure < minimumPressure)
        {
            return computeTwoRarefactionStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
        }
        else
        {
            return computeTwoShockStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
        }
    }
}

double HLLCSolver::computeStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                              EulerMaterialParameters material2Parameters)
{
    double leftPressure = leftStateVector.getPressure();
    double rightPressure = rightStateVector.getPressure();

    double leftAdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double leftStiffeningParameter = material1Parameters.getStiffeningParameter();

    double rightAdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double rightStiffeningParameter = material2Parameters.getStiffeningParameter();

    double minimumPressure = min(leftPressure, rightPressure);
    double maximumPressure = max(leftPressure, rightPressure);
    double maximumPressureRatio = maximumPressure / minimumPressure;

    double PVRSPressure = computePVRSStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (((maximumPressureRatio <= 2.0) && (minimumPressure <= PVRSPressure) && (PVRSPressure <= maximumPressure)) || abs(leftAdiabaticIndex - rightAdiabaticIndex) > pow(10.0, -8.0) ||
            leftStiffeningParameter > 0.0 || rightStiffeningParameter > 0.0)
    {
        return PVRSPressure;
    }
    else
    {
        if (PVRSPressure < minimumPressure)
        {
            return computeTwoRarefactionStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
        }
        else
        {
            return computeTwoShockStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
        }
    }
}

double HLLCSolver::computeStarRegionYVelocity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                              EulerMaterialParameters material2Parameters)
{
    double topPressure = topStateVector.getPressure();
    double bottomPressure = bottomStateVector.getPressure();

    double topAdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double topStiffeningParameter = material1Parameters.getStiffeningParameter();

    double bottomAdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double bottomStiffeningParameter = material2Parameters.getStiffeningParameter();

    double minimumPressure = min(topPressure, bottomPressure);
    double maximumPressure = max(topPressure, bottomPressure);
    double maximumPressureRatio = maximumPressure / minimumPressure;

    double PVRSPressure = computePVRSStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (((maximumPressureRatio <= 2.0) && (minimumPressure <= PVRSPressure) && (PVRSPressure <= maximumPressure)) || abs(topAdiabaticIndex - bottomAdiabaticIndex) > pow(10.0, -8.0) ||
            topStiffeningParameter > 0.0 || bottomStiffeningParameter > 0.0)
    {
        return computePVRSStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    }
    else
    {
        if (PVRSPressure < minimumPressure)
        {
            return computeTwoRarefactionStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
        }
        else
        {
            return computeTwoShockStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
        }
    }
}

double HLLCSolver::computeStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                              EulerMaterialParameters material2Parameters)
{
    double topPressure = topStateVector.getPressure();
    double bottomPressure = bottomStateVector.getPressure();

    double topAdiabaticIndex = material1Parameters.getAdiabaticIndex();
    double topStiffeningParameter = material1Parameters.getStiffeningParameter();

    double bottomAdiabaticIndex = material2Parameters.getAdiabaticIndex();
    double bottomStiffeningParameter = material2Parameters.getStiffeningParameter();

    double minimumPressure = min(topPressure, bottomPressure);
    double maximumPressure = max(topPressure, bottomPressure);
    double maximumPressureRatio = maximumPressure / minimumPressure;

    double PVRSPressure = computePVRSStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (((maximumPressureRatio <= 2.0) && (minimumPressure <= PVRSPressure) && (PVRSPressure <= maximumPressure)) || abs(topAdiabaticIndex - bottomAdiabaticIndex) > pow(10.0, -8.0) ||
            topStiffeningParameter > 0.0 || bottomStiffeningParameter > 0.0)
    {
        return PVRSPressure;
    }
    else
    {
        if (PVRSPressure < minimumPressure)
        {
            return computeTwoRarefactionStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
        }
        else
        {
            return computeTwoShockStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
        }
    }
}

double HLLCSolver::computeWaveSpeedWeighting(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters)
{
    double pressure = stateVector.getPressure();

    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    if (starRegionPressure <= pressure)
    {
        return 1.0;
    }
    else
    {
        return sqrt(1.0 + ((adiabaticIndex + 1.0) / (2.0 * adiabaticIndex)) * ((starRegionPressure + stiffeningParameter) / (pressure + stiffeningParameter) - 1.0));
    }
}

double HLLCSolver::computeLeftWaveSpeed(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                        EulerMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftSoundSpeed = leftStateVector.computeSoundSpeed(material1Parameters);

    double starRegionPressure = computeStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    return leftXVelocity - (computeWaveSpeedWeighting(starRegionPressure, leftStateVector, material1Parameters) * leftSoundSpeed);
}

double HLLCSolver::computeRightWaveSpeed(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                         EulerMaterialParameters material2Parameters)
{
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightSoundSpeed = rightStateVector.computeSoundSpeed(material2Parameters);

    double starRegionPressure = computeStarRegionXPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    return rightXVelocity + (computeWaveSpeedWeighting(starRegionPressure, rightStateVector, material2Parameters) * rightSoundSpeed);
}

double HLLCSolver::computeTopWaveSpeed(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                       EulerMaterialParameters material2Parameters)
{
    double topYVelocity = topStateVector.getYVelocity();
    double topSoundSpeed = topStateVector.computeSoundSpeed(material1Parameters);

    double starRegionPressure = computeStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    return topYVelocity - (computeWaveSpeedWeighting(starRegionPressure, topStateVector, material1Parameters) * topSoundSpeed);
}

double HLLCSolver::computeBottomWaveSpeed(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                          EulerMaterialParameters material2Parameters)
{
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomSoundSpeed = bottomStateVector.computeSoundSpeed(material2Parameters);

    double starRegionPressure = computeStarRegionYPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    return bottomYVelocity + (computeWaveSpeedWeighting(starRegionPressure, bottomStateVector, material2Parameters) * bottomSoundSpeed);
}
