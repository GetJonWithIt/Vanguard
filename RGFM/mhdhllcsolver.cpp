#include "mhdhllcsolver.h"

MHDHLLCSolver::MHDHLLCSolver()
{
}

double MHDHLLCSolver::computeStarRegionDensity(double density, double waveSpeed, double velocity, double starRegionWaveSpeed)
{
    return density * ((waveSpeed - velocity) / (waveSpeed - starRegionWaveSpeed));
}

double MHDHLLCSolver::computeLeftStarRegionDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();

    return computeStarRegionDensity(leftDensity, leftWaveSpeed, leftXVelocity, starRegionWaveSpeed);
}

double MHDHLLCSolver::computeRightStarRegionDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                    MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();

    return computeStarRegionDensity(rightDensity, rightWaveSpeed, rightXVelocity, starRegionWaveSpeed);
}

vector<double> MHDHLLCSolver::computeLeftStarRegionConservedVariableVector(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                                           MHDMaterialParameters material2Parameters)
{
    vector<double> starRegionConservedVariableVector(9);

    double starRegionDensity = computeLeftStarRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionTotalEnergy = computeLeftStarRegionTotalEnergy(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionYMomentum = computeLeftStarRegionYMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZMomentum = computeLeftStarRegionZMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeXHLLYMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeXHLLZMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftAuxiliaryField = leftStateVector.getAuxiliaryField();

    starRegionConservedVariableVector[0] = starRegionDensity;

    starRegionConservedVariableVector[1] = starRegionDensity * starRegionWaveSpeed;
    starRegionConservedVariableVector[2] = starRegionYMomentum;
    starRegionConservedVariableVector[3] = starRegionZMomentum;

    starRegionConservedVariableVector[4] = starRegionTotalEnergy;

    starRegionConservedVariableVector[5] = starRegionXMagneticField;
    starRegionConservedVariableVector[6] = starRegionYMagneticField;
    starRegionConservedVariableVector[7] = starRegionZMagneticField;

    starRegionConservedVariableVector[8] = leftAuxiliaryField;

    return starRegionConservedVariableVector;
}

vector<double> MHDHLLCSolver::computeRightStarRegionConservedVariableVector(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                                            MHDMaterialParameters material2Parameters)
{
    vector<double> starRegionConservedVariableVector(9);

    double starRegionDensity = computeRightStarRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionTotalEnergy = computeRightStarRegionTotalEnergy(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionYMomentum = computeRightStarRegionYMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZMomentum = computeRightStarRegionZMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeXHLLYMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeXHLLZMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightAuxiliaryField = rightStateVector.getAuxiliaryField();

    starRegionConservedVariableVector[0] = starRegionDensity;

    starRegionConservedVariableVector[1] = starRegionDensity * starRegionWaveSpeed;
    starRegionConservedVariableVector[2] = starRegionYMomentum;
    starRegionConservedVariableVector[3] = starRegionZMomentum;

    starRegionConservedVariableVector[4] = starRegionTotalEnergy;

    starRegionConservedVariableVector[5] = starRegionXMagneticField;
    starRegionConservedVariableVector[6] = starRegionYMagneticField;
    starRegionConservedVariableVector[7] = starRegionZMagneticField;

    starRegionConservedVariableVector[8] = rightAuxiliaryField;

    return starRegionConservedVariableVector;
}

MHDStateVector MHDHLLCSolver::solveX(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftStateVector;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightStateVector;
    }
    else if (0.0 <= starRegionWaveSpeed)
    {
        vector<double> starRegionConservedVariableVector = computeLeftStarRegionConservedVariableVector(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

        MHDStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material1Parameters);

        return solution;
    }
    else
    {
        vector<double> starRegionConservedVariableVector = computeRightStarRegionConservedVariableVector(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

        MHDStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material2Parameters);

        return solution;
    }
}

double MHDHLLCSolver::computeTotalPressure(MHDStateVector stateVector)
{
    double pressure = stateVector.getPressure();

    double xMagneticField = stateVector.getXMagneticField();
    double yMagneticField = stateVector.getYMagneticField();
    double zMagneticField = stateVector.getZMagneticField();

    double magneticFieldSquared = (xMagneticField * xMagneticField) + (yMagneticField * yMagneticField) + (zMagneticField * zMagneticField);

    return pressure + (0.5 * magneticFieldSquared);
}

double MHDHLLCSolver::computeStarRegionXVelocity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftXMagneticField = leftStateVector.getXMagneticField();

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightXMagneticField = rightStateVector.getXMagneticField();

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftTotalPressure = computeTotalPressure(leftStateVector);
    double rightTotalPressure = computeTotalPressure(rightStateVector);

    double leftCoefficient = leftDensity * (leftWaveSpeed - leftXVelocity);
    double rightCoefficient = rightDensity * (rightWaveSpeed - rightXVelocity);

    double numerator = (rightCoefficient * rightXVelocity) - (leftCoefficient * leftXVelocity) + leftTotalPressure - rightTotalPressure - (leftXMagneticField * leftXMagneticField) +
            (rightXMagneticField * rightXMagneticField);
    double denominator = rightCoefficient - leftCoefficient;

    return numerator / denominator;
}

double MHDHLLCSolver::computeLeftStarRegionTotalPressure(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                         MHDMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftXMagneticField = leftStateVector.getXMagneticField();

    double leftTotalPressure = computeTotalPressure(leftStateVector);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    return (leftDensity * (leftWaveSpeed - leftXVelocity) * (starRegionWaveSpeed - leftXVelocity)) + leftTotalPressure - (leftXMagneticField * leftXMagneticField) +
            (starRegionXMagneticField * starRegionXMagneticField);
}

double MHDHLLCSolver::computeRightStarRegionTotalPressure(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                          MHDMaterialParameters material2Parameters)
{
    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightXMagneticField = rightStateVector.getXMagneticField();

    double rightTotalPressure = computeTotalPressure(rightStateVector);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    return (rightDensity * (rightWaveSpeed - rightXVelocity) * (starRegionWaveSpeed - rightXVelocity)) + rightTotalPressure - (rightXMagneticField * rightXMagneticField) +
            (starRegionXMagneticField * starRegionXMagneticField);
}

double MHDHLLCSolver::computeLeftStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                     MHDMaterialParameters material2Parameters)
{
    double starRegionDensity = computeLeftStarRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeXHLLYMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftYVelocity = leftStateVector.getYVelocity();
    double leftXMagneticField = leftStateVector.getXMagneticField();
    double leftYMagneticField = leftStateVector.getYMagneticField();

    return (starRegionDensity * leftYVelocity) - (((starRegionXMagneticField * starRegionYMagneticField) - (leftXMagneticField * leftYMagneticField)) /
                                                  (leftWaveSpeed - starRegionWaveSpeed));
}

double MHDHLLCSolver::computeRightStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                      MHDMaterialParameters material2Parameters)
{
    double starRegionDensity = computeRightStarRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeXHLLYMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightYVelocity = rightStateVector.getYVelocity();
    double rightXMagneticField = rightStateVector.getXMagneticField();
    double rightYMagneticField = rightStateVector.getYMagneticField();

    return (starRegionDensity * rightYVelocity) - (((starRegionXMagneticField * starRegionYMagneticField) - (rightXMagneticField * rightYMagneticField)) /
                                                   (rightWaveSpeed - starRegionWaveSpeed));
}

double MHDHLLCSolver::computeLeftStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                     MHDMaterialParameters material2Parameters)
{
    double starRegionDensity = computeLeftStarRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeXHLLZMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftZVelocity = leftStateVector.getZVelocity();
    double leftXMagneticField = leftStateVector.getXMagneticField();
    double leftZMagneticField = leftStateVector.getZMagneticField();

    return (starRegionDensity * leftZVelocity) - (((starRegionXMagneticField * starRegionZMagneticField) - (leftXMagneticField * leftZMagneticField)) /
                                                  (leftWaveSpeed - starRegionWaveSpeed));
}

double MHDHLLCSolver::computeRightStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                      MHDMaterialParameters material2Parameters)
{
    double starRegionDensity = computeRightStarRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeXHLLZMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightZVelocity = rightStateVector.getZVelocity();
    double rightXMagneticField = rightStateVector.getXMagneticField();
    double rightZMagneticField = rightStateVector.getZMagneticField();

    return (starRegionDensity * rightZVelocity) - (((starRegionXMagneticField * starRegionZMagneticField) - (rightXMagneticField * rightZMagneticField)) /
                                                   (rightWaveSpeed - starRegionWaveSpeed));
}

double MHDHLLCSolver::computeXStarRegionXMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double leftXMagneticField = leftStateVector.getXMagneticField();
    double rightXMagneticField = rightStateVector.getXMagneticField();

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftFlux = leftStateVector.getAuxiliaryField();
    double rightFlux = rightStateVector.getAuxiliaryField();

    return ((rightWaveSpeed * rightXMagneticField) - (leftWaveSpeed * leftXMagneticField) - (rightFlux - leftFlux)) / (rightWaveSpeed - leftWaveSpeed);
}

double MHDHLLCSolver::computeXStarRegionYMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftYVelocity = leftStateVector.getYVelocity();

    double leftXMagneticField = leftStateVector.getXMagneticField();
    double leftYMagneticField = leftStateVector.getYMagneticField();

    double rightXVelocity = rightStateVector.getXVelocity();
    double rightYVelocity = rightStateVector.getYVelocity();

    double rightXMagneticField = rightStateVector.getXMagneticField();
    double rightYMagneticField = rightStateVector.getYMagneticField();

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftFlux = (leftYMagneticField * leftXVelocity) - (leftXMagneticField * leftYVelocity);
    double rightFlux = (rightYMagneticField * rightXVelocity) - (rightXMagneticField * rightYVelocity);

    return ((rightWaveSpeed * rightYMagneticField) - (leftWaveSpeed * leftYMagneticField) - (rightFlux - leftFlux)) / (rightWaveSpeed - leftWaveSpeed);
}

double MHDHLLCSolver::computeXStarRegionZMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftZVelocity = leftStateVector.getZVelocity();

    double leftXMagneticField = leftStateVector.getXMagneticField();
    double leftZMagneticField = leftStateVector.getZMagneticField();

    double rightXVelocity = rightStateVector.getXVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    double rightXMagneticField = rightStateVector.getXMagneticField();
    double rightZMagneticField = rightStateVector.getZMagneticField();

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftFlux = (leftZMagneticField * leftXVelocity) - (leftXMagneticField * leftZVelocity);
    double rightFlux = (rightZMagneticField * rightXVelocity) - (rightXMagneticField * rightZVelocity);

    return ((rightWaveSpeed * rightZMagneticField) - (leftWaveSpeed * leftZMagneticField) - (rightFlux - leftFlux)) / (rightWaveSpeed - leftWaveSpeed);
}

double MHDHLLCSolver::computeStarRegionXMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftXVelocity = leftStateVector.getXVelocity();
    double rightXVelocity = rightStateVector.getXVelocity();

    double leftPressure = leftStateVector.getPressure();
    double rightPressure = rightStateVector.getPressure();

    double leftXMagneticField = leftStateVector.getXMagneticField();
    double leftYMagneticField = leftStateVector.getYMagneticField();
    double leftZMagneticField = leftStateVector.getZMagneticField();

    double rightXMagneticField = rightStateVector.getXMagneticField();
    double rightYMagneticField = rightStateVector.getYMagneticField();
    double rightZMagneticField = rightStateVector.getZMagneticField();

    double leftMagneticFieldSquared = (leftXMagneticField * leftXMagneticField) + (leftYMagneticField * leftYMagneticField) + (leftZMagneticField * leftZMagneticField);
    double rightMagneticFieldSquared = (rightXMagneticField * rightXMagneticField) + (rightYMagneticField * rightYMagneticField) + (rightZMagneticField * rightZMagneticField);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftFlux = (leftDensity * (leftXVelocity * leftXVelocity)) + leftPressure + (0.5 * leftMagneticFieldSquared) - (leftXMagneticField * leftXMagneticField);
    double rightFlux = (rightDensity * (rightXVelocity * rightXVelocity)) + rightPressure + (0.5 * rightMagneticFieldSquared) - (rightXMagneticField * rightXMagneticField);

    return ((rightWaveSpeed * rightDensity * rightXVelocity) - (leftWaveSpeed * leftDensity * leftXVelocity) - (rightFlux - leftFlux)) / (rightWaveSpeed - leftWaveSpeed);
}

double MHDHLLCSolver::computeStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftXVelocity = leftStateVector.getXVelocity();
    double leftYVelocity = leftStateVector.getYVelocity();

    double rightXVelocity = rightStateVector.getXVelocity();
    double rightYVelocity = rightStateVector.getYVelocity();

    double leftXMagneticField = leftStateVector.getXMagneticField();
    double leftYMagneticField = leftStateVector.getYMagneticField();

    double rightXMagneticField = rightStateVector.getXMagneticField();
    double rightYMagneticField = rightStateVector.getYMagneticField();

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftFlux = (leftDensity * (leftXVelocity * leftYVelocity)) - (leftXMagneticField * leftYMagneticField);
    double rightFlux = (rightDensity * (rightXVelocity * rightYVelocity)) - (rightXMagneticField * rightYMagneticField);

    return ((rightWaveSpeed * rightDensity * rightYVelocity) - (leftWaveSpeed * leftDensity * leftYVelocity) - (rightFlux - leftFlux)) / (rightWaveSpeed - leftWaveSpeed);
}

double MHDHLLCSolver::computeStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftXVelocity = leftStateVector.getXVelocity();
    double leftZVelocity = leftStateVector.getZVelocity();

    double rightXVelocity = rightStateVector.getXVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    double leftXMagneticField = leftStateVector.getXMagneticField();
    double leftZMagneticField = leftStateVector.getZMagneticField();

    double rightXMagneticField = rightStateVector.getXMagneticField();
    double rightZMagneticField = rightStateVector.getZMagneticField();

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftFlux = (leftDensity * (leftXVelocity * leftZVelocity)) - (leftXMagneticField * leftZMagneticField);
    double rightFlux = (rightDensity * (rightXVelocity * rightZVelocity)) - (rightXMagneticField * rightZMagneticField);

    return ((rightWaveSpeed * rightDensity * rightZVelocity) - (leftWaveSpeed * leftDensity * leftZVelocity) - (rightFlux - leftFlux)) / (rightWaveSpeed - leftWaveSpeed);
}

double MHDHLLCSolver::computeXStarRegionDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftXVelocity = leftStateVector.getXVelocity();
    double rightXVelocity = rightStateVector.getXVelocity();

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftFlux = leftDensity * leftXVelocity;
    double rightFlux = rightDensity * rightXVelocity;

    return ((rightWaveSpeed * rightDensity) - (leftWaveSpeed * leftDensity) - (rightFlux - leftFlux)) / (rightWaveSpeed - leftWaveSpeed);
}

double MHDHLLCSolver::computeXHLLDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                         MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double starRegionDensity = computeXStarRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftDensity;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightDensity;
    }
    else
    {
        return starRegionDensity;
    }
}

double MHDHLLCSolver::computeXHLLXMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftXMagneticField = leftStateVector.getXMagneticField();
    double rightXMagneticField = rightStateVector.getXMagneticField();

    double starRegionXMagneticField = computeXStarRegionXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftXMagneticField;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightXMagneticField;
    }
    else
    {
        return starRegionXMagneticField;
    }
}

double MHDHLLCSolver::computeXHLLYMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftYMagneticField = leftStateVector.getYMagneticField();
    double rightYMagneticField = rightStateVector.getYMagneticField();

    double starRegionYMagneticField = computeXStarRegionYMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftYMagneticField;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightYMagneticField;
    }
    else
    {
        return starRegionYMagneticField;
    }
}

double MHDHLLCSolver::computeXHLLZMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftZMagneticField = leftStateVector.getZMagneticField();
    double rightZMagneticField = rightStateVector.getZMagneticField();

    double starRegionZMagneticField = computeXStarRegionZMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftZMagneticField;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightZMagneticField;
    }
    else
    {
        return starRegionZMagneticField;
    }
}

double MHDHLLCSolver::computeHLLXMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftXVelocity = leftStateVector.getXVelocity();
    double rightXVelocity = rightStateVector.getXVelocity();

    double starRegionXMomentum = computeStarRegionXMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftDensity * leftXVelocity;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightDensity * rightXVelocity;
    }
    else
    {
        return starRegionXMomentum;
    }
}

double MHDHLLCSolver::computeHLLYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftYVelocity = leftStateVector.getYVelocity();
    double rightYVelocity = rightStateVector.getYVelocity();

    double starRegionYMomentum = computeStarRegionYMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftDensity * leftYVelocity;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightDensity * rightYVelocity;
    }
    else
    {
        return starRegionYMomentum;
    }
}

double MHDHLLCSolver::computeHLLZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftZVelocity = leftStateVector.getZVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    double starRegionZMomentum = computeStarRegionZMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftDensity * leftZVelocity;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightDensity * rightZVelocity;
    }
    else
    {
        return starRegionZMomentum;
    }
}

double MHDHLLCSolver::computeLeftStarRegionTotalEnergy(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftXVelocity = leftStateVector.getXVelocity();
    double leftYVelocity = leftStateVector.getYVelocity();
    double leftZVelocity = leftStateVector.getZVelocity();

    double leftXMagneticField = leftStateVector.getXMagneticField();
    double leftYMagneticField = leftStateVector.getYMagneticField();
    double leftZMagneticField = leftStateVector.getZMagneticField();

    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeXHLLYMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeXHLLZMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionDensity = computeXHLLDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionXVelocity = computeHLLXMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionYVelocity = computeHLLYMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionZVelocity = computeHLLZMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionMagneticFieldVelocityVectorProduct = (starRegionXMagneticField * starRegionXVelocity) + (starRegionYMagneticField * starRegionYVelocity) +
            (starRegionZMagneticField * starRegionZVelocity);

    double magneticFieldVelocityVectorProduct = (leftXMagneticField * leftXVelocity) + (leftYMagneticField * leftYVelocity) + (leftZMagneticField * leftZVelocity);

    double leftTotalPressure = computeTotalPressure(leftStateVector);
    double starRegionTotalPressure = computeLeftStarRegionTotalPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftTotalEnergy = leftStateVector.computeTotalEnergy(material1Parameters);

    double firstTerm = leftTotalEnergy * ((leftWaveSpeed - leftXVelocity) / (leftWaveSpeed - starRegionWaveSpeed));
    double secondTerm = ((starRegionTotalPressure * starRegionWaveSpeed) - (leftTotalPressure * leftXVelocity) -
                         ((starRegionXMagneticField * starRegionMagneticFieldVelocityVectorProduct) - (leftXMagneticField * magneticFieldVelocityVectorProduct))) /
            (leftWaveSpeed - starRegionWaveSpeed);

    return firstTerm + secondTerm;
}

double MHDHLLCSolver::computeRightStarRegionTotalEnergy(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                        MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightXVelocity = rightStateVector.getXVelocity();
    double rightYVelocity = rightStateVector.getYVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    double rightXMagneticField = rightStateVector.getXMagneticField();
    double rightYMagneticField = rightStateVector.getYMagneticField();
    double rightZMagneticField = rightStateVector.getZMagneticField();

    double starRegionXMagneticField = computeXHLLXMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeXHLLYMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeXHLLZMagneticField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionDensity = computeXHLLDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionXVelocity = computeHLLXMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionYVelocity = computeHLLYMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionZVelocity = computeHLLZMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionMagneticFieldVelocityVectorProduct = (starRegionXMagneticField * starRegionXVelocity) + (starRegionYMagneticField * starRegionYVelocity) +
            (starRegionZMagneticField * starRegionZVelocity);

    double magneticFieldVelocityVectorProduct = (rightXMagneticField * rightXVelocity) + (rightYMagneticField * rightYVelocity) + (rightZMagneticField * rightZVelocity);

    double rightTotalPressure = computeTotalPressure(rightStateVector);
    double starRegionTotalPressure = computeRightStarRegionTotalPressure(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightTotalEnergy = rightStateVector.computeTotalEnergy(material2Parameters);

    double firstTerm = rightTotalEnergy * ((rightWaveSpeed - rightXVelocity) / (rightWaveSpeed - starRegionWaveSpeed));
    double secondTerm = ((starRegionTotalPressure * starRegionWaveSpeed) - (rightTotalPressure * rightXVelocity) -
                         ((starRegionXMagneticField * starRegionMagneticFieldVelocityVectorProduct) - (rightXMagneticField * magneticFieldVelocityVectorProduct))) /
            (rightWaveSpeed - starRegionWaveSpeed);

    return firstTerm + secondTerm;
}

double MHDHLLCSolver::computeLeftWaveSpeed(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double rightXVelocity = rightStateVector.getXVelocity();

    double leftXFastMagnetoAcousticSpeed = leftStateVector.computeXFastMagnetoAcousticSpeed(material1Parameters);
    double rightXFastMagnetoAcousticSpeed = rightStateVector.computeXFastMagnetoAcousticSpeed(material2Parameters);

    return min(leftXVelocity, rightXVelocity) - max(leftXFastMagnetoAcousticSpeed, rightXFastMagnetoAcousticSpeed);
}

double MHDHLLCSolver::computeRightWaveSpeed(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double rightXVelocity = rightStateVector.getXVelocity();

    double leftXFastMagnetoAcousticSpeed = leftStateVector.computeXFastMagnetoAcousticSpeed(material1Parameters);
    double rightXFastMagnetoAcousticSpeed = rightStateVector.computeXFastMagnetoAcousticSpeed(material2Parameters);

    return max(leftXVelocity, rightXVelocity) + max(leftXFastMagnetoAcousticSpeed, rightXFastMagnetoAcousticSpeed);
}
