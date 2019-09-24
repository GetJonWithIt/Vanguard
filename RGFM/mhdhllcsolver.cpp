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

double MHDHLLCSolver::computeTopStarRegionDensity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();

    return computeStarRegionDensity(topDensity, topWaveSpeed, topYVelocity, starRegionWaveSpeed);
}

double MHDHLLCSolver::computeBottomStarRegionDensity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                     MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    return computeStarRegionDensity(bottomDensity, bottomWaveSpeed, bottomYVelocity, starRegionWaveSpeed);
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

    double starRegionAuxiliaryField = computeXHLLAuxiliaryField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    starRegionConservedVariableVector[0] = starRegionDensity;

    starRegionConservedVariableVector[1] = starRegionDensity * starRegionWaveSpeed;
    starRegionConservedVariableVector[2] = starRegionYMomentum;
    starRegionConservedVariableVector[3] = starRegionZMomentum;

    starRegionConservedVariableVector[4] = starRegionTotalEnergy;

    starRegionConservedVariableVector[5] = starRegionXMagneticField;
    starRegionConservedVariableVector[6] = starRegionYMagneticField;
    starRegionConservedVariableVector[7] = starRegionZMagneticField;

    starRegionConservedVariableVector[8] = starRegionAuxiliaryField;

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

    double starRegionAuxiliaryField = computeXHLLAuxiliaryField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    starRegionConservedVariableVector[0] = starRegionDensity;

    starRegionConservedVariableVector[1] = starRegionDensity * starRegionWaveSpeed;
    starRegionConservedVariableVector[2] = starRegionYMomentum;
    starRegionConservedVariableVector[3] = starRegionZMomentum;

    starRegionConservedVariableVector[4] = starRegionTotalEnergy;

    starRegionConservedVariableVector[5] = starRegionXMagneticField;
    starRegionConservedVariableVector[6] = starRegionYMagneticField;
    starRegionConservedVariableVector[7] = starRegionZMagneticField;

    starRegionConservedVariableVector[8] = starRegionAuxiliaryField;

    return starRegionConservedVariableVector;
}

vector<double> MHDHLLCSolver::computeTopStarRegionConservedVariableVector(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                                          MHDMaterialParameters material2Parameters)
{
    vector<double> starRegionConservedVariableVector(9);

    double starRegionDensity = computeTopStarRegionDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionTotalEnergy = computeTopStarRegionTotalEnergy(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionXMomentum = computeTopStarRegionXMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionZMomentum = computeTopStarRegionZMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionXMagneticField = computeYHLLXMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeYHLLZMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionAuxiliaryField = computeYHLLAuxiliaryField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    starRegionConservedVariableVector[0] = starRegionDensity;

    starRegionConservedVariableVector[1] = starRegionXMomentum;
    starRegionConservedVariableVector[2] = starRegionDensity * starRegionWaveSpeed;
    starRegionConservedVariableVector[3] = starRegionZMomentum;

    starRegionConservedVariableVector[4] = starRegionTotalEnergy;

    starRegionConservedVariableVector[5] = starRegionXMagneticField;
    starRegionConservedVariableVector[6] = starRegionYMagneticField;
    starRegionConservedVariableVector[7] = starRegionZMagneticField;

    starRegionConservedVariableVector[8] = starRegionAuxiliaryField;

    return starRegionConservedVariableVector;
}

vector<double> MHDHLLCSolver::computeBottomStarRegionConservedVariableVector(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                                             MHDMaterialParameters material2Parameters)
{
    vector<double> starRegionConservedVariableVector(9);

    double starRegionDensity = computeBottomStarRegionDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionTotalEnergy = computeBottomStarRegionTotalEnergy(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionXMomentum = computeBottomStarRegionXMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionZMomentum = computeBottomStarRegionZMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionXMagneticField = computeYHLLXMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeYHLLZMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionAuxiliaryField = computeYHLLAuxiliaryField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    starRegionConservedVariableVector[0] = starRegionDensity;

    starRegionConservedVariableVector[1] = starRegionXMomentum;
    starRegionConservedVariableVector[2] = starRegionDensity * starRegionWaveSpeed;
    starRegionConservedVariableVector[3] = starRegionZMomentum;

    starRegionConservedVariableVector[4] = starRegionTotalEnergy;

    starRegionConservedVariableVector[5] = starRegionXMagneticField;
    starRegionConservedVariableVector[6] = starRegionYMagneticField;
    starRegionConservedVariableVector[7] = starRegionZMagneticField;

    starRegionConservedVariableVector[8] = starRegionAuxiliaryField;

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

MHDStateVector MHDHLLCSolver::solveY(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topStateVector;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomStateVector;
    }
    else if (0.0 <= starRegionWaveSpeed)
    {
        vector<double> starRegionConservedVariableVector = computeTopStarRegionConservedVariableVector(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

        MHDStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material1Parameters);

        return solution;
    }
    else
    {
        vector<double> starRegionConservedVariableVector = computeBottomStarRegionConservedVariableVector(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

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

double MHDHLLCSolver::computeStarRegionYVelocity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();
    double topYMagneticField = topStateVector.getYMagneticField();

    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topTotalPressure = computeTotalPressure(topStateVector);
    double bottomTotalPressure = computeTotalPressure(bottomStateVector);

    double topCoefficient = topDensity * (topWaveSpeed - topYVelocity);
    double bottomCoefficient = bottomDensity * (bottomWaveSpeed - bottomYVelocity);

    double numerator = (bottomCoefficient * bottomYVelocity) - (topCoefficient * topYVelocity) + topTotalPressure - bottomTotalPressure - (topYMagneticField * topYMagneticField) +
            (bottomYMagneticField * bottomYMagneticField);
    double denominator = bottomCoefficient - topCoefficient;

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

double MHDHLLCSolver::computeTopStarRegionTotalPressure(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                        MHDMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double topYVelocity = topStateVector.getYVelocity();
    double topYMagneticField = topStateVector.getYMagneticField();

    double topTotalPressure = computeTotalPressure(topStateVector);
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    return (topDensity * (topWaveSpeed - topYVelocity) * (starRegionWaveSpeed - topYVelocity)) + topTotalPressure - (topYMagneticField * topYMagneticField) +
            (starRegionYMagneticField * starRegionYMagneticField);
}

double MHDHLLCSolver::computeBottomStarRegionTotalPressure(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                           MHDMaterialParameters material2Parameters)
{
    double bottomDensity = bottomStateVector.getDensity();
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();

    double bottomTotalPressure = computeTotalPressure(bottomStateVector);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    return (bottomDensity * (bottomWaveSpeed - bottomYVelocity) * (starRegionWaveSpeed - bottomYVelocity)) + bottomTotalPressure - (bottomYMagneticField * bottomYMagneticField) +
            (starRegionYMagneticField * starRegionYMagneticField);
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

double MHDHLLCSolver::computeTopStarRegionXMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                    MHDMaterialParameters material2Parameters)
{
    double starRegionDensity = computeTopStarRegionDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionXMagneticField = computeYHLLXMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topXVelocity = topStateVector.getXVelocity();
    double topXMagneticField = topStateVector.getXMagneticField();
    double topYMagneticField = topStateVector.getYMagneticField();

    return (starRegionDensity * topXVelocity) - (((starRegionYMagneticField * starRegionXMagneticField) - (topYMagneticField * topXMagneticField)) /
                                                 (topWaveSpeed - starRegionWaveSpeed));
}

double MHDHLLCSolver::computeBottomStarRegionXMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double starRegionDensity = computeBottomStarRegionDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionXMagneticField = computeYHLLXMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double bottomXVelocity = bottomStateVector.getXVelocity();
    double bottomXMagneticField = bottomStateVector.getXMagneticField();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();

    return (starRegionDensity * bottomXVelocity) - (((starRegionYMagneticField * starRegionXMagneticField) - (bottomYMagneticField * bottomXMagneticField)) /
                                                    (bottomWaveSpeed - starRegionWaveSpeed));
}

double MHDHLLCSolver::computeTopStarRegionZMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                    MHDMaterialParameters material2Parameters)
{
    double starRegionDensity = computeTopStarRegionDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeYHLLZMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topZVelocity = topStateVector.getZVelocity();
    double topYMagneticField = topStateVector.getYMagneticField();
    double topZMagneticField = topStateVector.getZMagneticField();

    return (starRegionDensity * topZVelocity) - (((starRegionYMagneticField * starRegionZMagneticField) - (topYMagneticField * topZMagneticField)) /
                                                 (topWaveSpeed - starRegionWaveSpeed));
}

double MHDHLLCSolver::computeBottomStarRegionZMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double starRegionDensity = computeBottomStarRegionDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeYHLLZMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double bottomZVelocity = bottomStateVector.getZVelocity();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();
    double bottomZMagneticField = bottomStateVector.getZMagneticField();

    return (starRegionDensity * bottomZVelocity) - (((starRegionYMagneticField * starRegionZMagneticField) - (bottomYMagneticField * bottomZMagneticField)) /
                                                    (bottomWaveSpeed - starRegionWaveSpeed));
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

double MHDHLLCSolver::computeYStarRegionXMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double topXVelocity = topStateVector.getXVelocity();
    double topYVelocity = topStateVector.getYVelocity();

    double topXMagneticField = topStateVector.getXMagneticField();
    double topYMagneticField = topStateVector.getYMagneticField();

    double bottomXVelocity = bottomStateVector.getXVelocity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    double bottomXMagneticField = bottomStateVector.getXMagneticField();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topFlux = (topXMagneticField * topYVelocity) - (topYMagneticField * topXVelocity);
    double bottomFlux = (bottomXMagneticField * bottomYVelocity) - (bottomYMagneticField * bottomXVelocity);

    return ((bottomWaveSpeed * bottomXMagneticField) - (topWaveSpeed * topXMagneticField) - (bottomFlux - topFlux)) / (bottomWaveSpeed - topWaveSpeed);
}

double MHDHLLCSolver::computeYStarRegionYMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double topYMagneticField = topStateVector.getYMagneticField();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topFlux = topStateVector.getAuxiliaryField();
    double bottomFlux = bottomStateVector.getAuxiliaryField();

    return ((bottomWaveSpeed * bottomYMagneticField) - (topWaveSpeed * topYMagneticField) - (bottomFlux - topFlux)) / (bottomWaveSpeed - topWaveSpeed);
}

double MHDHLLCSolver::computeYStarRegionZMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double topYVelocity = topStateVector.getYVelocity();
    double topZVelocity = topStateVector.getZVelocity();

    double topYMagneticField = topStateVector.getYMagneticField();
    double topZMagneticField = topStateVector.getZMagneticField();

    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomZVelocity = bottomStateVector.getZVelocity();

    double bottomYMagneticField = bottomStateVector.getYMagneticField();
    double bottomZMagneticField = bottomStateVector.getZMagneticField();

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topFlux = (topZMagneticField * topYVelocity) - (topYMagneticField * topZVelocity);
    double bottomFlux = (bottomZMagneticField * bottomYVelocity) - (bottomYMagneticField * bottomZVelocity);

    return ((bottomWaveSpeed * bottomZMagneticField) - (topWaveSpeed * topZMagneticField) - (bottomFlux - topFlux)) / (bottomWaveSpeed - topWaveSpeed);
}

double MHDHLLCSolver::computeXStarRegionAuxiliaryField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double leftAuxiliaryField = leftStateVector.getAuxiliaryField();
    double rightAuxiliaryField = rightStateVector.getAuxiliaryField();

    double leftXMagneticField = leftStateVector.getXMagneticField();
    double rightXMagneticField = rightStateVector.getXMagneticField();

    double leftHyperbolicWaveSpeedSquared = material1Parameters.computeHyperbolicWaveSpeedSquared();
    double rightHyperbolicWaveSpeedSquared = material2Parameters.computeHyperbolicWaveSpeedSquared();

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftFlux = leftHyperbolicWaveSpeedSquared * leftXMagneticField;
    double rightFlux = rightHyperbolicWaveSpeedSquared * rightXMagneticField;

    return ((rightWaveSpeed * rightAuxiliaryField) - (leftWaveSpeed * leftAuxiliaryField) - (rightFlux - leftFlux)) / (rightWaveSpeed - leftWaveSpeed);
}

double MHDHLLCSolver::computeYStarRegionAuxiliaryField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters)
{
    double topAuxiliaryField = topStateVector.getAuxiliaryField();
    double bottomAuxiliaryField = bottomStateVector.getAuxiliaryField();

    double topYMagneticField = topStateVector.getYMagneticField();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();

    double topHyperbolicWaveSpeedSquared = material1Parameters.computeHyperbolicWaveSpeedSquared();
    double bottomHyperbolicWaveSpeedSquared = material2Parameters.computeHyperbolicWaveSpeedSquared();

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topFlux = topHyperbolicWaveSpeedSquared * topYMagneticField;
    double bottomFlux = bottomHyperbolicWaveSpeedSquared * bottomYMagneticField;

    return ((bottomWaveSpeed * bottomAuxiliaryField) - (topWaveSpeed * topAuxiliaryField) - (bottomFlux - topFlux)) / (bottomWaveSpeed - topWaveSpeed);
}

double MHDHLLCSolver::computeXStarRegionXMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
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

double MHDHLLCSolver::computeXStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
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

double MHDHLLCSolver::computeXStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
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

double MHDHLLCSolver::computeYStarRegionXMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double bottomDensity = bottomStateVector.getDensity();

    double topXVelocity = topStateVector.getXVelocity();
    double topYVelocity = topStateVector.getYVelocity();

    double bottomXVelocity = bottomStateVector.getXVelocity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    double topXMagneticField = topStateVector.getXMagneticField();
    double topYMagneticField = topStateVector.getYMagneticField();

    double bottomXMagneticField = bottomStateVector.getXMagneticField();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topFlux = (topDensity * (topYVelocity * topXVelocity)) - (topYMagneticField * topXMagneticField);
    double bottomFlux = (bottomDensity * (bottomYVelocity * bottomXVelocity)) - (bottomYMagneticField * bottomXMagneticField);

    return ((bottomWaveSpeed * bottomDensity * bottomXVelocity) - (topWaveSpeed * topDensity * topXVelocity) - (bottomFlux - topFlux)) / (bottomWaveSpeed - topWaveSpeed);
}

double MHDHLLCSolver::computeYStarRegionYMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double bottomDensity = bottomStateVector.getDensity();

    double topYVelocity = topStateVector.getYVelocity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    double topPressure = topStateVector.getPressure();
    double bottomPressure = bottomStateVector.getPressure();

    double topXMagneticField = topStateVector.getXMagneticField();
    double topYMagneticField = topStateVector.getYMagneticField();
    double topZMagneticField = topStateVector.getZMagneticField();

    double bottomXMagneticField = bottomStateVector.getXMagneticField();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();
    double bottomZMagneticField = bottomStateVector.getZMagneticField();

    double topMagneticFieldSquared = (topXMagneticField * topXMagneticField) + (topYMagneticField * topYMagneticField) + (topZMagneticField * topZMagneticField);
    double bottomMagneticFieldSquared = (bottomXMagneticField * bottomXMagneticField) + (bottomYMagneticField * bottomYMagneticField) + (bottomZMagneticField * bottomZMagneticField);

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topFlux = (topDensity * (topYVelocity * topYVelocity)) + topPressure + (0.5 * topMagneticFieldSquared) - (topYMagneticField * topYMagneticField);
    double bottomFlux = (bottomDensity * (bottomYVelocity * bottomYVelocity)) + bottomPressure + (0.5 * bottomMagneticFieldSquared) - (bottomYMagneticField * bottomYMagneticField);

    return ((bottomWaveSpeed * bottomDensity * bottomYVelocity) - (topWaveSpeed * topDensity * topYVelocity) - (bottomFlux - topFlux)) / (bottomWaveSpeed - topWaveSpeed);
}

double MHDHLLCSolver::computeYStarRegionZMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double bottomDensity = bottomStateVector.getDensity();

    double topYVelocity = topStateVector.getYVelocity();
    double topZVelocity = topStateVector.getZVelocity();

    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomZVelocity = bottomStateVector.getZVelocity();

    double topYMagneticField = topStateVector.getYMagneticField();
    double topZMagneticField = topStateVector.getZMagneticField();

    double bottomYMagneticField = bottomStateVector.getYMagneticField();
    double bottomZMagneticField = bottomStateVector.getZMagneticField();

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topFlux = (topDensity * (topYVelocity * topZVelocity)) - (topYMagneticField * topZMagneticField);
    double bottomFlux = (bottomDensity * (bottomYVelocity * bottomZVelocity)) - (bottomYMagneticField * bottomZMagneticField);

    return ((bottomWaveSpeed * bottomDensity * bottomZVelocity) - (topWaveSpeed * topDensity * topZVelocity) - (bottomFlux - topFlux)) / (bottomWaveSpeed - topWaveSpeed);
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

double MHDHLLCSolver::computeYStarRegionDensity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double topDensity = topStateVector.getDensity();
    double bottomDensity = bottomStateVector.getDensity();

    double topYVelocity = topStateVector.getYVelocity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topFlux = topDensity * topYVelocity;
    double bottomFlux = bottomDensity * bottomYVelocity;

    return ((bottomWaveSpeed * bottomDensity) - (topWaveSpeed * topDensity) - (bottomFlux - topFlux)) / (bottomWaveSpeed - topWaveSpeed);
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

double MHDHLLCSolver::computeYHLLDensity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                         MHDMaterialParameters material2Parameters)
{
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topDensity = topStateVector.getDensity();
    double bottomDensity = bottomStateVector.getDensity();

    double starRegionDensity = computeYStarRegionDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topDensity;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomDensity;
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

double MHDHLLCSolver::computeYHLLXMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topXMagneticField = topStateVector.getXMagneticField();
    double bottomXMagneticField = bottomStateVector.getXMagneticField();

    double starRegionXMagneticField = computeYStarRegionXMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topXMagneticField;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomXMagneticField;
    }
    else
    {
        return starRegionXMagneticField;
    }
}

double MHDHLLCSolver::computeYHLLYMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topYMagneticField = topStateVector.getYMagneticField();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();

    double starRegionYMagneticField = computeYStarRegionYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topYMagneticField;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomYMagneticField;
    }
    else
    {
        return starRegionYMagneticField;
    }
}

double MHDHLLCSolver::computeYHLLZMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topZMagneticField = topStateVector.getZMagneticField();
    double bottomZMagneticField = bottomStateVector.getZMagneticField();

    double starRegionZMagneticField = computeYStarRegionZMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topZMagneticField;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomZMagneticField;
    }
    else
    {
        return starRegionZMagneticField;
    }
}

double MHDHLLCSolver::computeXHLLAuxiliaryField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftAuxiliaryField = leftStateVector.getAuxiliaryField();
    double rightAuxiliaryField = rightStateVector.getAuxiliaryField();

    double starRegionAuxiliaryField = computeXStarRegionAuxiliaryField(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    if (0.0 <= leftWaveSpeed)
    {
        return leftAuxiliaryField;
    }
    else if (rightWaveSpeed <= 0.0)
    {
        return rightAuxiliaryField;
    }
    else
    {
        return starRegionAuxiliaryField;
    }
}

double MHDHLLCSolver::computeYHLLAuxiliaryField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters)
{
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topAuxiliaryField = topStateVector.getAuxiliaryField();
    double bottomAuxiliaryField = bottomStateVector.getAuxiliaryField();

    double starRegionAuxiliaryField = computeYStarRegionAuxiliaryField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topAuxiliaryField;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomAuxiliaryField;
    }
    else
    {
        return starRegionAuxiliaryField;
    }
}

double MHDHLLCSolver::computeXHLLXMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftXVelocity = leftStateVector.getXVelocity();
    double rightXVelocity = rightStateVector.getXVelocity();

    double starRegionXMomentum = computeXStarRegionXMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

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

double MHDHLLCSolver::computeXHLLYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftYVelocity = leftStateVector.getYVelocity();
    double rightYVelocity = rightStateVector.getYVelocity();

    double starRegionYMomentum = computeXStarRegionYMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

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

double MHDHLLCSolver::computeXHLLZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double rightDensity = rightStateVector.getDensity();

    double leftZVelocity = leftStateVector.getZVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    double starRegionZMomentum = computeXStarRegionZMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

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

double MHDHLLCSolver::computeYHLLXMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topDensity = topStateVector.getDensity();
    double bottomDensity = bottomStateVector.getDensity();

    double topXVelocity = topStateVector.getXVelocity();
    double bottomXVelocity = bottomStateVector.getXVelocity();

    double starRegionXMomentum = computeYStarRegionXMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topDensity * topXVelocity;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomDensity * bottomXVelocity;
    }
    else
    {
        return starRegionXMomentum;
    }
}

double MHDHLLCSolver::computeYHLLYMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topDensity = topStateVector.getDensity();
    double bottomDensity = bottomStateVector.getDensity();

    double topYVelocity = topStateVector.getYVelocity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    double starRegionYMomentum = computeYStarRegionYMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topDensity * topYVelocity;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomDensity * bottomYVelocity;
    }
    else
    {
        return starRegionYMomentum;
    }
}

double MHDHLLCSolver::computeYHLLZMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topDensity = topStateVector.getDensity();
    double bottomDensity = bottomStateVector.getDensity();

    double topZVelocity = topStateVector.getZVelocity();
    double bottomZVelocity = bottomStateVector.getZVelocity();

    double starRegionZMomentum = computeYStarRegionZMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    if (0.0 <= topWaveSpeed)
    {
        return topDensity * topZVelocity;
    }
    else if (bottomWaveSpeed <= 0.0)
    {
        return bottomDensity * bottomZVelocity;
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
    double starRegionXVelocity = computeXHLLXMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionYVelocity = computeXHLLYMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionZVelocity = computeXHLLZMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
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
    double starRegionXVelocity = computeXHLLXMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionYVelocity = computeXHLLYMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionZVelocity = computeXHLLZMomentum(leftStateVector, rightStateVector, material1Parameters, material2Parameters) / starRegionDensity;
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

double MHDHLLCSolver::computeTopStarRegionTotalEnergy(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                      MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double topWaveSpeed = computeTopWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topXVelocity = topStateVector.getXVelocity();
    double topYVelocity = topStateVector.getYVelocity();
    double topZVelocity = topStateVector.getZVelocity();

    double topXMagneticField = topStateVector.getXMagneticField();
    double topYMagneticField = topStateVector.getYMagneticField();
    double topZMagneticField = topStateVector.getZMagneticField();

    double starRegionXMagneticField = computeYHLLXMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeYHLLZMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionDensity = computeYHLLDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionXVelocity = computeYHLLXMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionYVelocity = computeYHLLYMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionZVelocity = computeYHLLZMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionMagneticFieldVelocityVectorProduct = (starRegionXMagneticField * starRegionXVelocity) + (starRegionYMagneticField * starRegionYVelocity) +
            (starRegionZMagneticField * starRegionZVelocity);

    double magneticFieldVelocityVectorProduct = (topXMagneticField * topXVelocity) + (topYMagneticField * topYVelocity) + (topZMagneticField * topZVelocity);

    double topTotalPressure = computeTotalPressure(topStateVector);
    double starRegionTotalPressure = computeTopStarRegionTotalPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double topTotalEnergy = topStateVector.computeTotalEnergy(material1Parameters);

    double firstTerm = topTotalEnergy * ((topWaveSpeed - topYVelocity) / (topWaveSpeed - starRegionWaveSpeed));
    double secondTerm = ((starRegionTotalPressure * starRegionWaveSpeed) - (topTotalPressure * topYVelocity) -
                         ((starRegionYMagneticField * starRegionMagneticFieldVelocityVectorProduct) - (topYMagneticField * magneticFieldVelocityVectorProduct))) /
            (topWaveSpeed - starRegionWaveSpeed);

    return firstTerm + secondTerm;
}

double MHDHLLCSolver::computeBottomStarRegionTotalEnergy(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                         MHDMaterialParameters material2Parameters)
{
    double starRegionWaveSpeed = computeStarRegionYVelocity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double bottomWaveSpeed = computeBottomWaveSpeed(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double bottomXVelocity = bottomStateVector.getXVelocity();
    double bottomYVelocity = bottomStateVector.getYVelocity();
    double bottomZVelocity = bottomStateVector.getZVelocity();

    double bottomXMagneticField = bottomStateVector.getXMagneticField();
    double bottomYMagneticField = bottomStateVector.getYMagneticField();
    double bottomZMagneticField = bottomStateVector.getZMagneticField();

    double starRegionXMagneticField = computeYHLLXMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionYMagneticField = computeYHLLYMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionZMagneticField = computeYHLLZMagneticField(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double starRegionDensity = computeYHLLDensity(topStateVector, bottomStateVector, material1Parameters, material2Parameters);
    double starRegionXVelocity = computeYHLLXMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionYVelocity = computeYHLLYMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionZVelocity = computeYHLLZMomentum(topStateVector, bottomStateVector, material1Parameters, material2Parameters) / starRegionDensity;
    double starRegionMagneticFieldVelocityVectorProduct = (starRegionXMagneticField * starRegionXVelocity) + (starRegionYMagneticField * starRegionYVelocity) +
            (starRegionZMagneticField * starRegionZVelocity);

    double magneticFieldVelocityVectorProduct = (bottomXMagneticField * bottomXVelocity) + (bottomYMagneticField * bottomYVelocity) + (bottomZMagneticField * bottomZVelocity);

    double bottomTotalPressure = computeTotalPressure(bottomStateVector);
    double starRegionTotalPressure = computeBottomStarRegionTotalPressure(topStateVector, bottomStateVector, material1Parameters, material2Parameters);

    double bottomTotalEnergy = bottomStateVector.computeTotalEnergy(material2Parameters);

    double firstTerm = bottomTotalEnergy * ((bottomWaveSpeed - bottomYVelocity) / (bottomWaveSpeed - starRegionWaveSpeed));
    double secondTerm = ((starRegionTotalPressure * starRegionWaveSpeed) - (bottomTotalPressure * bottomYVelocity) -
                         ((starRegionYMagneticField * starRegionMagneticFieldVelocityVectorProduct) - (bottomYMagneticField * magneticFieldVelocityVectorProduct))) /
            (bottomWaveSpeed - starRegionWaveSpeed);

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

double MHDHLLCSolver::computeTopWaveSpeed(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double topYVelocity = topStateVector.getYVelocity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    double topYFastMagnetoAcousticSpeed = topStateVector.computeYFastMagnetoAcousticSpeed(material1Parameters);
    double bottomYFastMagnetoAcousticSpeed = bottomStateVector.computeYFastMagnetoAcousticSpeed(material2Parameters);

    return min(topYVelocity, bottomYVelocity) - max(topYFastMagnetoAcousticSpeed, bottomYFastMagnetoAcousticSpeed);
}

double MHDHLLCSolver::computeBottomWaveSpeed(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double topYVelocity = topStateVector.getYVelocity();
    double bottomYVelocity = bottomStateVector.getYVelocity();

    double topYFastMagnetoAcousticSpeed = topStateVector.computeYFastMagnetoAcousticSpeed(material1Parameters);
    double bottomYFastMagnetoAcousticSpeed = bottomStateVector.computeYFastMagnetoAcousticSpeed(material2Parameters);

    return max(topYVelocity, bottomYVelocity) + max(topYFastMagnetoAcousticSpeed, bottomYFastMagnetoAcousticSpeed);
}
