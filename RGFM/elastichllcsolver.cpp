#include "elastichllcsolver.h"

ElasticHLLCSolver::ElasticHLLCSolver()
{
}

double ElasticHLLCSolver::computeTildeRegionDensity(double density, double waveSpeed, double velocity, double tildeRegionWaveSpeed)
{
    return density * ((waveSpeed - velocity) / (waveSpeed - tildeRegionWaveSpeed));
}

double ElasticHLLCSolver::computeLeftTildeRegionDensity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                        HyperelasticMaterialParameters material2Parameters)
{
    double tildeRegionWaveSpeed = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();

    return computeTildeRegionDensity(leftDensity, leftWaveSpeed, leftXVelocity, tildeRegionWaveSpeed);
}

double ElasticHLLCSolver::computeRightTildeRegionDensity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                         HyperelasticMaterialParameters material2Parameters)
{
    double tildeRegionWaveSpeed = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();

    return computeTildeRegionDensity(rightDensity, rightWaveSpeed, rightXVelocity, tildeRegionWaveSpeed);
}

vector<double> ElasticHLLCSolver::computeLeftTildeRegionConservedVariableVector(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> tildeRegionConservedVariableVector(14);

    double tildeRegionDensity = computeLeftTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > tildeRegionDistortionTensor = computeLeftTildeRegionDistortionTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionTotalEnergy = computeLeftTildeRegionTotalEnergy(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftYVelocity = leftStateVector.getYVelocity();
    double leftZVelocity = leftStateVector.getZVelocity();

    tildeRegionConservedVariableVector[0] = tildeRegionDensity;

    tildeRegionConservedVariableVector[1] = tildeRegionDensity * tildeRegionXVelocity;
    tildeRegionConservedVariableVector[2] = tildeRegionDensity * leftYVelocity;
    tildeRegionConservedVariableVector[3] = tildeRegionDensity * leftZVelocity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            tildeRegionConservedVariableVector[4 + (i * 3) + j] = tildeRegionDistortionTensor[i][j];
        }
    }

    tildeRegionConservedVariableVector[13] = tildeRegionDensity * tildeRegionTotalEnergy;

    return tildeRegionConservedVariableVector;
}

vector<double> ElasticHLLCSolver::computeRightTildeRegionConservedVariableVector(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                 HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> tildeRegionConservedVariableVector(14);

    double tildeRegionDensity = computeRightTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > tildeRegionDistortionTensor = computeRightTildeRegionDistortionTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionTotalEnergy = computeRightTildeRegionTotalEnergy(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightYVelocity = rightStateVector.getYVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    tildeRegionConservedVariableVector[0] = tildeRegionDensity;

    tildeRegionConservedVariableVector[1] = tildeRegionDensity * tildeRegionXVelocity;
    tildeRegionConservedVariableVector[2] = tildeRegionDensity * rightYVelocity;
    tildeRegionConservedVariableVector[3] = tildeRegionDensity * rightZVelocity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            tildeRegionConservedVariableVector[4 + (i * 3) + j] = tildeRegionDistortionTensor[i][j];
        }
    }

    tildeRegionConservedVariableVector[13] = tildeRegionDensity * tildeRegionTotalEnergy;

    return tildeRegionConservedVariableVector;
}

vector<double> ElasticHLLCSolver::computeLeftStarRegionConservedVariableVector(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                               HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> starRegionConservedVariableVector(14);

    double tildeRegionDensity = computeLeftTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > starRegionDistortionTensor = computeLeftStarRegionDistortionTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionTotalEnergy = computeLeftStarRegionTotalEnergy(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionYVelocity = computeXStarRegionYVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZVelocity = computeXStarRegionZVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    starRegionConservedVariableVector[0] = tildeRegionDensity;

    starRegionConservedVariableVector[1] = tildeRegionDensity * tildeRegionXVelocity;
    starRegionConservedVariableVector[2] = tildeRegionDensity * starRegionYVelocity;
    starRegionConservedVariableVector[3] = tildeRegionDensity * starRegionZVelocity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            starRegionConservedVariableVector[4 + (i * 3) + j] = starRegionDistortionTensor[i][j];
        }
    }

    starRegionConservedVariableVector[13] = tildeRegionDensity * starRegionTotalEnergy;

    return starRegionConservedVariableVector;
}

vector<double> ElasticHLLCSolver::computeRightStarRegionConservedVariableVector(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> starRegionConservedVariableVector(14);

    double tildeRegionDensity = computeRightTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > starRegionDistortionTensor = computeRightStarRegionDistortionTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionTotalEnergy = computeRightStarRegionTotalEnergy(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionYVelocity = computeXStarRegionYVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZVelocity = computeXStarRegionZVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    starRegionConservedVariableVector[0] = tildeRegionDensity;

    starRegionConservedVariableVector[1] = tildeRegionDensity * tildeRegionXVelocity;
    starRegionConservedVariableVector[2] = tildeRegionDensity * starRegionYVelocity;
    starRegionConservedVariableVector[3] = tildeRegionDensity * starRegionZVelocity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            starRegionConservedVariableVector[4 + (i * 3) + j] = starRegionDistortionTensor[i][j];
        }
    }

    starRegionConservedVariableVector[13] = tildeRegionDensity * starRegionTotalEnergy;

    return starRegionConservedVariableVector;
}

ElasticStateVector ElasticHLLCSolver::solveTildeX(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                  HyperelasticMaterialParameters material2Parameters)
{
    double tildeRegionWaveSpeed = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
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
    else if (0.0 <= tildeRegionWaveSpeed)
    {
        vector<double> tildeRegionConservedVariableVector = computeLeftTildeRegionConservedVariableVector(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

        ElasticStateVector solution;
        solution.setConservedVariableVector(tildeRegionConservedVariableVector, material1Parameters);

        return solution;
    }
    else
    {
        vector<double> tildeRegionConservedVariableVector = computeRightTildeRegionConservedVariableVector(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

        ElasticStateVector solution;
        solution.setConservedVariableVector(tildeRegionConservedVariableVector, material2Parameters);

        return solution;
    }
}

ElasticStateVector ElasticHLLCSolver::solveStarX(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                 HyperelasticMaterialParameters material2Parameters)
{
    double tildeRegionWaveSpeed = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
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
    else if (0.0 <= tildeRegionWaveSpeed)
    {
        vector<double> starRegionConservedVariableVector = computeLeftStarRegionConservedVariableVector(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

        ElasticStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material1Parameters);

        return solution;
    }
    else
    {
        vector<double> starRegionConservedVariableVector = computeRightStarRegionConservedVariableVector(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

        ElasticStateVector solution;
        solution.setConservedVariableVector(starRegionConservedVariableVector, material2Parameters);

        return solution;
    }
}

double ElasticHLLCSolver::computeTildeRegionXVelocity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                      HyperelasticMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftEntropy = leftStateVector.getEntropy();

    vector<vector<double> > leftDistortionTensor = leftStateVector.getDistortionTensor();
    vector<vector<double> > leftTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(leftDensity, leftDistortionTensor, leftEntropy, material1Parameters);

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightEntropy = rightStateVector.getEntropy();

    vector<vector<double> > rightDistortionTensor = rightStateVector.getDistortionTensor();
    vector<vector<double> > rightTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(rightDensity, rightDistortionTensor, rightEntropy, material2Parameters);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftCoefficient = leftDensity * (leftWaveSpeed - leftXVelocity);
    double rightCoefficient = rightDensity * (rightWaveSpeed - rightXVelocity);

    double numerator = (rightCoefficient * rightXVelocity) - (leftCoefficient * leftXVelocity) + rightTotalStressTensor[0][0] - leftTotalStressTensor[0][0];
    double denominator = rightCoefficient - leftCoefficient;

    return numerator / denominator;
}

vector<vector<double> > ElasticHLLCSolver::computeLeftTildeRegionTotalStressTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                   HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftEntropy = leftStateVector.getEntropy();

    vector<vector<double> > leftDistortionTensor = leftStateVector.getDistortionTensor();
    vector<vector<double> > leftTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(leftDensity, leftDistortionTensor, leftEntropy, material1Parameters);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > tildeRegionTotalStressTensor = leftTotalStressTensor;
    tildeRegionTotalStressTensor[0][0] = leftTotalStressTensor[0][0] + (leftDensity * (leftWaveSpeed - leftXVelocity) * (leftXVelocity - tildeRegionXVelocity));

    return tildeRegionTotalStressTensor;
}

vector<vector<double> > ElasticHLLCSolver::computeRightTildeRegionTotalStressTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                    HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightEntropy = rightStateVector.getEntropy();

    vector<vector<double> > rightDistortionTensor = rightStateVector.getDistortionTensor();
    vector<vector<double> > rightTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(rightDensity, rightDistortionTensor, rightEntropy, material2Parameters);

    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > tildeRegionTotalStressTensor = rightTotalStressTensor;
    tildeRegionTotalStressTensor[0][0] = rightTotalStressTensor[0][0] + (rightDensity * (rightWaveSpeed - rightXVelocity) * (rightXVelocity - tildeRegionXVelocity));

    return tildeRegionTotalStressTensor;
}

vector<vector<double> > ElasticHLLCSolver::computeLeftTildeRegionDistortionTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                  HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    vector<vector<double> > leftDistortionTensor = leftStateVector.getDistortionTensor();

    double leftTildeRegionDensity = computeLeftTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > leftDeformationGradientTensor = MatrixAlgebra::computeInverseMatrix(leftDistortionTensor);
    vector<vector<double> > tildeRegionDeformationGradientTensor = leftDeformationGradientTensor;

    for (int i = 0; i < 3; i++)
    {
        tildeRegionDeformationGradientTensor[0][i] = (leftDensity / leftTildeRegionDensity) * leftDeformationGradientTensor[0][i];
    }

    return MatrixAlgebra::computeInverseMatrix(tildeRegionDeformationGradientTensor);
}

vector<vector<double> > ElasticHLLCSolver::computeRightTildeRegionDistortionTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                   HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double rightDensity = rightStateVector.getDensity();
    vector<vector<double> > rightDistortionTensor = rightStateVector.getDistortionTensor();

    double rightTildeRegionDensity = computeRightTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > rightDeformationGradientTensor = MatrixAlgebra::computeInverseMatrix(rightDistortionTensor);
    vector<vector<double> > tildeRegionDeformationGradientTensor = rightDeformationGradientTensor;

    for (int i = 0; i < 3; i++)
    {
        tildeRegionDeformationGradientTensor[0][i] = (rightDensity / rightTildeRegionDensity) * rightDeformationGradientTensor[0][i];
    }

    return MatrixAlgebra::computeInverseMatrix(tildeRegionDeformationGradientTensor);
}

double ElasticHLLCSolver::computeLeftTildeRegionTotalEnergy(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                            HyperelasticMaterialParameters material2Parameters)
{
    double leftDensity = leftStateVector.getDensity();
    double leftEntropy = leftStateVector.getEntropy();

    double leftXVelocity = leftStateVector.getXVelocity();
    double leftYVelocity = leftStateVector.getYVelocity();
    double leftZVelocity = leftStateVector.getZVelocity();

    vector<vector<double> > leftDistortionTensor = leftStateVector.getDistortionTensor();
    vector<vector<double> > leftTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(leftDensity, leftDistortionTensor, leftEntropy, material1Parameters);

    double leftTotalEnergy = ElasticEquationOfState::computeTotalEnergy(leftDistortionTensor, leftEntropy, leftXVelocity, leftYVelocity, leftZVelocity, material1Parameters);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    vector<vector<double> > leftTildeRegionTotalStressTensor = computeLeftTildeRegionTotalStressTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<double> leftVelocityVector(3);
    leftVelocityVector[0] = leftXVelocity;
    leftVelocityVector[1] = leftYVelocity;
    leftVelocityVector[2] = leftZVelocity;

    vector<double> leftTildeRegionVelocityVector(3);
    leftTildeRegionVelocityVector[0] = tildeRegionXVelocity;
    leftTildeRegionVelocityVector[1] = leftYVelocity;
    leftTildeRegionVelocityVector[2] = leftZVelocity;

    double numerator = (leftDensity * leftTotalEnergy) * (leftWaveSpeed - leftXVelocity) + VectorAlgebra::computeDotProduct(leftTotalStressTensor[0], leftVelocityVector) -
            VectorAlgebra::computeDotProduct(leftTildeRegionTotalStressTensor[0], leftTildeRegionVelocityVector);
    double denominator = leftDensity * (leftWaveSpeed - leftXVelocity);

    return numerator / denominator;
}

double ElasticHLLCSolver::computeRightTildeRegionTotalEnergy(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                             HyperelasticMaterialParameters material2Parameters)
{
    double rightDensity = rightStateVector.getDensity();
    double rightEntropy = rightStateVector.getEntropy();

    double rightXVelocity = rightStateVector.getXVelocity();
    double rightYVelocity = rightStateVector.getYVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    vector<vector<double> > rightDistortionTensor = rightStateVector.getDistortionTensor();
    vector<vector<double> > rightTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(rightDensity, rightDistortionTensor, rightEntropy, material2Parameters);

    double rightTotalEnergy = ElasticEquationOfState::computeTotalEnergy(rightDistortionTensor, rightEntropy, rightXVelocity, rightYVelocity, rightZVelocity, material2Parameters);

    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    vector<vector<double> > rightTildeRegionTotalStressTensor = computeRightTildeRegionTotalStressTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<double> rightVelocityVector(3);
    rightVelocityVector[0] = rightXVelocity;
    rightVelocityVector[1] = rightYVelocity;
    rightVelocityVector[2] = rightZVelocity;

    vector<double> rightTildeRegionVelocityVector(3);
    rightTildeRegionVelocityVector[0] = tildeRegionXVelocity;
    rightTildeRegionVelocityVector[1] = rightYVelocity;
    rightTildeRegionVelocityVector[2] = rightZVelocity;

    double numerator = (rightDensity * rightTotalEnergy) - (rightWaveSpeed - rightXVelocity) + VectorAlgebra::computeDotProduct(rightTotalStressTensor[0], rightVelocityVector) -
            VectorAlgebra::computeDotProduct(rightTildeRegionTotalStressTensor[0], rightTildeRegionVelocityVector);
    double denominator = rightDensity * (rightWaveSpeed - rightXVelocity);

    return numerator / denominator;
}

double ElasticHLLCSolver::computeXStarRegionYVelocity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                      HyperelasticMaterialParameters material2Parameters)
{
    double leftTildeRegionDensity = computeLeftTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightTildeRegionDensity = computeRightTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double leftEntropy = leftStateVector.getEntropy();
    double leftYVelocity = leftStateVector.getYVelocity();

    vector<vector<double> > leftDistortionTensor = leftStateVector.getDistortionTensor();
    vector<vector<double> > leftTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(leftDensity, leftDistortionTensor, leftEntropy, material1Parameters);

    double rightDensity = rightStateVector.getDensity();
    double rightEntropy = rightStateVector.getEntropy();
    double rightYVelocity = rightStateVector.getYVelocity();

    vector<vector<double> > rightDistortionTensor = rightStateVector.getDistortionTensor();
    vector<vector<double> > rightTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(rightDensity, rightDistortionTensor, rightEntropy, material2Parameters);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double numerator = leftTotalStressTensor[0][1] - rightTotalStressTensor[0][1] + (leftTildeRegionDensity * leftYVelocity) * (leftWaveSpeed - tildeRegionXVelocity) -
            (rightTildeRegionDensity * rightYVelocity) * (rightWaveSpeed - tildeRegionXVelocity);
    double denominator = leftTildeRegionDensity * (leftWaveSpeed - tildeRegionXVelocity) - rightTildeRegionDensity * (rightWaveSpeed - tildeRegionXVelocity);

    return numerator / denominator;
}

double ElasticHLLCSolver::computeXStarRegionZVelocity(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                      HyperelasticMaterialParameters material2Parameters)
{
    double leftTildeRegionDensity = computeLeftTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightTildeRegionDensity = computeRightTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftDensity = leftStateVector.getDensity();
    double leftEntropy = leftStateVector.getEntropy();
    double leftZVelocity = leftStateVector.getZVelocity();

    vector<vector<double> > leftDistortionTensor = leftStateVector.getDistortionTensor();
    vector<vector<double> > leftTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(leftDensity, leftDistortionTensor, leftEntropy, material1Parameters);

    double rightDensity = rightStateVector.getDensity();
    double rightEntropy = rightStateVector.getEntropy();
    double rightZVelocity = rightStateVector.getZVelocity();

    vector<vector<double> > rightDistortionTensor = rightStateVector.getDistortionTensor();
    vector<vector<double> > rightTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(rightDensity, rightDistortionTensor, rightEntropy, material2Parameters);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double numerator = leftTotalStressTensor[0][2] - rightTotalStressTensor[0][2] + (leftTildeRegionDensity * leftZVelocity) * (leftWaveSpeed - tildeRegionXVelocity) -
            (rightTildeRegionDensity * rightZVelocity) * (rightWaveSpeed - tildeRegionXVelocity);
    double denominator = leftTildeRegionDensity * (leftWaveSpeed - tildeRegionXVelocity) - rightTildeRegionDensity * (rightWaveSpeed - tildeRegionXVelocity);

    return numerator / denominator;
}

vector<vector<double> > ElasticHLLCSolver::computeLeftStarRegionTotalStressTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                  HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double leftTildeRegionDensity = computeLeftTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftYVelocity = leftStateVector.getYVelocity();
    double leftZVelocity = leftStateVector.getZVelocity();

    vector<vector<double> > leftTildeRegionTotalStressTensor = computeLeftTildeRegionTotalStressTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionYVelocity = computeXStarRegionYVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZVelocity = computeXStarRegionZVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > starRegionTotalStressTensor = leftTildeRegionTotalStressTensor;

    double firstTerm1 = (leftTildeRegionDensity * leftYVelocity) * (leftWaveSpeed - tildeRegionXVelocity) + leftTildeRegionTotalStressTensor[0][1];
    double secondTerm1 = (starRegionYVelocity * leftTildeRegionDensity) * (leftWaveSpeed - tildeRegionXVelocity);
    starRegionTotalStressTensor[0][1] = firstTerm1 - secondTerm1;

    double firstTerm2 = (leftTildeRegionDensity * leftZVelocity) * (leftWaveSpeed - tildeRegionXVelocity) + leftTildeRegionTotalStressTensor[0][2];
    double secondTerm2 = (starRegionZVelocity * leftTildeRegionDensity) * (leftWaveSpeed - tildeRegionXVelocity);
    starRegionTotalStressTensor[0][2] = firstTerm2 - secondTerm2;

    return starRegionTotalStressTensor;
}

vector<vector<double> > ElasticHLLCSolver::computeRightStarRegionTotalStressTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                   HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double rightTildeRegionDensity = computeRightTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightYVelocity = rightStateVector.getYVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    vector<vector<double> > rightTildeRegionTotalStressTensor = computeRightTildeRegionTotalStressTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double starRegionYVelocity = computeXStarRegionYVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZVelocity = computeXStarRegionZVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > starRegionTotalStressTensor = rightTildeRegionTotalStressTensor;

    double firstTerm1 = (rightTildeRegionDensity * rightYVelocity) * (rightWaveSpeed - tildeRegionXVelocity) + rightTildeRegionTotalStressTensor[0][1];
    double secondTerm1 = (starRegionYVelocity * rightTildeRegionDensity) * (rightWaveSpeed - tildeRegionXVelocity);
    starRegionTotalStressTensor[0][1] = firstTerm1 - secondTerm1;

    double firstTerm2 = (rightTildeRegionDensity * rightZVelocity) * (rightWaveSpeed - tildeRegionXVelocity) + rightTildeRegionTotalStressTensor[0][2];
    double secondTerm2 = (starRegionZVelocity * rightTildeRegionDensity) * (rightWaveSpeed - tildeRegionXVelocity);
    starRegionTotalStressTensor[0][2] = firstTerm2 - secondTerm2;

    return starRegionTotalStressTensor;
}

vector<vector<double> > ElasticHLLCSolver::computeLeftStarRegionDistortionTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                 HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<vector<double> > leftTildeRegionDistortionTensor = computeLeftTildeRegionDistortionTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    vector<vector<double> > leftTildeRegionDeformationGradientTensor = MatrixAlgebra::computeInverseMatrix(leftTildeRegionDistortionTensor);

    double leftTildeRegionDensity = computeLeftTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftYVelocity = leftStateVector.getYVelocity();
    double leftZVelocity = leftStateVector.getZVelocity();

    double starRegionYVelocity = computeXStarRegionYVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZVelocity = computeXStarRegionZVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double denominator = leftTildeRegionDensity * (leftWaveSpeed - tildeRegionXVelocity);
    vector<vector<double> > leftStarRegionDeformationGradientTensor = leftTildeRegionDeformationGradientTensor;

    double numerator;
    for (int i = 0; i < 3; i++)
    {
        numerator = (leftTildeRegionDensity * leftTildeRegionDeformationGradientTensor[1][i]) * (leftWaveSpeed - tildeRegionXVelocity) +
                (leftTildeRegionDensity * leftYVelocity * leftTildeRegionDeformationGradientTensor[0][i]) -
                (leftTildeRegionDensity * starRegionYVelocity * leftTildeRegionDeformationGradientTensor[0][i]);

        leftStarRegionDeformationGradientTensor[1][i] = numerator / denominator;
    }

    for (int i = 0; i < 3; i++)
    {
        numerator = (leftTildeRegionDensity * leftTildeRegionDeformationGradientTensor[2][i]) * (leftWaveSpeed - tildeRegionXVelocity) +
                (leftTildeRegionDensity * leftZVelocity * leftTildeRegionDeformationGradientTensor[0][i]) -
                (leftTildeRegionDensity * starRegionZVelocity * leftTildeRegionDeformationGradientTensor[0][i]);

        leftStarRegionDeformationGradientTensor[2][i] = numerator / denominator;
    }

    return MatrixAlgebra::computeInverseMatrix(leftStarRegionDeformationGradientTensor);
}

vector<vector<double> > ElasticHLLCSolver::computeRightStarRegionDistortionTensor(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                                                                  HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<vector<double> > rightTildeRegionDistortionTensor = computeRightTildeRegionDistortionTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    vector<vector<double> > rightTildeRegionDeformationGradientTensor = MatrixAlgebra::computeInverseMatrix(rightTildeRegionDistortionTensor);

    double rightTildeRegionDensity = computeRightTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightYVelocity = rightStateVector.getYVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    double starRegionYVelocity = computeXStarRegionYVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZVelocity = computeXStarRegionZVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double denominator = rightTildeRegionDensity * (rightWaveSpeed - tildeRegionXVelocity);
    vector<vector<double> > rightStarRegionDeformationGradientTensor = rightTildeRegionDeformationGradientTensor;

    double numerator;
    for (int i = 0; i < 3; i++)
    {
        numerator = (rightTildeRegionDensity * rightTildeRegionDeformationGradientTensor[1][i]) * (rightWaveSpeed - tildeRegionXVelocity) +
                (rightTildeRegionDensity * rightYVelocity * rightTildeRegionDeformationGradientTensor[0][i]) -
                (rightTildeRegionDensity * starRegionYVelocity * rightTildeRegionDeformationGradientTensor[0][i]);

        rightStarRegionDeformationGradientTensor[1][i] = numerator / denominator;
    }

    for (int i = 0; i < 3; i++)
    {
        numerator = (rightTildeRegionDensity * rightTildeRegionDeformationGradientTensor[2][i]) * (rightWaveSpeed - tildeRegionXVelocity) +
                (rightTildeRegionDensity * rightZVelocity * rightTildeRegionDeformationGradientTensor[0][i]) -
                (rightTildeRegionDensity * starRegionZVelocity * rightTildeRegionDeformationGradientTensor[0][i]);

        rightStarRegionDeformationGradientTensor[2][i] = numerator / denominator;
    }

    return MatrixAlgebra::computeInverseMatrix(rightStarRegionDeformationGradientTensor);
}

double ElasticHLLCSolver::computeLeftStarRegionTotalEnergy(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                           HyperelasticMaterialParameters material2Parameters)
{
    double leftTildeRegionDensity = computeLeftTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double leftTildeRegionTotalEnergy = computeLeftTildeRegionTotalEnergy(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftYVelocity = leftStateVector.getYVelocity();
    double leftZVelocity = leftStateVector.getZVelocity();

    double starRegionYVelocity = computeXStarRegionYVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZVelocity = computeXStarRegionZVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > leftTildeRegionTotalStressTensor = computeLeftTildeRegionTotalStressTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    vector<vector<double> > leftStarRegionTotalStressTensor = computeLeftStarRegionTotalStressTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double leftWaveSpeed = computeLeftWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<double> leftTildeRegionVelocityVector(3);
    leftTildeRegionVelocityVector[0] = tildeRegionXVelocity;
    leftTildeRegionVelocityVector[1] = leftYVelocity;
    leftTildeRegionVelocityVector[2] = leftZVelocity;

    vector<double> leftStarRegionVelocityVector(3);
    leftStarRegionVelocityVector[0] = tildeRegionXVelocity;
    leftStarRegionVelocityVector[1] = starRegionYVelocity;
    leftStarRegionVelocityVector[2] = starRegionZVelocity;

    double numerator = (leftTildeRegionDensity * leftTildeRegionTotalEnergy) * (leftWaveSpeed - tildeRegionXVelocity) +
            VectorAlgebra::computeDotProduct(leftTildeRegionTotalStressTensor[0], leftTildeRegionVelocityVector) -
            VectorAlgebra::computeDotProduct(leftStarRegionTotalStressTensor[0], leftStarRegionVelocityVector);
    double denominator = leftTildeRegionDensity * (leftWaveSpeed - tildeRegionXVelocity);

    return numerator / denominator;
}

double ElasticHLLCSolver::computeRightStarRegionTotalEnergy(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                            HyperelasticMaterialParameters material2Parameters)
{
    double rightTildeRegionDensity = computeRightTildeRegionDensity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double rightTildeRegionTotalEnergy = computeRightTildeRegionTotalEnergy(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightYVelocity = rightStateVector.getYVelocity();
    double rightZVelocity = rightStateVector.getZVelocity();

    double starRegionYVelocity = computeXStarRegionYVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double starRegionZVelocity = computeXStarRegionZVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<vector<double> > rightTildeRegionTotalStressTensor = computeRightTildeRegionTotalStressTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    vector<vector<double> > rightStarRegionTotalStressTensor = computeRightStarRegionTotalStressTensor(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    double rightWaveSpeed = computeRightWaveSpeed(leftStateVector, rightStateVector, material1Parameters, material2Parameters);
    double tildeRegionXVelocity = computeTildeRegionXVelocity(leftStateVector, rightStateVector, material1Parameters, material2Parameters);

    vector<double> rightTildeRegionVelocityVector(3);
    rightTildeRegionVelocityVector[0] = tildeRegionXVelocity;
    rightTildeRegionVelocityVector[1] = rightYVelocity;
    rightTildeRegionVelocityVector[2] = rightZVelocity;

    vector<double> rightStarRegionVelocityVector(3);
    rightStarRegionVelocityVector[0] = tildeRegionXVelocity;
    rightStarRegionVelocityVector[1] = starRegionYVelocity;
    rightStarRegionVelocityVector[2] = starRegionZVelocity;

    double numerator = (rightTildeRegionDensity * rightTildeRegionTotalEnergy) * (rightWaveSpeed - tildeRegionXVelocity) +
            VectorAlgebra::computeDotProduct(rightTildeRegionTotalStressTensor[0], rightTildeRegionVelocityVector) -
            VectorAlgebra::computeDotProduct(rightStarRegionTotalStressTensor[0], rightStarRegionVelocityVector);
    double denominator = rightTildeRegionDensity * (rightWaveSpeed - tildeRegionXVelocity);

    return numerator / denominator;
}

double ElasticHLLCSolver::computeLeftWaveSpeed(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                               HyperelasticMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double rightXVelocity = rightStateVector.getXVelocity();

    double leftXSoundSpeed = leftStateVector.computeSoundSpeed(material1Parameters, 0);
    double rightXSoundSpeed = rightStateVector.computeSoundSpeed(material2Parameters, 0);

    return min(leftXVelocity, rightXVelocity) - max(leftXSoundSpeed, rightXSoundSpeed);
}

double ElasticHLLCSolver::computeRightWaveSpeed(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, HyperelasticMaterialParameters material1Parameters,
                                                HyperelasticMaterialParameters material2Parameters)
{
    double leftXVelocity = leftStateVector.getXVelocity();
    double rightXVelocity = rightStateVector.getXVelocity();

    double leftXSoundSpeed = leftStateVector.computeSoundSpeed(material1Parameters, 0);
    double rightXSoundSpeed = rightStateVector.computeSoundSpeed(material2Parameters, 0);

    return max(leftXVelocity, rightXVelocity) + max(leftXSoundSpeed, rightXSoundSpeed);
}
