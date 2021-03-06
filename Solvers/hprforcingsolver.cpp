#include "hprforcingsolver.h"

HPRForcingSolver::HPRForcingSolver()
{
}

vector<double> HPRForcingSolver::evolveConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                               double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    HPRStateVector leftStateVector;
    HPRStateVector middleStateVector;
    HPRStateVector rightStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, materialParameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, materialParameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, HPRStateVector::computeSourceTermVector(middleConservedVariableVector, materialParameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, materialParameters);
    HPRStateVector middleStateVectorFirstStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                    slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, HPRStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved, materialParameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, materialParameters);
    HPRStateVector middleStateVectorSecondStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                     slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, HPRStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved, materialParameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, materialParameters);
    HPRStateVector middleStateVectorThirdStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(1.0, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                    slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, HPRStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, materialParameters));

    return computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

vector<double> HPRForcingSolver::evolveIntermediateConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                           vector<double> rightConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                           HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRIntermediateStateVector leftStateVector;
    HPRIntermediateStateVector middleStateVector;
    HPRIntermediateStateVector rightStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, material1Parameters, material2Parameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, HPRIntermediateStateVector::computeSourceTermVector(middleConservedVariableVector, material1Parameters,
                                                                                                                           material2Parameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, material1Parameters, material2Parameters);
    HPRIntermediateStateVector middleStateVectorFirstStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                                timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, HPRIntermediateStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved, material1Parameters,
                                                                                                                            material2Parameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, material1Parameters, material2Parameters);
    HPRIntermediateStateVector middleStateVectorSecondStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                                 timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, HPRIntermediateStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved, material1Parameters,
                                                                                                                           material2Parameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, material1Parameters, material2Parameters);
    HPRIntermediateStateVector middleStateVectorThirdStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(1.0, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                                timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, HPRIntermediateStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, material1Parameters,
                                                                                                                            material2Parameters));

    return computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

vector<double> HPRForcingSolver::evolveReducedConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                      vector<double> rightConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                      HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRReducedStateVector leftStateVector;
    HPRReducedStateVector middleStateVector;
    HPRReducedStateVector rightStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, material1Parameters, material2Parameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, HPRReducedStateVector::computeSourceTermVector(middleConservedVariableVector, material1Parameters, material2Parameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorFirstStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                           bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, HPRReducedStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved, material1Parameters,
                                                                                                                       material2Parameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorSecondStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                            bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, HPRReducedStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved, material1Parameters,
                                                                                                                      material2Parameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorThirdStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(1.0, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                           bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, HPRReducedStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, material1Parameters,
                                                                                                                       material2Parameters));

    return computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);

}

vector<double> HPRForcingSolver::evolveConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                                 vector<double> topConservedVariableVector, vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep,
                                                                 double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    HPRStateVector leftStateVector;
    HPRStateVector middleStateVector;
    HPRStateVector rightStateVector;

    HPRStateVector topStateVector;
    HPRStateVector bottomStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, materialParameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, materialParameters);

    topStateVector.setConservedVariableVector(topConservedVariableVector, materialParameters);
    bottomStateVector.setConservedVariableVector(bottomConservedVariableVector, materialParameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, HPRStateVector::computeSourceTermVector(middleConservedVariableVector, materialParameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, materialParameters);

    HPRStateVector middleStateVectorFirstStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                     slopeLimiter, materialParameters);
    HPRStateVector middleStateVectorFirstStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorFirstStepEvolved1, bottomStateVector, cellSpacing,
                                                                                                    timeStep, bias, slopeLimiter, materialParameters);
    HPRStateVector middleStateVectorFirstStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorFirstStepEvolved2, rightStateVector, cellSpacing,
                                                                                                    timeStep, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, HPRStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved, materialParameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, materialParameters);

    HPRStateVector middleStateVectorSecondStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                      slopeLimiter, materialParameters);
    HPRStateVector middleStateVectorSecondStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorSecondStepEvolved1, bottomStateVector, cellSpacing,
                                                                                                     timeStep, bias, slopeLimiter, materialParameters);
    HPRStateVector middleStateVectorSecondStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorSecondStepEvolved2, rightStateVector, cellSpacing,
                                                                                                     timeStep, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, HPRStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved, materialParameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, materialParameters);

    HPRStateVector middleStateVectorThirdStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                     slopeLimiter, materialParameters);
    HPRStateVector middleStateVectorThirdStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(1.0, topStateVector, middleStateVectorThirdStepEvolved1, bottomStateVector, cellSpacing,
                                                                                                    timeStep, bias, slopeLimiter, materialParameters);
    HPRStateVector middleStateVectorThirdStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVectorThirdStepEvolved2, rightStateVector, cellSpacing,
                                                                                                    timeStep, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, HPRStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, materialParameters));

    return computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

vector<double> HPRForcingSolver::evolveIntermediateConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                             vector<double> rightConservedVariableVector, vector<double> topConservedVariableVector,
                                                                             vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                             HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRIntermediateStateVector leftStateVector;
    HPRIntermediateStateVector middleStateVector;
    HPRIntermediateStateVector rightStateVector;

    HPRIntermediateStateVector topStateVector;
    HPRIntermediateStateVector bottomStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, material1Parameters, material2Parameters);

    topStateVector.setConservedVariableVector(topConservedVariableVector, material1Parameters, material2Parameters);
    bottomStateVector.setConservedVariableVector(bottomConservedVariableVector, material1Parameters, material2Parameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, HPRIntermediateStateVector::computeSourceTermVector(middleConservedVariableVector, material1Parameters,
                                                                                                                           material2Parameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, material1Parameters, material2Parameters);

    HPRIntermediateStateVector middleStateVectorFirstStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                                 timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRIntermediateStateVector middleStateVectorFirstStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorFirstStepEvolved1, bottomStateVector,
                                                                                                                 cellSpacing, timeStep, bias, slopeLimiter, material1Parameters,
                                                                                                                 material2Parameters);
    HPRIntermediateStateVector middleStateVectorFirstStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorFirstStepEvolved2, rightStateVector,
                                                                                                                cellSpacing, timeStep, bias, slopeLimiter, material1Parameters,
                                                                                                                material2Parameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, HPRIntermediateStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved,
                                                                                                                            material1Parameters, material2Parameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, material1Parameters, material2Parameters);

    HPRIntermediateStateVector middleStateVectorSecondStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                                  timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRIntermediateStateVector middleStateVectorSecondStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorSecondStepEvolved1,
                                                                                                                  bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                                  material1Parameters, material2Parameters);
    HPRIntermediateStateVector middleStateVectorSecondStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorSecondStepEvolved2,
                                                                                                                 rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                                 material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, HPRIntermediateStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved,
                                                                                                                           material1Parameters, material2Parameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, material1Parameters, material2Parameters);

    HPRIntermediateStateVector middleStateVectorThirdStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                                 timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRIntermediateStateVector middleStateVectorThirdStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(1.0, topStateVector, middleStateVectorThirdStepEvolved1,
                                                                                                                 bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                                 material1Parameters, material2Parameters);
    HPRIntermediateStateVector middleStateVectorThirdStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVectorThirdStepEvolved2, rightStateVector,
                                                                                                                cellSpacing, timeStep, bias, slopeLimiter, material1Parameters,
                                                                                                                material2Parameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, HPRIntermediateStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, material1Parameters,
                                                                                                                            material2Parameters));

    return computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

vector<double> HPRForcingSolver::evolveReducedConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                        vector<double> rightConservedVariableVector, vector<double> topConservedVariableVector,
                                                                        vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                        HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRReducedStateVector leftStateVector;
    HPRReducedStateVector middleStateVector;
    HPRReducedStateVector rightStateVector;

    HPRReducedStateVector topStateVector;
    HPRReducedStateVector bottomStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, material1Parameters, material2Parameters);

    topStateVector.setConservedVariableVector(topConservedVariableVector, material1Parameters, material2Parameters);
    bottomStateVector.setConservedVariableVector(bottomConservedVariableVector, material1Parameters, material2Parameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, HPRReducedStateVector::computeSourceTermVector(middleConservedVariableVector, material1Parameters, material2Parameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, material1Parameters, material2Parameters);

    HPRReducedStateVector middleStateVectorFirstStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                            bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorFirstStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorFirstStepEvolved1, bottomStateVector,
                                                                                                            cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorFirstStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorFirstStepEvolved2, rightStateVector,
                                                                                                           cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, HPRReducedStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved, material1Parameters,
                                                                                                                       material2Parameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, material1Parameters, material2Parameters);

    HPRReducedStateVector middleStateVectorSecondStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                             bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorSecondStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorSecondStepEvolved1, bottomStateVector,
                                                                                                             cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorSecondStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorSecondStepEvolved2, rightStateVector,
                                                                                                            cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, HPRReducedStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved, material1Parameters,
                                                                                                                      material2Parameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, material1Parameters, material2Parameters);

    HPRReducedStateVector middleStateVectorThirdStepEvolved1 = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                            timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorThirdStepEvolved2 = HPRSolvers::evolveStateByFractionalYTimeStep(1.0, topStateVector, middleStateVectorThirdStepEvolved1, bottomStateVector,
                                                                                                            cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    HPRReducedStateVector middleStateVectorThirdStepEvolved = HPRSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVectorThirdStepEvolved2, rightStateVector,
                                                                                                           cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, HPRReducedStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, material1Parameters,
                                                                                                                       material2Parameters));

    return computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

vector<double> HPRForcingSolver::computeEvolvedConservedVariableVector(vector<double> middleConservedVariableVector, vector<double> firstStep, vector<double> secondStep, vector<double> thirdStep,
                                                                       vector<double> fourthStep)
{
    return VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(
                                         (1.0 / 6.0), VectorAlgebra::addVectors(firstStep, VectorAlgebra::addVectors(
                                                                                    VectorAlgebra::multiplyVector(2.0, secondStep), VectorAlgebra::addVectors(
                                                                                        VectorAlgebra::multiplyVector(2.0, thirdStep), fourthStep)))));
}

void HPRForcingSolver::computeRungeKuttaTimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                                 int slopeLimiter, HPRMaterialParameters materialParameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> leftConservedVariableVector = currentCellsWithBoundary[i].computeConservedVariableVector(materialParameters);
        vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1].computeConservedVariableVector(materialParameters);
        vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 2].computeConservedVariableVector(materialParameters);

        currentCells[i].setConservedVariableVector(evolveConservedVariableVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, cellSpacing,
                                                                                 timeStep, bias, slopeLimiter, materialParameters), materialParameters);
    }
}

void HPRForcingSolver::computeRungeKuttaTimeStep(vector<HPRIntermediateStateVector> & currentCells, vector<HPRIntermediateStateVector> & currentCellsWithBoundary, double cellSpacing,
                                                 double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> leftConservedVariableVector = currentCellsWithBoundary[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 2].computeConservedVariableVector(material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(evolveIntermediateConservedVariableVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, cellSpacing,
                                                                                             timeStep, bias, slopeLimiter, material1Parameters, material2Parameters), material1Parameters,
                                                   material2Parameters);
    }
}

void HPRForcingSolver::computeRungeKuttaTimeStep(vector<HPRReducedStateVector> & currentCells, vector<HPRReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                 double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> leftConservedVariableVector = currentCellsWithBoundary[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 2].computeConservedVariableVector(material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(evolveReducedConservedVariableVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, cellSpacing,
                                                                                        timeStep, bias, slopeLimiter, material1Parameters, material2Parameters), material1Parameters,
                                                   material2Parameters);
    }
}

void HPRForcingSolver::computeRungeKuttaTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                                   double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> leftConservedVariableVector = currentCellsWithBoundary[i + 1][j].computeConservedVariableVector(materialParameters);
            vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(materialParameters);
            vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 1][j + 2].computeConservedVariableVector(materialParameters);

            vector<double> topConservedVariableVector = currentCellsWithBoundary[i][j + 1].computeConservedVariableVector(materialParameters);
            vector<double> bottomConservedVariableVector = currentCellsWithBoundary[i + 2][j + 1].computeConservedVariableVector(materialParameters);

            currentCells[i][j].setConservedVariableVector(evolveConservedVariableVector2D(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector,
                                                                                          topConservedVariableVector, bottomConservedVariableVector, cellSpacing, timeStep, bias,
                                                                                          slopeLimiter, materialParameters), materialParameters);
        }
    }
}

void HPRForcingSolver::computeRungeKuttaTimeStep2D(vector<vector<HPRIntermediateStateVector> > & currentCells, vector<vector<HPRIntermediateStateVector> > & currentCellsWithBoundary,
                                                   double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters,
                                                   HPRMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> leftConservedVariableVector = currentCellsWithBoundary[i + 1][j].computeConservedVariableVector(material1Parameters, material2Parameters);
            vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(material1Parameters, material2Parameters);
            vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 1][j + 2].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> topConservedVariableVector = currentCellsWithBoundary[i][j + 1].computeConservedVariableVector(material1Parameters, material2Parameters);
            vector<double> bottomConservedVariableVector = currentCellsWithBoundary[i + 2][j + 1].computeConservedVariableVector(material1Parameters, material2Parameters);

            currentCells[i][j].setConservedVariableVector(evolveIntermediateConservedVariableVector2D(leftConservedVariableVector, middleConservedVariableVector,
                                                                                                      rightConservedVariableVector, topConservedVariableVector,
                                                                                                      bottomConservedVariableVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                      material1Parameters, material2Parameters), material1Parameters,
                                                          material2Parameters);
        }
    }
}

void HPRForcingSolver::computeRungeKuttaTimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, vector<vector<HPRReducedStateVector> > & currentCellsWithBoundary,
                                                   double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters,
                                                   HPRMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            vector<double> leftConservedVariableVector = currentCellsWithBoundary[i + 1][j].computeConservedVariableVector(material1Parameters, material2Parameters);
            vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1][j + 1].computeConservedVariableVector(material1Parameters, material2Parameters);
            vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 1][j + 2].computeConservedVariableVector(material1Parameters, material2Parameters);

            vector<double> topConservedVariableVector = currentCellsWithBoundary[i][j + 1].computeConservedVariableVector(material1Parameters, material2Parameters);
            vector<double> bottomConservedVariableVector = currentCellsWithBoundary[i + 2][j + 1].computeConservedVariableVector(material1Parameters, material2Parameters);

            currentCells[i][j].setConservedVariableVector(evolveReducedConservedVariableVector2D(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector,
                                                                                                 topConservedVariableVector, bottomConservedVariableVector, cellSpacing, timeStep, bias,
                                                                                                 slopeLimiter, material1Parameters, material2Parameters), material1Parameters,
                                                          material2Parameters);
        }
    }
}
