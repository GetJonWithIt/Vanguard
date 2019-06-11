#include "mhdforcingsolver.h"

MHDForcingSolver::MHDForcingSolver()
{
}

vector<double> MHDForcingSolver::evolveConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                               double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters)
{
    MHDStateVector leftStateVector;
    MHDStateVector middleStateVector;
    MHDStateVector rightStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, materialParameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, materialParameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, MHDStateVector::computeSourceTermVector(middleConservedVariableVector, materialParameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, materialParameters);
    MHDStateVector middleStateVectorFirstStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                    slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, MHDStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved, materialParameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, materialParameters);
    MHDStateVector middleStateVectorSecondStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                     slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, MHDStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved, materialParameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, materialParameters);
    MHDStateVector middleStateVectorThirdStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(1.0, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                    slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, MHDStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, materialParameters));

    return HPRForcingSolver::computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

vector<double> MHDForcingSolver::evolveReducedConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                      vector<double> rightConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                      MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    MHDReducedStateVector leftStateVector;
    MHDReducedStateVector middleStateVector;
    MHDReducedStateVector rightStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, material1Parameters, material2Parameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, MHDReducedStateVector::computeSourceTermVector(middleConservedVariableVector, material1Parameters, material2Parameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorFirstStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                           bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, MHDReducedStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved, material1Parameters,
                                                                                                                       material2Parameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorSecondStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                            bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, MHDReducedStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved, material1Parameters,
                                                                                                                      material2Parameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorThirdStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(1.0, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep,
                                                                                                           bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, MHDReducedStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, material1Parameters,
                                                                                                                       material2Parameters));

    return HPRForcingSolver::computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

vector<double> MHDForcingSolver::evolveConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                                 vector<double> topConservedVariableVector, vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep,
                                                                 double bias, int slopeLimiter, MHDMaterialParameters materialParameters)
{
    MHDStateVector leftStateVector;
    MHDStateVector middleStateVector;
    MHDStateVector rightStateVector;

    MHDStateVector topStateVector;
    MHDStateVector bottomStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, materialParameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, materialParameters);

    topStateVector.setConservedVariableVector(topConservedVariableVector, materialParameters);
    bottomStateVector.setConservedVariableVector(bottomConservedVariableVector, materialParameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, MHDStateVector::computeSourceTermVector(middleConservedVariableVector, materialParameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, materialParameters);

    MHDStateVector middleStateVectorFirstStepEvolved1 = MHDSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                     slopeLimiter, materialParameters);
    MHDStateVector middleStateVectorFirstStepEvolved2 = MHDSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorFirstStepEvolved1, bottomStateVector, cellSpacing,
                                                                                                     timeStep, bias, slopeLimiter, materialParameters);
    MHDStateVector middleStateVectorFirstStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorFirstStepEvolved2, rightStateVector, cellSpacing,
                                                                                                    timeStep, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, MHDStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved, materialParameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, materialParameters);

    MHDStateVector middleStateVectorSecondStepEvolved1 = MHDSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                      slopeLimiter, materialParameters);
    MHDStateVector middleStateVectorSecondStepEvolved2 = MHDSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorSecondStepEvolved1, bottomStateVector, cellSpacing,
                                                                                                      timeStep, bias, slopeLimiter, materialParameters);
    MHDStateVector middleStateVectorSecondStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorSecondStepEvolved2, rightStateVector, cellSpacing,
                                                                                                     timeStep, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, MHDStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved, materialParameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, materialParameters);

    MHDStateVector middleStateVectorThirdStepEvolved1 = MHDSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing, timeStep, bias,
                                                                                                     slopeLimiter, materialParameters);
    MHDStateVector middleStateVectorThirdStepEvolved2 = MHDSolvers::evolveStateByFractionalYTimeStep(1.0, topStateVector, middleStateVectorThirdStepEvolved1, bottomStateVector, cellSpacing,
                                                                                                     timeStep, bias, slopeLimiter, materialParameters);
    MHDStateVector middleStateVectorThirdStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVectorThirdStepEvolved2, rightStateVector, cellSpacing,
                                                                                                    timeStep, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(materialParameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, MHDStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved, materialParameters));

    return HPRForcingSolver::computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

vector<double> MHDForcingSolver::evolveReducedConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                        vector<double> rightConservedVariableVector, vector<double> topConservedVariableVector,
                                                                        vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep, double bias,
                                                                        int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    MHDReducedStateVector leftStateVector;
    MHDReducedStateVector middleStateVector;
    MHDReducedStateVector rightStateVector;

    MHDReducedStateVector topStateVector;
    MHDReducedStateVector bottomStateVector;

    leftStateVector.setConservedVariableVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    rightStateVector.setConservedVariableVector(rightConservedVariableVector, material1Parameters, material2Parameters);

    topStateVector.setConservedVariableVector(topConservedVariableVector, material1Parameters, material2Parameters);
    bottomStateVector.setConservedVariableVector(bottomConservedVariableVector, material1Parameters, material2Parameters);

    vector<double> firstStep = VectorAlgebra::multiplyVector(timeStep, MHDReducedStateVector::computeSourceTermVector(middleConservedVariableVector, material1Parameters,
                                                                                                                      material2Parameters));

    vector<double> middleConservedVariableVectorFirstStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, firstStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorFirstStep, material1Parameters, material2Parameters);

    MHDReducedStateVector middleStateVectorFirstStepEvolved1 = MHDSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                            timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorFirstStepEvolved2 = MHDSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorFirstStepEvolved1,
                                                                                                            bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                            material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorFirstStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorFirstStepEvolved2,
                                                                                                           rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                           material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorFirstStepEvolved = middleStateVectorFirstStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> secondStep = VectorAlgebra::multiplyVector(timeStep, MHDReducedStateVector::computeSourceTermVector(middleConservedVariableVectorFirstStepEvolved,
                                                                                                                       material1Parameters, material2Parameters));

    vector<double> middleConservedVariableVectorSecondStep = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, secondStep));
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorSecondStep, material1Parameters, material2Parameters);

    MHDReducedStateVector middleStateVectorSecondStepEvolved1 = MHDSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                             timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorSecondStepEvolved2 = MHDSolvers::evolveStateByFractionalYTimeStep(0.5, topStateVector, middleStateVectorSecondStepEvolved1,
                                                                                                             bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                             material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorSecondStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.25, leftStateVector, middleStateVectorSecondStepEvolved2,
                                                                                                            rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                            material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorSecondStepEvolved = middleStateVectorSecondStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> thirdStep = VectorAlgebra::multiplyVector(timeStep, MHDReducedStateVector::computeSourceTermVector(middleConservedVariableVectorSecondStepEvolved,
                                                                                                                      material1Parameters, material2Parameters));

    vector<double> middleConservedVariableVectorThirdStep = VectorAlgebra::addVectors(middleConservedVariableVector, thirdStep);
    middleStateVector.setConservedVariableVector(middleConservedVariableVectorThirdStep, material1Parameters, material2Parameters);

    MHDReducedStateVector middleStateVectorThirdStepEvolved1 = MHDSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVector, rightStateVector, cellSpacing,
                                                                                                            timeStep, bias, slopeLimiter, material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorThirdStepEvolved2 = MHDSolvers::evolveStateByFractionalYTimeStep(1.0, topStateVector, middleStateVectorThirdStepEvolved1,
                                                                                                            bottomStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                            material1Parameters, material2Parameters);
    MHDReducedStateVector middleStateVectorThirdStepEvolved = MHDSolvers::evolveStateByFractionalXTimeStep(0.5, leftStateVector, middleStateVectorThirdStepEvolved2,
                                                                                                           rightStateVector, cellSpacing, timeStep, bias, slopeLimiter,
                                                                                                           material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVectorThirdStepEvolved = middleStateVectorThirdStepEvolved.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> fourthStep = VectorAlgebra::multiplyVector(timeStep, MHDReducedStateVector::computeSourceTermVector(middleConservedVariableVectorThirdStepEvolved,
                                                                                                                       material1Parameters, material2Parameters));

    return HPRForcingSolver::computeEvolvedConservedVariableVector(middleConservedVariableVector, firstStep, secondStep, thirdStep, fourthStep);
}

void MHDForcingSolver::computeRungeKuttaTimeStep(vector<MHDStateVector> & currentCells, vector<MHDStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                                 int slopeLimiter, MHDMaterialParameters materialParameters)
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

void MHDForcingSolver::computeRungeKuttaTimeStep(vector<MHDReducedStateVector> & currentCells, vector<MHDReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                                 double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        vector<double> leftConservedVariableVector = currentCellsWithBoundary[i].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> middleConservedVariableVector = currentCellsWithBoundary[i + 1].computeConservedVariableVector(material1Parameters, material2Parameters);
        vector<double> rightConservedVariableVector = currentCellsWithBoundary[i + 2].computeConservedVariableVector(material1Parameters, material2Parameters);

        currentCells[i].setConservedVariableVector(evolveReducedConservedVariableVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector,
                                                                                        cellSpacing, timeStep, bias, slopeLimiter, material1Parameters, material2Parameters),
                                                   material1Parameters, material2Parameters);
    }
}

void MHDForcingSolver::computeRungeKuttaTimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                                   double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters)
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

void MHDForcingSolver::computeRungeKuttaTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary,
                                                   double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters)
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
                                                                                                 slopeLimiter, material1Parameters, material2Parameters),
                                                          material1Parameters, material2Parameters);
        }
    }
}
