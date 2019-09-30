#include "elasticrgfmtests.h"

ElasticRGFMTests::ElasticRGFMTests()
{
}

void ElasticRGFMTests::solveBartonTest1Tilde(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    vector<ElasticStateVector> material1Cells(cellCount);
    vector<ElasticStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    HyperelasticMaterialParameters materialParameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationGradientTensor(3, vector<double>(3));
    leftDeformationGradientTensor[0][0] = 0.98;
    leftDeformationGradientTensor[0][1] = 0.0;
    leftDeformationGradientTensor[0][2] = 0.0;
    leftDeformationGradientTensor[1][0] = 0.02;
    leftDeformationGradientTensor[1][1] = 1.0;
    leftDeformationGradientTensor[1][2] = 0.1;
    leftDeformationGradientTensor[2][0] = 0.0;
    leftDeformationGradientTensor[2][1] = 0.0;
    leftDeformationGradientTensor[2][2] = 1.0;

    vector<vector<double> > rightDeformationGradientTensor(3, vector<double>(3));
    rightDeformationGradientTensor[0][0] = 1.0;
    rightDeformationGradientTensor[0][1] = 0.0;
    rightDeformationGradientTensor[0][2] = 0.0;
    rightDeformationGradientTensor[1][0] = 0.0;
    rightDeformationGradientTensor[1][1] = 1.0;
    rightDeformationGradientTensor[1][2] = 0.1;
    rightDeformationGradientTensor[2][0] = 0.0;
    rightDeformationGradientTensor[2][1] = 0.0;
    rightDeformationGradientTensor[2][2] = 1.0;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeInverseMatrix(leftDeformationGradientTensor);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeInverseMatrix(rightDeformationGradientTensor);

    for (int i = 0; i < cellCount; i++)
    {
        material1Cells[i] = ElasticStateVector(0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), materialParameters);
        material2Cells[i] = ElasticStateVector(0.0, 0.0, 0.0, rightDistortionTensor, 0.0, materialParameters);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(ElasticRGFMSolver::solveTilde(ElasticMultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.06, 0.0, 0, 0, 0, materialParameters,
                                                 materialParameters));
}

void ElasticRGFMTests::solveBartonTest1Star(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    vector<ElasticStateVector> material1Cells(cellCount);
    vector<ElasticStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    HyperelasticMaterialParameters materialParameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationGradientTensor(3, vector<double>(3));
    leftDeformationGradientTensor[0][0] = 0.98;
    leftDeformationGradientTensor[0][1] = 0.0;
    leftDeformationGradientTensor[0][2] = 0.0;
    leftDeformationGradientTensor[1][0] = 0.02;
    leftDeformationGradientTensor[1][1] = 1.0;
    leftDeformationGradientTensor[1][2] = 0.1;
    leftDeformationGradientTensor[2][0] = 0.0;
    leftDeformationGradientTensor[2][1] = 0.0;
    leftDeformationGradientTensor[2][2] = 1.0;

    vector<vector<double> > rightDeformationGradientTensor(3, vector<double>(3));
    rightDeformationGradientTensor[0][0] = 1.0;
    rightDeformationGradientTensor[0][1] = 0.0;
    rightDeformationGradientTensor[0][2] = 0.0;
    rightDeformationGradientTensor[1][0] = 0.0;
    rightDeformationGradientTensor[1][1] = 1.0;
    rightDeformationGradientTensor[1][2] = 0.1;
    rightDeformationGradientTensor[2][0] = 0.0;
    rightDeformationGradientTensor[2][1] = 0.0;
    rightDeformationGradientTensor[2][2] = 1.0;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeInverseMatrix(leftDeformationGradientTensor);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeInverseMatrix(rightDeformationGradientTensor);

    for (int i = 0; i < cellCount; i++)
    {
        material1Cells[i] = ElasticStateVector(0.0, 0.5, 1.0, leftDistortionTensor, pow(10.0, -3.0), materialParameters);
        material2Cells[i] = ElasticStateVector(0.0, 0.0, 0.0, rightDistortionTensor, 0.0, materialParameters);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(ElasticRGFMSolver::solveStar(ElasticMultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.06, 0.0, 0, 0, 0, materialParameters,
                                                materialParameters));
}

void ElasticRGFMTests::solveZhangTestTilde(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    vector<ElasticStateVector> material1Cells(cellCount);
    vector<ElasticStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    HyperelasticMaterialParameters material1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HyperelasticMaterialParameters material2Parameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationGradientTensor(3, vector<double>(3));
    leftDeformationGradientTensor[0][0] = 1.0;
    leftDeformationGradientTensor[0][1] = 0.0;
    leftDeformationGradientTensor[0][2] = 0.0;
    leftDeformationGradientTensor[1][0] = -0.01;
    leftDeformationGradientTensor[1][1] = 0.95;
    leftDeformationGradientTensor[1][2] = 0.02;
    leftDeformationGradientTensor[2][0] = -0.015;
    leftDeformationGradientTensor[2][1] = 0.0;
    leftDeformationGradientTensor[2][2] = 0.9;

    vector<vector<double> > rightDeformationGradientTensor(3, vector<double>(3));
    rightDeformationGradientTensor[0][0] = 1.0;
    rightDeformationGradientTensor[0][1] = 0.0;
    rightDeformationGradientTensor[0][2] = 0.0;
    rightDeformationGradientTensor[1][0] = 0.015;
    rightDeformationGradientTensor[1][1] = 0.95;
    rightDeformationGradientTensor[1][2] = 0.0;
    rightDeformationGradientTensor[2][0] = -0.01;
    rightDeformationGradientTensor[2][1] = 0.0;
    rightDeformationGradientTensor[2][2] = 0.9;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeInverseMatrix(leftDeformationGradientTensor);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeInverseMatrix(rightDeformationGradientTensor);

    for (int i = 0; i < cellCount; i++)
    {
        material1Cells[i] = ElasticStateVector(2.0, 0.0, 0.1, leftDistortionTensor, 0.0, material1Parameters);
        material2Cells[i] = ElasticStateVector(0.0, -0.03, -0.01, rightDistortionTensor, 0.0, material2Parameters);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(ElasticRGFMSolver::solveTilde(ElasticMultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.05, 0.0, 0, 0, 0, material1Parameters,
                                                 material2Parameters));
}

void ElasticRGFMTests::solveZhangTestStar(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    vector<ElasticStateVector> material1Cells(cellCount);
    vector<ElasticStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    HyperelasticMaterialParameters material1Parameters(2.71, 6.22, 3.16, 9.0 * pow(10.0, -4.0), 300.0, 1.0, 3.577, 2.088);
    HyperelasticMaterialParameters material2Parameters(8.9, 4.6, 2.1, 3.9 * pow(10.0, -4.0), 300.0, 1.0, 3.0, 2.0);

    vector<vector<double> > leftDeformationGradientTensor(3, vector<double>(3));
    leftDeformationGradientTensor[0][0] = 1.0;
    leftDeformationGradientTensor[0][1] = 0.0;
    leftDeformationGradientTensor[0][2] = 0.0;
    leftDeformationGradientTensor[1][0] = -0.01;
    leftDeformationGradientTensor[1][1] = 0.95;
    leftDeformationGradientTensor[1][2] = 0.02;
    leftDeformationGradientTensor[2][0] = -0.015;
    leftDeformationGradientTensor[2][1] = 0.0;
    leftDeformationGradientTensor[2][2] = 0.9;

    vector<vector<double> > rightDeformationGradientTensor(3, vector<double>(3));
    rightDeformationGradientTensor[0][0] = 1.0;
    rightDeformationGradientTensor[0][1] = 0.0;
    rightDeformationGradientTensor[0][2] = 0.0;
    rightDeformationGradientTensor[1][0] = 0.015;
    rightDeformationGradientTensor[1][1] = 0.95;
    rightDeformationGradientTensor[1][2] = 0.0;
    rightDeformationGradientTensor[2][0] = -0.01;
    rightDeformationGradientTensor[2][1] = 0.0;
    rightDeformationGradientTensor[2][2] = 0.9;

    vector<vector<double> > leftDistortionTensor = MatrixAlgebra::computeInverseMatrix(leftDeformationGradientTensor);
    vector<vector<double> > rightDistortionTensor = MatrixAlgebra::computeInverseMatrix(rightDeformationGradientTensor);

    for (int i = 0; i < cellCount; i++)
    {
        material1Cells[i] = ElasticStateVector(2.0, 0.0, 0.1, leftDistortionTensor, 0.0, material1Parameters);
        material2Cells[i] = ElasticStateVector(0.0, -0.03, -0.01, rightDistortionTensor, 0.0, material2Parameters);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(ElasticRGFMSolver::solveStar(ElasticMultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.05, 0.0, 0, 0, 0, material1Parameters,
                                                material2Parameters));
}

void ElasticRGFMTests::outputSolution(ElasticMultimaterialSystem multimaterialSystem)
{
    vector<ElasticStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<ElasticStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    int cellCount = levelSetFunction.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream levelSetFunctionFile("levelSetFunction.dat");
    ofstream densityFile("density.dat");
    ofstream xVelocityFile("xVelocity.dat");
    ofstream yVelocityFile("yVelocity.dat");

    bool inGhostRegion = false;
    double levelSetValue = levelSetFunction[0];

    for (int i = 0; i < cellCount; i++)
    {
        if ((levelSetFunction[i] * levelSetValue) <= 0.0 && (levelSetFunction[i - 1] * levelSetValue) > 0.0)
        {
            inGhostRegion = !inGhostRegion;
        }

        if (!inGhostRegion)
        {
            densityFile << (cellSpacing * i) << " " << material1Cells[i].getDensity() << endl;
            xVelocityFile << (cellSpacing * i) << " " << material1Cells[i].getXVelocity() << endl;
            yVelocityFile << (cellSpacing * i) << " " << material1Cells[i].getYVelocity() << endl;
        }
        else
        {
            densityFile << (cellSpacing * i) << " " << material2Cells[i].getDensity() << endl;
            xVelocityFile << (cellSpacing * i) << " " << material2Cells[i].getXVelocity() << endl;
            yVelocityFile << (cellSpacing * i) << " " << material2Cells[i].getYVelocity() << endl;
        }

        levelSetFunctionFile << (cellSpacing * i) << " " << levelSetFunction[i] << endl;

        levelSetValue = levelSetFunction[i];
    }

    levelSetFunctionFile.close();
    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
}
