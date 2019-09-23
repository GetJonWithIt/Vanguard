#include "mhdrgfmtests.h"

MHDRGFMTests::MHDRGFMTests()
{
}

void MHDRGFMTests::solveDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    vector<MHDStateVector> material1Cells(cellCount);
    vector<MHDStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        material1Cells[i] = MHDStateVector(1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, 0.0);
        material2Cells[i] = MHDStateVector(0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(MHDRGFMSolver::solve(MHDMultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.1, 0.0, 0, 0, 0, materialParameters, materialParameters));
}

void MHDRGFMTests::solveDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    vector<MHDStateVector> material1Cells(cellCount);
    vector<MHDStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    MHDMaterialParameters material1Parameters(1.4);
    MHDMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        material1Cells[i] = MHDStateVector(1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, 0.0);
        material2Cells[i] = MHDStateVector(0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(MHDRGFMSolver::solve(MHDMultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.2, 0.0, 0, 0, 0, material1Parameters, material2Parameters));
}

void MHDRGFMTests::outputSolution(MHDMultimaterialSystem multimaterialSystem)
{
    vector<MHDStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<MHDStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    int cellCount = levelSetFunction.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream levelSetFunctionFile("levelSetFunction.dat");
    ofstream densityFile("density.dat");
    ofstream yMagneticFieldFile("yMagneticField.dat");

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
            yMagneticFieldFile << (cellSpacing * i) << " " << material1Cells[i].getYMagneticField() << endl;
        }
        else
        {
            densityFile << (cellSpacing * i) << " " << material2Cells[i].getDensity() << endl;
            yMagneticFieldFile << (cellSpacing * i) << " " << material2Cells[i].getYMagneticField() << endl;
        }

        levelSetFunctionFile << (cellSpacing * i) << " " << levelSetFunction[i] << endl;

        levelSetValue = levelSetFunction[i];
    }

    levelSetFunctionFile.close();
    densityFile.close();
    yMagneticFieldFile.close();
}
