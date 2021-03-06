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

void MHDRGFMTests::solve2DDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    vector<vector<MHDStateVector> > material1Cells(cellCount, vector<MHDStateVector>(cellCount));
    vector<vector<MHDStateVector> > material2Cells(cellCount, vector<MHDStateVector>(cellCount));
    vector<vector<double> > levelSetFunction(cellCount, vector<double>(cellCount));

    double cellSpacing = 1.0 / cellCount;
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            material1Cells[i][j] = MHDStateVector(1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, 0.0);
            material2Cells[i][j] = MHDStateVector(0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, 0.0);

            levelSetFunction[i][j] = (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) / cellCount) - 0.2;
        }
    }

    outputSolution2D(MHDRGFMSolver::solve2D(MHDMultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.05, 0.0, 0, subcyclingIterations, 0,
                                            materialParameters, materialParameters));
}

void MHDRGFMTests::solve2DDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    vector<vector<MHDStateVector> > material1Cells(cellCount, vector<MHDStateVector>(cellCount));
    vector<vector<MHDStateVector> > material2Cells(cellCount, vector<MHDStateVector>(cellCount));
    vector<vector<double> > levelSetFunction(cellCount, vector<double>(cellCount));

    double cellSpacing = 1.0 / cellCount;
    MHDMaterialParameters material1Parameters(1.4);
    MHDMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            material1Cells[i][j] = MHDStateVector(1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, 0.0);
            material2Cells[i][j] = MHDStateVector(0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, 0.0);

            levelSetFunction[i][j] = (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) / cellCount) - 0.2;
        }
    }

    outputSolution2D(MHDRGFMSolver::solve2D(MHDMultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.05, 0.0, 0, subcyclingIterations, 0,
                                            material1Parameters, material2Parameters));
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

void MHDRGFMTests::outputSolution2D(MHDMultimaterialSystem multimaterialSystem)
{
    vector<vector<MHDStateVector> > material1Cells = multimaterialSystem.getMaterial1Cells2D();
    vector<vector<MHDStateVector> > material2Cells = multimaterialSystem.getMaterial2Cells2D();
    vector<vector<double> > levelSetFunction = multimaterialSystem.getLevelSetFunction2D();

    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream levelSetFunctionFile("levelSetFunction.dat");
    ofstream densityFile("density.dat");
    ofstream yMagneticFieldFile("yMagneticField.dat");

    ofstream levelSetFunctionSliceFile("levelSetFunctionSlice.dat");
    ofstream densitySliceFile("densitySlice.dat");
    ofstream yMagneticFieldSliceFile("yMagneticFieldSlice.dat");

    bool inGhostRegion;
    double levelSetValue;

    for (int i = 0; i < rowCount; i++)
    {
        inGhostRegion = true;
        levelSetValue = levelSetFunction[i][0];

        for (int j = 0; j < columnCount; j++)
        {
            if ((levelSetFunction[i][j] * levelSetValue) <= 0.0 && (levelSetFunction[i][j - 1] * levelSetValue) > 0.0)
            {
                inGhostRegion = !inGhostRegion;
            }

            if (!inGhostRegion)
            {
                densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << material1Cells[i][j].getDensity() << endl;
                yMagneticFieldFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << material1Cells[i][j].getYMagneticField() << endl;

                if (i == floor(0.5 * rowCount))
                {
                    densitySliceFile << (cellSpacing * j) << " " << material1Cells[i][j].getDensity() << endl;
                    yMagneticFieldSliceFile << (cellSpacing * j) << " " << material1Cells[i][j].getYMagneticField() << endl;
                }
            }
            else
            {
                densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << material2Cells[i][j].getDensity() << endl;
                yMagneticFieldFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << material2Cells[i][j].getYMagneticField() << endl;

                if (i == floor(0.5 * rowCount))
                {
                    densitySliceFile << (cellSpacing * j) << " " << material2Cells[i][j].getDensity() << endl;
                    yMagneticFieldSliceFile << (cellSpacing * j) << " " << material2Cells[i][j].getYMagneticField() << endl;
                }
            }

            levelSetFunctionFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << levelSetFunction[i][j] << endl;

            if (i == floor(0.5 * rowCount))
            {
                levelSetFunctionSliceFile << (cellSpacing * j) << levelSetFunction[i][j] << endl;
            }

            levelSetValue = levelSetFunction[i][j];
        }
    }

    levelSetFunctionFile.close();
    densityFile.close();
    yMagneticFieldFile.close();

    levelSetFunctionSliceFile.close();
    densitySliceFile.close();
    yMagneticFieldSliceFile.close();
}
