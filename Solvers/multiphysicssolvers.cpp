#include "multiphysicssolvers.h"

MultiphysicsSolvers::MultiphysicsSolvers()
{
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<EulerMultiphysicsStateVector> & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    int interfaceLocation = 0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        if (currentCells[i].getMaterial1VolumeFraction() < 0.5 && currentCells[i - 1].getMaterial1VolumeFraction() >= 0.5)
        {
            interfaceLocation = i;
        }
    }

    for (int i = interfaceLocation - 1; i < interfaceLocation + 2; i++)
    {
        currentCells[i].relaxTotalDensity();
        currentCells[i].relaxTotalPressure(material1Parameters, material2Parameters);

        if (i < interfaceLocation)
        {
            currentCells[i].setMaterial1VolumeFraction(0.999);
        }
        else
        {
            currentCells[i].setMaterial1VolumeFraction(0.001);
        }
    }
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<vector<EulerMultiphysicsStateVector> > & currentCells, EulerMaterialParameters material1Parameters,
                                                     EulerMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        int interfaceLocation1 = 0;
        int interfaceLocation2 = 0;

        for (int j = 1; j < columnCount; j++)
        {
            if (currentCells[i][j].getMaterial1VolumeFraction() < 0.5 && currentCells[i][j - 1].getMaterial1VolumeFraction() >= 0.5)
            {
                interfaceLocation1 = j;
            }
            if (currentCells[i][j].getMaterial1VolumeFraction() > 0.5 && currentCells[i][j - 1].getMaterial1VolumeFraction() <= 0.5)
            {
                interfaceLocation2 = j;
            }
        }

        if (interfaceLocation1 != 0)
        {
            for (int j = interfaceLocation1 - 1; j < interfaceLocation1 + 2; j++)
            {
                currentCells[i][j].relaxTotalDensity();
                currentCells[i][j].relaxTotalPressure(material1Parameters, material2Parameters);

                if (j < interfaceLocation1)
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.999);
                }
                else
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.001);
                }
            }
        }
        if (interfaceLocation2 != 0)
        {
            for (int j = interfaceLocation2 - 2; j < interfaceLocation2 + 1; j++)
            {
                currentCells[i][j].relaxTotalDensity();
                currentCells[i][j].relaxTotalPressure(material1Parameters, material2Parameters);

                if (j < interfaceLocation2)
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.001);
                }
                else
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.999);
                }
            }
        }
    }

    for (int i = 0; i < columnCount; i++)
    {
        int interfaceLocation1 = 0;
        int interfaceLocation2 = 0;

        for (int j = 1; j < rowCount; j++)
        {
            if (currentCells[j][i].getMaterial1VolumeFraction() < 0.5 && currentCells[j - 1][i].getMaterial1VolumeFraction() >= 0.5)
            {
                interfaceLocation1 = j;
            }
            if (currentCells[j][i].getMaterial1VolumeFraction() > 0.5 && currentCells[j - 1][i].getMaterial1VolumeFraction() <= 0.5)
            {
                interfaceLocation2 = j;
            }
        }

        if (interfaceLocation1 != 0)
        {
            for (int j = interfaceLocation1 - 1; j < interfaceLocation1 + 2; j++)
            {
                currentCells[j][i].relaxTotalDensity();
                currentCells[j][i].relaxTotalPressure(material1Parameters, material2Parameters);

                if (j < interfaceLocation1)
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.999);
                }
                else
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.001);
                }
            }
        }
        if (interfaceLocation2 != 0)
        {
            for (int j = interfaceLocation2 - 2; j < interfaceLocation2 + 1; j++)
            {
                currentCells[j][i].relaxTotalDensity();
                currentCells[j][i].relaxTotalPressure(material1Parameters, material2Parameters);

                if (j < interfaceLocation2)
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.001);
                }
                else
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.999);
                }
            }
        }
    }
}
