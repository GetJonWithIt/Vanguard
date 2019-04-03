#include <QCoreApplication>

#include "Euler/eulertests.h"
#include "Euler/Multiphysics/eulermultiphysicstests.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    //EulerTests::solveToroTest1(400);
    //EulerTests::solve2DToroTest1(100);
    //EulerMultiphysicsTests::solveToroTest1(100, 3);
    //EulerMultiphysicsTests::solveFedkiwTest(100, 3);
    //EulerMultiphysicsTests::solveChinnayyaTest(100);

    //EulerMultiphysicsTests::solve2DToroTest1(100, 0);
    EulerMultiphysicsTests::solve2DFedkiwTest(100, 3);

    return a.exec();
}
