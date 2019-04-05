#include <QCoreApplication>

#include "Euler/eulertests.h"
#include "Euler/Multiphysics/eulermultiphysicstests.h"
#include "Elasticity/elastictests.h"
#include "Elasticity/Multiphysics/elasticmultiphysicstests.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    //EulerTests::solveToroTest1(400);
    //EulerTests::solve2DToroTest1(100);
    //EulerMultiphysicsTests::solveToroTest1(100, 3);
    //EulerMultiphysicsTests::solveFedkiwTest(100, 3);
    //EulerMultiphysicsTests::solveChinnayyaTest(100);

    //EulerMultiphysicsTests::solve2DToroTest1(100, 3);
    //EulerMultiphysicsTests::solve2DFedkiwTest(200, 3);

    //ElasticTests::solveBartonTest1(500);
    //ElasticTests::solveBartonTest2(200);
    //ElasticTests::solve2DBartonTest1(50);
    //ElasticTests::solve2DBartonTest2(20);

    //ElasticMultiphysicsTests::solveBartonTest1(100, 3);
    //ElasticMultiphysicsTests::solveBartonTest2(100, 3);
    ElasticMultiphysicsTests::solveZhangTest(100, 0);

    return a.exec();
}
