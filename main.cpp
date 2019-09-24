#include <QCoreApplication>

#include "Euler/eulertests.h"
#include "Euler/Multiphysics/eulermultiphysicstests.h"
#include "Euler/Multiphysics/eulerreducedtests.h"
#include "Elasticity/elastictests.h"
#include "Elasticity/Multiphysics/elasticmultiphysicstests.h"
#include "Elasticity/Multiphysics/elasticreducedtests.h"
#include "MHD/mhdtests.h"
#include "MHD/Multiphysics/mhdmultiphysicstests.h"
#include "MHD/Multiphysics/mhdintermediatetests.h"
#include "MHD/Multiphysics/mhdreducedtests.h"
#include "HPR/hprtests.h"
#include "HPR/Multiphysics/hprmultiphysicstests.h"
#include "HPR/Multiphysics/hprintermediatetests.h"
#include "HPR/Multiphysics/hprreducedtests.h"

#include "RGFM/eulerrgfmtests.h"
#include "RGFM/mhdrgfmtests.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    //EulerTests::solveToroTest1(1000);
    //EulerTests::solve2DToroTest1(100);
    //EulerMultiphysicsTests::solveToroTest1(400, 0);
    //EulerMultiphysicsTests::solveFedkiwTest(400, 0);
    //EulerMultiphysicsTests::solveChinnayyaTest(100);

    //EulerMultiphysicsTests::solve2DToroTest1(100, 3);
    //EulerMultiphysicsTests::solve2DFedkiwTest(100, 3);

    //EulerReducedTests::solveToroTest1(100, 3);
    //EulerReducedTests::solveFedkiwTest(200, 0);
    //EulerReducedTests::solve2DToroTest1(100, 3);
    //EulerReducedTests::solve2DFedkiwTest(100, 3);

    //ElasticTests::solveBartonTest1(100);
    //ElasticTests::solveBartonTest2(200);
    //ElasticTests::solve2DBartonTest1(100);
    //ElasticTests::solve2DBartonTest2(20);

    //ElasticMultiphysicsTests::solveBartonTest1(100, 3);
    //ElasticMultiphysicsTests::solveBartonTest2(100, 3);
    //ElasticMultiphysicsTests::solveZhangTest(200, 0);

    //ElasticReducedTests::solveBartonTest1(100, 3);
    //ElasticReducedTests::solveZhangTest(800, 3);
    //ElasticReducedTests::solve2DBartonTest1(100, 3);
    //ElasticReducedTests::solve2DZhangTest(100, 3);
    //ElasticReducedTests::solve2DZhangTest2(100, 3);

    //MHDTests::solveDumbserTest2(100);
    //MHDTests::solveDumbserTest1(400, 1);
    //MHDTests::solveDumbserTest3(100);
    //MHDTests::solve2DDumbserTest1(100, 1);
    //MHDTests::solve2DDumbserTest2(100, 1);

    //MHDMultiphysicsTests::solveDumbserTest1(100, 0, 3);
    //MHDMultiphysicsTests::solveDumbserMultimaterialTest1(400, 0, 0);

    //MHDIntermediateTests::solveDumbserTest1(800, 1, 3);
    //MHDIntermediateTests::solveDumbserMultimaterialTest1(800, 1, 3);
    //MHDIntermediateTests::solve2DDumbserTest1(100, 1, 3);
    //MHDIntermediateTests::solve2DDumbserMultimaterialTest1(100, 1, 3);

    //MHDReducedTests::solveDumbserTest1(400, 1, 3);
    //MHDReducedTests::solveDumbserMultimaterialTest1(400, 1, 3);
    //MHDReducedTests::solve2DDumbserTest1(100, 1, 0);
    //MHDReducedTests::solve2DDumbserMultimaterialTest1(100, 1, 3);

    //HPRTests::solveStokesFirstProblem(800, 100);
    //HPRTests::solveHeatConductionProblem(800, 20);
    //HPRTests::solveBartonTest1(100, 0);
    //HPRTests::solveZhangElasticPlasticTest(800, 1);
    //HPRTests::solve2DBartonTest1(100, 0);
    //HPRTests::solve2DZhangElasticPlasticTest(100, 1);

    //HPRMultiphysicsTests::solveBartonTest1(100, 3);
    //HPRMultiphysicsTests::solveZhangTest(100, 0);

    //HPRIntermediateTests::solveBartonTest1(100, 3);
    //HPRIntermediateTests::solveZhangTest(100, 3);
    //HPRIntermediateTests::solveZhangElasticPlasticTest(100, 1, 3);
    //HPRIntermediateTests::solve2DBartonTest1(100, 3);
    //HPRIntermediateTests::solve2DZhangTest(100, 3);
    //HPRIntermediateTests::solve2DZhangTest2(100, 3);
    //HPRIntermediateTests::solve2DZhangElasticPlasticTest(100, 1, 3);
    //HPRIntermediateTests::solve2DZhangElasticPlasticTest2(100, 1, 3);
    //HPRIntermediateTests::solve2DUdayKumarTest(100, 1, 3);

    //HPRReducedTests::solveBartonTest1(800, 3);
    //HPRReducedTests::solveZhangTest(800, 3);
    //HPRReducedTests::solveZhangElasticPlasticTest(800, 1, 3);
    //HPRReducedTests::solve2DBartonTest1(100, 3);
    //HPRReducedTests::solve2DZhangTest(100, 3);
    //HPRReducedTests::solve2DZhangTest2(100, 3);
    //HPRReducedTests::solve2DZhangElasticPlasticTest(100, 1, 3);
    //HPRReducedTests::solve2DZhangElasticPlasticTest2(100, 1, 3);
    //HPRReducedTests::solve2DUdayKumarTest(100, 1, 3);

    //EulerRGFMTests::solveToroTest1Exact(400, 0);
    //EulerRGFMTests::solveToroTest1HLLC(100, 0);
    //EulerRGFMTests::solveFedkiwTestExact(400, 0);
    //EulerRGFMTests::solveFedkiwTestHLLC(400, 0);

    //EulerRGFMTests::solve2DToroTest1Exact(100, 0);
    //EulerRGFMTests::solve2DToroTest1HLLC(200, 0);
    //EulerRGFMTests::solve2DFedkiwTestExact(100, 0);
    //EulerRGFMTests::solve2DFedkiwTestHLLC(100, 0);

    //MHDRGFMTests::solveDumbserTest1(400, 0, 0);
    //MHDRGFMTests::solveDumbserMultimaterialTest1(400, 0, 0);

    //MHDRGFMTests::solve2DDumbserTest1(100, 1, 0);
    MHDRGFMTests::solve2DDumbserMultimaterialTest1(100, 1, 0);

    return a.exec();
}
