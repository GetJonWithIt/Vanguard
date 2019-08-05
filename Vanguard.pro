QT += core
QT -= gui

TARGET = Vanguard
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

SOURCES += main.cpp \
    Solvers/firstordersolver.cpp \
    Solvers/secondordersolver.cpp \
    Mathematics/vectoralgebra.cpp \
    Solvers/solvers.cpp \
    Mathematics/matrixalgebra.cpp \
    Mathematics/tensoralgebra.cpp \
    Elasticity/elasticequationofstate.cpp \
    Elasticity/elasticstatevector.cpp \
    Elasticity/elasticacoustictensor.cpp \
    Elasticity/hyperelasticmaterialparameters.cpp \
    Solvers/slopelimiters.cpp \
    Elasticity/elastictests.cpp \
    Elasticity/Multiphysics/elasticmultiphysicsstatevector.cpp \
    Euler/eulermaterialparameters.cpp \
    Euler/eulerequationofstate.cpp \
    Euler/eulerstatevector.cpp \
    Euler/eulertests.cpp \
    Euler/Multiphysics/eulermultiphysicsstatevector.cpp \
    Euler/Multiphysics/eulermultiphysicstests.cpp \
    Solvers/multiphysicssolvers.cpp \
    MHD/mhdstatevector.cpp \
    MHD/mhdmaterialparameters.cpp \
    Elasticity/Multiphysics/elasticmultiphysicstests.cpp \
    Elasticity/Multiphysics/elasticreducedstatevector.cpp \
    Elasticity/Multiphysics/elasticreducedtests.cpp \
    Euler/Multiphysics/eulerreducedstatevector.cpp \
    Euler/Multiphysics/eulerreducedtests.cpp \
    Solvers/elasticsolvers.cpp \
    Solvers/elasticfirstordersolver.cpp \
    Solvers/elasticsecondordersolver.cpp \
    MHD/mhdequationofstate.cpp \
    MHD/mhdwavespeeds.cpp \
    Solvers/mhdsolvers.cpp \
    Solvers/mhdfirstordersolver.cpp \
    MHD/mhdtests.cpp \
    Solvers/mhdsecondordersolver.cpp \
    HPR/hprmaterialparameters.cpp \
    HPR/hprmiegruneisen.cpp \
    HPR/hprderivatives.cpp \
    HPR/hprequationofstate.cpp \
    HPR/hprsourceterms.cpp \
    HPR/hprstatevector.cpp \
    HPR/hprwavespeeds.cpp \
    HPR/hpracoustictensor.cpp \
    Solvers/hprsolvers.cpp \
    Solvers/hprfirstordersolver.cpp \
    HPR/hprtests.cpp \
    Solvers/hprsecondordersolver.cpp \
    Solvers/hprforcingsolver.cpp \
    HPR/Multiphysics/hprreducedstatevector.cpp \
    HPR/Multiphysics/hprreducedacoustictensor.cpp \
    HPR/Multiphysics/hprreducedtests.cpp \
    HPR/Multiphysics/hprintermediatestatevector.cpp \
    HPR/Multiphysics/hprintermediateacoustictensor.cpp \
    HPR/Multiphysics/hprintermediatetests.cpp \
    HPR/Multiphysics/hprmultiphysicsstatevector.cpp \
    HPR/Multiphysics/hprmultiphysicsacoustictensor.cpp \
    HPR/Multiphysics/hprmultiphysicstests.cpp \
    Solvers/mhdforcingsolver.cpp \
    MHD/Multiphysics/mhdreducedstatevector.cpp \
    MHD/Multiphysics/mhdreducedtests.cpp \
    MHD/Multiphysics/mhdintermediatestatevector.cpp \
    MHD/Multiphysics/mhdintermediatetests.cpp \
    MHD/Multiphysics/mhdmultiphysicsstatevector.cpp \
    MHD/Multiphysics/mhdmultiphysicstests.cpp \
    RGFM/exactsolver.cpp \
    RGFM/multimaterialsystem.cpp \
    RGFM/rgfmsolver.cpp \
    RGFM/eulerrgfmtests.cpp \
    RGFM/hllcsolver.cpp

HEADERS += \
    Solvers/firstordersolver.h \
    Solvers/secondordersolver.h \
    Mathematics/vectoralgebra.h \
    Solvers/solvers.h \
    Mathematics/matrixalgebra.h \
    Mathematics/tensoralgebra.h \
    Elasticity/elasticequationofstate.h \
    Elasticity/elasticstatevector.h \
    Elasticity/elasticacoustictensor.h \
    Elasticity/hyperelasticmaterialparameters.h \
    Solvers/slopelimiters.h \
    Elasticity/elastictests.h \
    Elasticity/Multiphysics/elasticmultiphysicsstatevector.h \
    Euler/eulermaterialparameters.h \
    Euler/eulerequationofstate.h \
    Euler/eulerstatevector.h \
    Euler/eulertests.h \
    Euler/Multiphysics/eulermultiphysicsstatevector.h \
    Euler/Multiphysics/eulermultiphysicstests.h \
    Solvers/multiphysicssolvers.h \
    MHD/mhdstatevector.h \
    MHD/mhdmaterialparameters.h \
    Elasticity/Multiphysics/elasticmultiphysicstests.h \
    Elasticity/Multiphysics/elasticreducedstatevector.h \
    Elasticity/Multiphysics/elasticreducedtests.h \
    Euler/Multiphysics/eulerreducedstatevector.h \
    Euler/Multiphysics/eulerreducedtests.h \
    Solvers/elasticsolvers.h \
    Solvers/elasticfirstordersolver.h \
    Solvers/elasticsecondordersolver.h \
    MHD/mhdequationofstate.h \
    MHD/mhdwavespeeds.h \
    Solvers/mhdsolvers.h \
    Solvers/mhdfirstordersolver.h \
    MHD/mhdtests.h \
    Solvers/mhdsecondordersolver.h \
    HPR/hprmaterialparameters.h \
    HPR/hprmiegruneisen.h \
    HPR/hprderivatives.h \
    HPR/hprequationofstate.h \
    HPR/hprsourceterms.h \
    HPR/hprstatevector.h \
    HPR/hprwavespeeds.h \
    HPR/hpracoustictensor.h \
    Solvers/hprsolvers.h \
    Solvers/hprfirstordersolver.h \
    HPR/hprtests.h \
    Solvers/hprsecondordersolver.h \
    Solvers/hprforcingsolver.h \
    HPR/Multiphysics/hprreducedstatevector.h \
    HPR/Multiphysics/hprreducedacoustictensor.h \
    HPR/Multiphysics/hprreducedtests.h \
    HPR/Multiphysics/hprintermediatestatevector.h \
    HPR/Multiphysics/hprintermediateacoustictensor.h \
    HPR/Multiphysics/hprintermediatetests.h \
    HPR/Multiphysics/hprmultiphysicsstatevector.h \
    HPR/Multiphysics/hprmultiphysicsacoustictensor.h \
    HPR/Multiphysics/hprmultiphysicstests.h \
    Solvers/mhdforcingsolver.h \
    MHD/Multiphysics/mhdreducedstatevector.h \
    MHD/Multiphysics/mhdreducedtests.h \
    MHD/Multiphysics/mhdintermediatestatevector.h \
    MHD/Multiphysics/mhdintermediatetests.h \
    MHD/Multiphysics/mhdmultiphysicsstatevector.h \
    MHD/Multiphysics/mhdmultiphysicstests.h \
    RGFM/exactsolver.h \
    RGFM/multimaterialsystem.h \
    RGFM/rgfmsolver.h \
    RGFM/eulerrgfmtests.h \
    RGFM/hllcsolver.h

