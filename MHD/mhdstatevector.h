#ifndef MHDSTATEVECTOR_H
#define MHDSTATEVECTOR_H


class MHDStateVector
{
public:
    MHDStateVector();

    void setDensity(double newDensity);
    void setXVelocity(double newXVelocity);
    void setYVelocity(double newYVelocity);
    void setZVelocity(double newZVelocity);
    void setPressure(double newPressure);

    void setXMagneticField(double newXMagneticField);
    void setYMagneticField(double newYMagneticField);
    void setZMagneticField(double newZMagneticField);

    void setAuxiliaryField(double newAuxiliaryField);

    double getDensity();
    double getXVelocity();
    double getYVelocity();
    double getZVelocity();
    double getPressure();

    double getXMagneticField();
    double getYMagneticField();
    double getZMagneticField();

    double getAuxiliaryField();

private:
    double density;
    double xVelocity;
    double yVelocity;
    double zVelocity;
    double pressure;

    double xMagneticField;
    double yMagneticField;
    double zMagneticField;

    double auxiliaryField;
};

#endif // MHDSTATEVECTOR_H
