/**
 * @file BatchMandelCalculator.h
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */
#ifndef BATCHMANDELCALCULATOR_H
#define BATCHMANDELCALCULATOR_H

#include <BaseMandelCalculator.h>

class BatchMandelCalculator : public BaseMandelCalculator
{
public:
    BatchMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~BatchMandelCalculator();
    int * calculateMandelbrot();

private:
    // buffer for output data (iteration count for each point)
    int *data;
    // real and image values for each point in a row
    float *imag;
    float *real;
    // indicates if point in row went over the limiting z_n absolute value
    bool *hitOverMax;
    
    const int batchBasicSize = 128;

    // real values of complex numbers, precounted for speedup
    float *realPrecounted;
};

#endif