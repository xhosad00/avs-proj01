/**
 * @file LineMandelCalculator.h
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */

#include <BaseMandelCalculator.h>

class LineMandelCalculator : public BaseMandelCalculator
{
public:
    LineMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~LineMandelCalculator();
    int *calculateMandelbrot();

private:
    // buffer for output data (iteration count for each point)
    int *data;
    // real and image values for each point in a row
    float *imag;
    float *real;
    // indicates if point in row went over the limiting z_n absolute value
    bool *hitOverMax;

    // real values of complex numbers, precounted for speedup
    float *realPrecounted;
};