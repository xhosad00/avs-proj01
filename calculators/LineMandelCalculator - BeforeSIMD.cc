/**
 * @file LineMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	// data = (int *)(malloc(height * width * sizeof(int)));

	const size_t alignment = 16;
	// buffer for output data (iteration count for each point)
	data = (int *)(aligned_alloc(alignment, height * width * sizeof(int)));
	// std::fill(data, data + (height * width), limit);

	// real and image values for each point in a row
	imag = (float *)(aligned_alloc(alignment, width * sizeof(float)));
	real = (float *)(aligned_alloc(alignment, width * sizeof(float)));
}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	data = NULL;
}

template <typename T>
static inline int mandelbrot(T real, T imag, int limit)
{
	T zReal = real;
	T zImag = imag;

	for (int i = 0; i < limit; ++i)
	{
		T r2 = zReal * zReal;
		T i2 = zImag * zImag;

		if (r2 + i2 > 4.0f)
			return i;

		zImag = 2.0f * zReal * zImag + imag;
		zReal = r2 - i2 + real;
	}
	return limit;
}

int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			float x = x_start + j * dx; // current real value
			float y = y_start + i * dy; // current imaginary value


			// int value = mandelbrot(x, y, limit);
			// *(pdata++) = value;

			
			float zReal = x;
			float zImag = y;

			int iterCnt = 100;
			for (int l = 0; l < limit; ++l)
			{
				float r2 = zReal * zReal;
				float i2 = zImag * zImag;

				if (r2 + i2 > 4.0f)
				{
					iterCnt = l;
					break;
				}

				zImag = 2.0f * zReal * zImag + y;
				zReal = r2 - i2 + x;
			}
			*(pdata++) = iterCnt;			
		}
	}
	return data;
}
