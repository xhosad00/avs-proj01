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
#include <stdio.h>


#include "LineMandelCalculator.h"
const size_t ALIGNMENT = 64;

LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
	// allocate buffers (aligned)
	data = (int *)(aligned_alloc(ALIGNMENT, height * width * sizeof(int)));
	#pragma omp simd aligned(data : ALIGNMENT)
	for (int i = 0; i < width * height; i++)
		data[i] = limit;
	
	imag = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
	real = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
	hitOverMax = (bool *)(aligned_alloc(ALIGNMENT, width * sizeof(bool)));
	realPrecounted = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	data = NULL;

	free(imag);
	imag = NULL;

	free(real);	
	real = NULL;
	
	free(hitOverMax);	
	hitOverMax = NULL;

	free(realPrecounted);	
	realPrecounted = NULL;
}


int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *zReal = this->real;
	float *zImag = this->imag;
	bool *phitOverMax = this->hitOverMax;
	float *realPre = this->realPrecounted;	

	// calculate only half of height, the rest can be copied due to symmetry
	for (int row = 0; row < height / 2; row++)
	{
		const float rowImag = y_start + row * dy; // current imaginary value for this row

		// set init z_n values of all points in row
		#pragma omp simd aligned(zReal : ALIGNMENT)
		for (int col = 0; col < width; col++)
		{
			float x = x_start + col * dx;
			zReal[col] = x;
			realPre[col] = x;
		}
		
		std::fill(zImag, zImag + width, rowImag);
		std::fill(hitOverMax, hitOverMax + width, false); // reset hit buffer

		bool rowFinished = false;

		// one iteration over entire row
		for (int i = 0; i < limit; ++i)
		{
			int allFinished = 0; // variable for reduction
			// for each column in row
			#pragma omp simd aligned(pdata, zReal, zImag, phitOverMax, realPre : ALIGNMENT) reduction(+ : allFinished) // reduction(&& : allFinished) reduction did not compile correctly
			for (int col = 0; col < width; col++)
			{
				const int pdIdx = row * width + col; // index into data for this point
				if (phitOverMax[col] == false)
				{
					float zR = zReal[col];
					float zI = zImag[col];
					float x = realPre[col];
					float r2 = zR * zR;
					float i2 = zI * zI;

					if (r2 + i2 > 4.0f) // point will go to infinity
					{
						pdata[pdIdx] = i;
						phitOverMax[col] = true;
					}
					zImag[col] = 2.0f * zR * zI + rowImag;
					zReal[col] = r2 - i2 + x;
				}
				else
				{
					allFinished += 1;
				}
			}
			
			if (allFinished == width) // all points in row hit limit, break of of this row
			{
				rowFinished = true;
				break;
			}
		}
		if (rowFinished)
		{
			continue;
		}
	}

	const int halfHeight = height / 2;
	for (int row = halfHeight; row < height; row++)
	{
		#pragma omp simd aligned(pdata : ALIGNMENT)
		for (int col = 0; col < width; col++)
		{
			pdata[row * width + col] = pdata[(height - row - 1) * width + col];
		}
	}
	return data;
}
