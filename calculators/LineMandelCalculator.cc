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
	// data = (int *)(malloc(height * width * sizeof(int)));

	
	
	data = (int *)(aligned_alloc(ALIGNMENT, height * width * sizeof(int)));
	// std::fill(data, data + (height * width), limit);
	
	imag = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
	real = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
	// hitOverMax = (bool *)(aligned_alloc(ALIGNMENT, width * sizeof(bool)));
}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	data = NULL;

	free(imag);
	imag = NULL;

	free(real);	
	real = NULL;
	
	// free(hitOverMax);	
	// hitOverMax = NULL;
}


int * LineMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *zReal = this->real;
	float *zImag = this->imag;
	// bool *phitOverMax = this->hitOverMax;
	
	#pragma omp simd aligned(pdata : ALIGNMENT)
	for (int i = 0; i < width * height; i++)
		pdata[i] = limit;

	for (int row = 0; row < height / 2; row++)
	{
		const float rowImag = y_start + row * dy; // current imaginary value for this row

		// set init z_n values of all points in row
		#pragma omp simd aligned(zReal : ALIGNMENT)
		for (int col = 0; col < width; col++)
		{
			zReal[col] = x_start + col * dx;
		}
		
		std::fill(zImag, zImag + width, rowImag);
		// std::fill(hitOverMax, hitOverMax + width, false); // reset hit buffer

		bool allFinished = false;
		// one iteration over entire row
		for (int i = 0; i < limit; ++i)
		{
			// check if all points finished
			// #pragma omp parallel for reduction(&7 : allFinished)
			// for (int col = 0; col < width; col++)
			// {
			// 	allFinished &= phitOverMax[col];
			// }
			// if (allFinished) // stop row calculation if they do
			// {
			// 	break;
			// 	std::cout << "test";
			// }
			
			// for each column in row
			// #pragma omp simd 
			// #pragma omp simd aligned(data, zReal, zImag:ALIGNMENT)
			#pragma omp simd aligned(pdata, zReal, zImag : ALIGNMENT)

			for (int col = 0; col < width; col++)
			{
				const int pdIdx = row * width + col; // index into data for this point
				if (pdata[pdIdx] == limit)
				// if (phitOverMax[pdIdx] == false)
				{
					float x = x_start + col * dx;
					float r2 = zReal[col] * zReal[col];
					float i2 = zImag[col] * zImag[col];

					if (r2 + i2 > 4.0f)
					{
						pdata[pdIdx] = i;
						// phitOverMax[col] = true;
					}
					zImag[col] = 2.0f * zReal[col] * zImag[col] + rowImag;
					zReal[col] = r2 - i2 + x;
				}
			}
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


	// for (int i = 0; i < height; i++)
	// {
	// 	for (int j = 0; j < width; j++)
	// 	{
	// 		float x = x_start + j * dx; // current real value
	// 		float y = y_start + i * dy; // current imaginary value
	// 		// int value = mandelbrot(x, y, limit);
	// 		// *(pdata++) = value;

			
	// 		float zReal = x;
	// 		float zImag = y;

	// 		int iterCnt = limit;
	// 		for (int l = 0; l < limit; ++l)
	// 		{
	// 			float r2 = zReal * zReal;
	// 			float i2 = zImag * zImag;

	// 			if (r2 + i2 > 4.0f)
	// 			{
	// 				iterCnt = l;
	// 				break;
	// 			}

	// 			zImag = 2.0f * zReal * zImag + y;
	// 			zReal = r2 - i2 + x;
	// 		}
	// 		*(pdata++) = iterCnt;			
	// 	}
	// }
	return data;
}
