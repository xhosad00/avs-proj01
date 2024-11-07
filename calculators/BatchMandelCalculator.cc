/**
 * @file BatchMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data = (int *)(aligned_alloc(ALIGNMENT, height * width * sizeof(int)));
	
	imag = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
	real = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	free(data);
	data = NULL;

	free(imag);
	imag = NULL;

	free(real);	
	real = NULL;
}


int * BatchMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *zReal = this->real;
	float *zImag = this->imag;
	
	const int batchSize = (this->batchBasicSize < width) ? batchSize : width;

	
	#pragma omp simd aligned(pdata : ALIGNMENT)
	for (int i = 0; i < width * height; i++)
		pdata[i] = limit;

	for (int row = 0; row < height / 2; row++)
	{
		const float rowImag = y_start + row * dy; // current imaginary value for this row
		std::fill(zImag, zImag + width, rowImag);

		// iterate over batches
		for (int b = 0; b < width; b += batchSize)
		{
			// set init z_n values of all points in batch
			#pragma omp simd aligned(zReal : ALIGNMENT)
			for (int col = b; col < b + batchSize; col++)
			{
				zReal[col] = x_start + col * dx;
			}

			// one iteration over batch
			for (int i = 0; i < limit; ++i)
			{
				#pragma omp simd aligned(pdata, zReal, zImag : ALIGNMENT)
				for (int col = b; col < b + batchSize; col++)
				{					
					const int pdIdx = row * width + col; // index into data for this point
					if (pdata[pdIdx] == limit)
					{
						float x = x_start + col * dx;
						float r2 = zReal[col] * zReal[col];
						float i2 = zImag[col] * zImag[col];

						if (r2 + i2 > 4.0f)
						{
							pdata[pdIdx] = i;
						}
						zImag[col] = 2.0f * zReal[col] * zImag[col] + rowImag;
						zReal[col] = r2 - i2 + x;
					}
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
	return data;
}