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
const size_t ALIGNMENT = 64;

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
	data = (int *)(aligned_alloc(ALIGNMENT, height * width * sizeof(int)));
	#pragma omp simd aligned(data : ALIGNMENT)
	for (int i = 0; i < width * height; i++)
		data[i] = limit;
	
	imag = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
	real = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
	hitOverMax = (bool *)(aligned_alloc(ALIGNMENT, width * sizeof(bool)));

	realPrecounted = (float *)(aligned_alloc(ALIGNMENT, width * sizeof(float)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
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


int * BatchMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *zReal = this->real;
	float *zImag = this->imag;
	bool *phitOverMax = this->hitOverMax;
	float *realPre = this->realPrecounted;
	
	const int batchSize = (this->batchBasicSize < width) ? this->batchBasicSize : width;

	for (int row = 0; row < height / 2; row++)
	{
		const float rowImag = y_start + row * dy; // current imaginary value for this row
		std::fill(zImag, zImag + width, rowImag);
		std::fill(hitOverMax, hitOverMax + width, false); // reset hit buffer

		// iterate over batches
		for (int b = 0; b < width; b += batchSize)
		{
			// set init z_n real values of all points in batch
			#pragma omp simd aligned(zReal, realPre : ALIGNMENT)
			for (int col = b; col < b + batchSize; col++)
			{
				float x = x_start + col * dx;
				zReal[col] = x;
				realPre[col] = x;
			}

			bool batchFinished = false;
			// one iteration over batch
			for (int i = 0; i < limit; ++i)
			{
				int allFinished = 0;
				#pragma omp simd aligned(pdata, zReal, zImag : ALIGNMENT) reduction(+ : allFinished)
				for (int col = b; col < b + batchSize; col++)
				{					
					const int pdIdx = row * width + col; // index into data for this point
					if (phitOverMax[col] == false)
					{
						float x = realPre[col];
						float r2 = zReal[col] * zReal[col];
						float i2 = zImag[col] * zImag[col];

						if (r2 + i2 > 4.0f)
						{
							pdata[pdIdx] = i;
							phitOverMax[col] = true;
						}
						zImag[col] = 2.0f * zReal[col] * zImag[col] + rowImag;
						zReal[col] = r2 - i2 + x;
					}
					else
					{
						allFinished += 1;
					}
				
				}
				if (allFinished == batchSize)
				{
					batchFinished = true;
					break;
				}
			}
			if (batchFinished)
			{
				continue;
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