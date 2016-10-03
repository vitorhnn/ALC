/**
 *  @file
 */
/*
 * Copyright (c) 2016 Andrei Parente
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
**/

#ifndef MATRIX_NORMS_H
#define MATRIX_NORMS_H

#include "matrix.h"

/**
 *  Calcula a norma de frobenius de mat
 *  @author Andrei Parente
 */
static double frobenius_norm(matrix_t *mat)
{
    double norm = 0;

    size_t i;
    for (i = 0; i < mat->rows; ++i) {
        size_t j;
        for (j = 0; j < mat->columns; ++j) {
            norm += pow(matrix_get_at(mat, i, j), 2);
        }
    }

    return sqrt(norm);
}

/**
 *  Calcula a norma linha de mat
 *  @author Andrei Parente
 */
static double row_norm(matrix_t *mat)
{
    double sum = 0,
           norm = 0;

    size_t i;
    for (i = 0; i < mat->rows; ++i) {
        size_t j;
        for (j = 0; j < mat->columns; ++j) {
            if (matrix_get_at(mat, i, j) >= 0) {
                sum += matrix_get_at(mat, i, j);
            } else {
                sum -= matrix_get_at(mat, i, j);
            }

            if (sum > norm) {
                norm = sum;
            }
        }

        sum = 0;
    }

    return norm;
}

/**
 *  Calcula a norma coluna de mat
 *  @author Andrei Parente
 */
static double column_norm(matrix_t *mat)
{
    double sum = 0,
           norm = 0;

    size_t j;
    for (j = 0; j < mat->columns; ++j) {
        size_t i;
        for (i = 0; i < mat->rows; ++i) {
            if (matrix_get_at(mat, i, j) >= 0) {
                sum += matrix_get_at(mat, i, j);
            } else {
                sum -= matrix_get_at(mat, i, j);
            }

            if (sum > norm) {
                norm = sum;
            }
        }

        sum = 0;
    }

    return norm;
}

#endif
