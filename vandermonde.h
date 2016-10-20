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

#ifndef VANDERMONDE_H
#define VANDERMONDE_H

#include "matrix.h"

/** 
 *  Verifica se a matriz é do tipo Vandermonde
 *  @author Andrei Parente
 */
static int vandermonde_check(const matrix_t *mat)
{
    size_t i;
    for (i = 0; i < mat->rows; ++i) {
        double first = 0;

        size_t j;
        for (j = 0; j < mat->columns; ++j) {
            if (j == 0 && matrix_get_at(mat, i, j) != 1) {
                return 0;
            } else if (j == 1) {
                first = matrix_get_at(mat, i, j);
            } else if (matrix_get_at(mat, i, j) != pow(first, j)) {
                return 0;
            }
        }
    }

    return 1;
}

/**
 *  Calcula o determinante de mat, usando o processo simplificado para matrizes de Vandermonde
 *  @author Andrei Parente
 */
static double vandermonde_determinant(const matrix_t *mat)
{
    double determinant = 1;

    size_t i;

    assert(vandermonde_check(mat) == 1);

    for (i = 0; i < (mat->rows - 1); ++i) {
        size_t j;
        for (j = (mat->rows - 1); j > i; --j) {
            determinant *= (matrix_get_at(mat, j, 1) - matrix_get_at(mat, i, 1));
        } 
    }

    return determinant;
}

#endif
