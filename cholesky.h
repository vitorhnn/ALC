/**
 *  @file
 */
/*
 * Copyright (c) 2016 Victor Hermann "vitorhnn" Chiletto
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

#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "matrix.h"

/**
 *  Calcula o fator de Cholesky da matriz mat
 *  @return O fator de Cholesky de mat se esta é positiva definida, NULL caso contrário
 */
static matrix_t *cholesky_factor(matrix_t *mat)
{
    matrix_t *factor = matrix_new(mat->rows, mat->columns);

    size_t i;

    memset(factor->elements, 0, sizeof(double) * (factor->rows * factor->columns));

    for (i = 0; i < mat->rows; ++i) {
        size_t j, k;

        double current = matrix_get_at(mat, i, i);

        for (k = 0; k < i; ++k) {
            current -= pow(matrix_get_at(factor, k, i), 2);
        }

        if (current < 0) {
            matrix_free(factor);
            return NULL;
        }

        current = sqrt(current);

        matrix_set_at(factor, i, i, current);

        for (j = i + 1; j < mat->rows; ++j) {
            double current = matrix_get_at(mat, i, j);

            for (k = 0; k < i; ++k) {
                current -= matrix_get_at(factor, k, i) * matrix_get_at(factor, k, j);
            }

            current /= matrix_get_at(factor, i, i);

            matrix_set_at(factor, i, j, current);
        }
    }

    return factor;
}

/**
 *  Resolve o sistema Ax = b por fatoração de Cholesky
 *  @return O vetor x caso seja possível aplicar fatoração de Cholesky, NULL caso contrário.
 */
static matrix_t *cholesky_solve(matrix_t *A, matrix_t *b)
{
    matrix_t *factor = cholesky_factor(A);

    if (factor == NULL) {
        return NULL;
    }

    matrix_t *trans_factor = matrix_transpose(factor);

    matrix_t *y = forwards_substitution(trans_factor, b);
    
    matrix_t *x = backwards_substitution(factor, y);

    matrix_free(factor);
    matrix_free(trans_factor);
    matrix_free(y);

    return x;
}

#endif
