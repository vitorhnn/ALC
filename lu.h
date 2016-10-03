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

#ifndef LU_H
#define LU_H

#include "matrix.h"

/**
 *  Decompõe a matriz A em L, triangular inferior, e U, triangular superior.
 */
static void lu_decompose(matrix_t *A, matrix_t **L, matrix_t **U)
{
    /* doolittle reduction. */
    /* using the variables from watkin's book because this is confusing */
    size_t k;
    *L = matrix_new(A->rows, A->columns);
    matrix_load_identity(*L);
    *U = matrix_copy(*L);

    for (k = 0; k < A->rows; ++k) {
        size_t j, i;
        for (j = k; j < A->rows; ++j) {
            size_t m;
            double u = matrix_get_at(A, k, j);
            for (m = 0; m < k; ++m) {
                u -= matrix_get_at(*L, k, m) * matrix_get_at(*U, m, j);
            }
            matrix_set_at(*U, k, j, u);
        }

        for (i = k + 1; i < A->rows; ++i) {
            size_t m;
            double l = matrix_get_at(A, i, k);

            for (m = 0; m < k; ++m) {
                l -= matrix_get_at(*L, i, m) * matrix_get_at(*U, m, k);
            }

            l /= matrix_get_at(*U, k, k);

            matrix_set_at(*L, i, k, l);
        }
    }
}

/**
 *  Resolve o sistema Ax = b por decomposição LU
 *  @return O vetor x.
 */
static matrix_t *lu_solve(matrix_t *A, matrix_t *b)
{
    matrix_t *l, *u;
    lu_decompose(A, &l, &u);

    matrix_t *y = forwards_substitution(l, b);
    matrix_t *x = backwards_substitution(u, y);

    matrix_free(l);
    matrix_free(u);
    matrix_free(y);

    return x;
}

#endif
