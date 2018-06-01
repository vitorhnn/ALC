
/**
 *  @file
 */
/*
 * Copyright (c) 2018 Victor Hermann "vitorhnn" Chiletto
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

#ifndef GAUSS_H
#define GAUSS_H

#include "matrix.h"
#include "basic.h"

static void swap_rows(matrix_t *A, size_t m, size_t k)
{
    size_t n = A->columns;
    size_t i = 0;

    for (i = 0; i < n; ++i) {
        double swap = matrix_get_at(A, m, i);
        matrix_set_at(A, m, i, matrix_get_at(A, k, i));
        matrix_set_at(A, k, i, swap);
    }
}

static void gauss_eliminate(matrix_t *A, matrix_t *b, int *singular)
{
    size_t n;
    size_t k;

    n = A->rows;

    for (k = 0; k < n; ++k) {
        size_t i;
        size_t max = 0;
        double maxval = 0.0;
        for (i = k; i < n; ++i) {
            double Aik = matrix_get_at(A, i, k);
            if (fabs(Aik) > maxval) {
                max = i;
                maxval = fabs(Aik);
            }
        }

        if (maxval == 0.0) {
            *singular = 1;
            return;
        }

        if (max != k) {
            swap_rows(A, max, k);
            swap_rows(b, max, k);
        }

        for (i = k + 1; i < n; ++i) {
            double mult = matrix_get_at(A, i, k) / matrix_get_at(A, k, k);

            size_t j;
            for (j = k + 1; j < n; ++j) {
                double thing = matrix_get_at(A, i, j) - mult * matrix_get_at(A, k, j);
                matrix_set_at(A, i, j, thing);
            }

            double thing = matrix_get_at(b, i, 0) - mult * matrix_get_at(b, k, 0);
            matrix_set_at(b, i, 0, thing);
        }
    }
}

static matrix_t *gauss_solve(matrix_t *A, matrix_t *b)
{
    int singular = 0;
    gauss_eliminate(A, b, &singular);

    if (singular) {
        return NULL;
    }

    matrix_t *x = backwards_substitution(A, b);

    return x;
}

#endif
