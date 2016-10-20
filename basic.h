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

#ifndef BASIC_OPS_H
#define BASIC_OPS_H

#include "matrix.h"

/** 
 *  Resolve o sistema Ax = b por substituição para frente.
 *  @return O vetor x.
 */
static matrix_t *forwards_substitution(const matrix_t * restrict A, const matrix_t * restrict b)
{
    size_t i;
    matrix_t *x = matrix_new(b->rows, 1);

    for (i = 0; i < b->rows; ++i) {
        size_t j;
        double result = matrix_get_at(b, i, 0);

        for (j = 0; j < i; ++j) {
            result -= matrix_get_at(A, i, j) * matrix_get_at(x, j, 0);    
        }
        result /= matrix_get_at(A, i, i);
        matrix_set_at(x, i, 0, result);
    }

    return x;
}

/** 
 *  Resolve o sistema Ax = b por substituição para trás
 *  @return O vetor x.
 */
static matrix_t *backwards_substitution(const matrix_t * restrict A, const matrix_t * restrict b)
{
    size_t i;
    matrix_t *x = matrix_new(b->rows, 1);

    /* 
     * this loop is a bit funny.
     * it relies on the fact that unsigned overflow is defined
     * and will surely be larger (or equal to) A->rows
     * this is done to avoid i >= 0, which will always return true
     * due to overflow.
     */
    for (i = A->rows - 1; i < A->rows; --i) {
        size_t j;
        double result = matrix_get_at(b, i, 0);

        for (j = i + 1; j < b->rows; ++j) {
            result -= matrix_get_at(A, i, j) * matrix_get_at(x, j, 0);
        }
        result /= matrix_get_at(A, i, i);
        matrix_set_at(x, i, 0, result);
    }

    return x;
}

#endif
