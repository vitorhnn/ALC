/**
 *  @file
 */
/*
 * Copyright (c) 2016 Márcio Medeiros
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

#ifndef CONDEST_H
#define CONDEST_H

#include "matrix.h"
#include "matrix_norms.h"
#include "lu.h"

/**
 *  Calcula uma aproximação do numero condição de uma matriz A de ordem qualquer
 */
static double condest(const matrix_t *A, unsigned tests)
{
    double cond;
    matrix_t *x, *b;
    
    b = matrix_new(A->rows,1);

    double current, biggest = 0;
    unsigned k;
    for (k = 0; k < tests; k++) {
        size_t i;
        for (i = 0; i < b->rows; i++) {
            double element = rand() / 1000000000.0;
            matrix_set_at(b, i,0,element);
        }

        x = lu_solve(A, b);
            
        current = column_norm(x) / column_norm(b);

        if (current > biggest) {
            biggest = current;
        }
        
        matrix_free(x);
    }

    cond = biggest * column_norm(A);
    matrix_free(b);
    return cond;
}


#endif
