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

#ifndef MATRIX_ANGLE_OPS_H
#define MATRIX_ANGLE_OPS_H

#include "matrix.h"

/**
 *  Calcula o produto interno entre a linha row e a coluna column de mat
 *  @author Andrei Parente
 */
static double row_column_dot_product(matrix_t *mat, size_t row, size_t column)
{
    double product = 0;
    size_t i;
    
    assert(mat->rows == mat->columns);

    for (i = 0; i < mat->rows; ++i) {
        product += matrix_get_at(mat, row, i)* matrix_get_at(mat, i, column);
    }
    return product; 
}

/**
 *  Calcula a magnitude do vetor vector
 */
static double vector_length(matrix_t *vector)
{
    size_t i;
    double accumulator = 0;

    assert(vector->rows == 1);
    
    for (i = 0; i < vector->columns; ++i) {
        accumulator += pow(matrix_get_at(vector, 0, i), 2);
    }

    return sqrt(accumulator);
}

/**
 *  Calcula o ângulo entre a linha row e a coluna column de mat
 *  @author Andrei Parente
 */
static double row_column_angle(matrix_t *mat, size_t row, size_t column)
{
    size_t i;
    double row_length = 0,
           column_length = 0;

    matrix_t *temp = matrix_new(1, mat->columns);
    for (i = 0; i < mat->columns; ++i) {
        matrix_set_at(temp, 0, i, matrix_get_at(mat, row, i));
    }

    row_length = vector_length(temp);
    matrix_free(temp);

    temp = matrix_new(1, mat->rows);
    for (i = 0; i < mat->rows; ++i) {
        matrix_set_at(temp, 0, i, matrix_get_at(mat, i, column)); 
    }

    column_length = vector_length(temp);
    matrix_free(temp);

    double angle = row_column_dot_product(mat, row, column) / (row_length * column_length);

    return acos(angle);
}

#endif
