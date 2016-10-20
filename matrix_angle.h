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
#include "matrix_norms.h"

/**
 *  Calcula o produto interno entre a linha row e a coluna column de mat
 *  @author Andrei Parente
 */
static double row_column_dot_product(const matrix_t *mat, size_t row, size_t column)
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
static double vector_length(const matrix_t *vector)
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
static double row_column_angle(const matrix_t *mat, size_t row, size_t column)
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

/**
 *  Calcula a distância entre dois vetores
 *
 *  @param u e v, vetor em matriz linha
 *
 *  @author Pedro da Luz
 */
static double vector_distance_vector(const matrix_t * restrict u, const matrix_t * restrict v)
{
    double accumulator;

    size_t i;

    assert(u->rows == 1);
    assert(v->rows == 1);
    assert(u->columns == v->columns);

    for(i = 0; i < u->columns; i++)
        accumulator += pow(matrix_get_at(v, 0, i) - matrix_get_at(u, 0, i), 2);

    return sqrt(accumulator);
}

/**
 *  Calcula o ângulo interno entre dois vetores
 *
 *  @param u e v, vetores em matriz linha
 *  @return angle, valor do ângulo em graus
 *
 *  @author Pedro da Luz
 */
static double vector_angle_vector(const matrix_t * restrict u, const matrix_t * restrict v)
{
    double angle;

    double accumulator = 0;

    assert(u->rows == 1);
    assert(v->rows == 1);
    assert(u->columns == v->columns);

    size_t i;
    for(i = 0; i < u->columns; i++)
        accumulator += matrix_get_at(u, 0, i) * matrix_get_at(v, 0, i);
    
    angle = acos(accumulator / (vector_norm2(u) * vector_norm2(v))); // Arco-Cosseno em Radianos

    angle = angle * (180/M_PI); // Converte para graus

    return angle;
}

/**
 *  Calcula o produto interno entre dois vetores
 *
 *  @param u e v, vetores em matriz linha
 *  @return innner_product_space, produto interno
 *
 *  @author Pedro da Luz
 */
static double vector_innerProductSpace(const matrix_t * restrict u, const matrix_t * restrict v)
{
    double innner_product_space;

    innner_product_space = (vector_norm2(u) * vector_norm2(v)) * cos(vector_angle_vector(u, v) * (M_PI/180));
    
    return innner_product_space;    
}

#endif
