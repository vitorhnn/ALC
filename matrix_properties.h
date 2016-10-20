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

#ifndef MATRIX_PROPERTIES_H
#define MATRIX_PROPERTIES_H

#include "matrix.h"
#include "cholesky.h"
#include "matrix_norms.h"
#include "vandermonde.h"

/**
 *  Confere se mat é tridiagonal.
 *  @author Andrei Parente
 */
static int tridiagonal_check(const matrix_t *mat)
{
    size_t i;
    for (i = 0; i < mat->rows; ++i) {
        size_t j;
        for (j = 0; j < mat->columns; ++j) {
            if (i == j || i == j + 1 || j == i + 1) {
                if (matrix_get_at(mat, i, j) == 0) {
                    return 0;
                }
            } else if (matrix_get_at(mat, i, j) != 0) {
                return 0;
            }
        }
    }

    return 1;
}

/**
 *  Confere se mat é ortogonal.
 *  @author Andrei Parente
 */
static int orthogonal_check(const matrix_t *mat)
{
    matrix_t *id = matrix_new(mat->rows, mat->columns);
    id = matrix_load_identity(id);

    matrix_t *transposed = matrix_transpose(mat);

    matrix_t *mulled = matrix_mul_matrix(mat, transposed);

    int result = matrix_cmp(id, mulled);

    matrix_free(transposed);
    matrix_free(mulled);

    if (result == 1) {
        return 1;
    }

    return 0;
}

/**
 *  Confere se mat é simétrica.
 *  @author Andrei Parente
 */
static int symmetric_check(const matrix_t *mat)
{
    matrix_t *transposed = matrix_transpose(mat);

    int result = matrix_cmp(mat, transposed);

    matrix_free(transposed);

    return result;
}

/**
 *  Confere se mat é positiva definida.
 *  @author Andrei Parente
 */
static int positive_definite_check(const matrix_t *mat)
{
    matrix_t *factor = cholesky_factor(mat);

    if (factor == NULL) {
        return 0;
    }

    matrix_free(factor);
    return 1;
}

/**
 *  Verifica se a matriz é estritamente diagonal dominante
 *  @autor Pedro da Luz
 */
static int strictly_dominant_diagonal_check(const matrix_t *A)
{
    matrix_t* singleRow;
    singleRow = matrix_new(1, A->columns);

    size_t i, j;
    for(i = 0; i < A->rows; i++)
    {
        for(j = 0; j < A->columns; j++)
            matrix_set_at(singleRow, 0, j, matrix_get_at(A, i, j));

        double norm = vector_norm1(singleRow);

        for(j = 0; j < A->columns; j++)
        {
            if(fabs(matrix_get_at(A, i, j)) < norm)
                return 0;
        }
    }
    matrix_free(singleRow);
    return 1;
}

/**
 *  Verifica se os vetores da matriz são lineramente independentes
 *  Pelo método da determinante em matrizes de ordem NxN
 *
 *  @param A, matriz composta por vetores linha
 *
 *  @author Pedro da Luz
 */
static int vector_lineary_independence_det_check(const matrix_t* A)
{
    if(vandermonde_check(A) == 1)
    {
        if(vandermonde_determinant(A) != 0)
            return 1;
        else
            return 0;
    }
    else
    {
        if(matrix_get_determinant(A) != 0)
            return 1;
        else
            return 0;
    }
}

#endif
