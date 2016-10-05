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

#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h> /* memcpy */
#include <assert.h>
#include <float.h>

#ifndef EPSILON
#define EPSILON DBL_EPSILON
#endif

/**
 * Estrutura representando uma matriz. <br>
 * É uma boa ideia não mexer nela diretamente.
 */
typedef struct {
    double *elements;
    size_t rows;
    size_t columns;
    /** Usado internamente na hora de transpor a matriz **/
    int transposed;
} matrix_t;

/** 
 *  Seta um elemento na matriz de forma segura <br>
 *  Necessário porque ela é alocada como um bloco de tamanho linhas * colunas.
 *  @param mat A matriz a ser acessada
 *  @param row A linha do elemento a ser setado.
 *  @param column A coluna do elemento a ser setado.
 *  @param value O valor a ser setado no elemento.
 */
static void matrix_set_at(matrix_t *mat, size_t row, size_t column, double value)
{
    size_t pos;

    if (mat->transposed) {
        size_t swap = column;
        column = row;
        row = swap;
    }
    
    
    assert(row < mat->rows);
    assert(column < mat->columns);

    assert(!isnan(value));
    assert(isfinite(value));

    pos = (row * mat->columns) + column;

    mat->elements[pos] = value;
}

/**  
 *  Acessa a matriz de forma segura <br>
 *  Necessário porque ela é alocada como um bloco de tamanho linhas * colunas.
 *  @param mat A matriz a ser acessada.
 *  @param row A linha do elemento a ser acessado.
 *  @param column A coluna do element a ser acessado.
 *  @return O elemento nas posições informadas
 */
static double matrix_get_at(matrix_t *mat, size_t row, size_t column)
{
    size_t pos;

    if (mat->transposed) {
        size_t swap = column;
        column = row;
        row = swap;
    }
    
    assert(row < mat->rows);
    assert(column < mat->columns);

    pos = (row * mat->columns) + column;

    return mat->elements[pos];
}

/** 
 *  Carrega a identidade na matriz.
 *  @param mat A matriz a ser carregada.
 *  @return A matriz de entrada, carregada com a identidade.
 */
static matrix_t *matrix_load_identity(matrix_t *mat)
{
    size_t i;
    for (i = 0; i < mat->columns; ++i) {
        size_t j;
        for (j = 0; j < mat->rows; ++j) {
            if (i == j) {
                matrix_set_at(mat, i, j, 1);
            } else {
                matrix_set_at(mat, i, j, 0);
            }
        }
    }

    return mat;
}

/**  
 *   Cria uma nova matriz.
 *   @param rows Quantidade de linhas
 *   @param columns Quantidade de colunas
 *   @return A nova matriz.
 */
static matrix_t *matrix_new(size_t rows, size_t columns)
{
    matrix_t *mat = malloc(sizeof(matrix_t));
    if (!mat) {
        return NULL;
    }

    mat->elements = malloc(sizeof(double) * (rows * columns));
    if (!mat->elements) {
        free(mat);
        return NULL;
    }

    mat->rows = rows;
    mat->columns = columns;
    mat->transposed = 0;

    return mat;
}

/**
 *  Libera uma matriz.
 */
static void matrix_free(matrix_t *mat)
{
    free(mat->elements);
    free(mat);
}

/** 
 *  Gera uma cópia funda da matrix fornecida
 *  @param mat A matriz a ser copiada
 *  @return Uma cópia funda de mat
 */
static matrix_t *matrix_copy(matrix_t *mat)
{
    matrix_t *new_mat = 
            mat->transposed ? matrix_new(mat->columns, mat->rows) : 
                matrix_new(mat->rows, mat->columns);

    if (new_mat != NULL) {
        if (mat->transposed) {
            /* can't do memcpy, values would be fudged */
            size_t i;
            for (i = 0; i < new_mat->rows; ++i) {
                size_t j;
                for (j = 0; j < new_mat->columns; ++j) {
                    matrix_set_at(new_mat, i, j, matrix_get_at(mat, i, j)); 
                }
            }
        } else {
            memcpy(new_mat->elements, mat->elements, sizeof(double) * (mat->rows * mat->columns));
        }
    }

    return new_mat;
}

/**
 *  Transpõe uma matriz.
 *  @param mat A matriz a ser transposta
 *  @return mat transposta
 */
static matrix_t *matrix_transpose(matrix_t *mat)
{
    matrix_t *new_mat;

    mat->transposed = 1;
    new_mat = matrix_copy(mat);
    mat->transposed = 0;

    return new_mat;
}

/**
 *  Compara duas doubles de acordo com o epsilon
 */
static int doublecmp(double a, double b)
{
    return fabs(a - b) < EPSILON;
}

static int matrix_cmp(matrix_t *a, matrix_t *b)
{
    size_t i;

    if (a->rows != b->rows || a->columns != b->columns) {
        return 0;
    }

    for (i = 0; i < a->rows; ++i) {
        size_t j;
        for (j = 0; j < a->columns; ++j) {
            if (!doublecmp(matrix_get_at(a, i, j), matrix_get_at(b, i, j))) {
                return 0;
            }
        }
    }


    return 1;
}

/**
 *  Adiciona um escalar a a.
 *  @return a, após a soma.
 *  \warning A matriz retornada não é uma cópia de a. Caso deseje uma cópia, use matrix_copy em a antes.
 */
static matrix_t *matrix_add_scalar(matrix_t *a, double scalar)
{
    size_t i;
    for (i = 0; i < a->rows; ++i) {
        size_t j;
        for (j = 0; j < a->columns; ++j) {
            matrix_set_at(a, i, j, matrix_get_at(a, i, j) + scalar);
        }
    }

    return a;
}

/** 
 *  Soma b à a.
 *  @return a, após a soma. 
 *  \warning A matriz retornada não é uma cópia de a. Caso deseje uma cópia, use matrix_copy em a antes.
 */
static matrix_t *matrix_add_matrix(matrix_t *a, matrix_t *b)
{
    size_t i;
    
    assert(a->rows == b->rows);
    assert(a->columns == b->columns);

    for (i = 0; i < a->rows; ++i) {
        size_t j;
        for (j = 0; j < a->columns; ++j) {
            matrix_set_at(a, i, j, matrix_get_at(a, i, j) + matrix_get_at(b, i, j));
        }
    }

    return a;
}

/**
 * Subtrai um escalar de a.
 * @return a, após a subtração.
 * \warning A matriz retornada não é uma cópia de a. Caso deseje uma cópia, use matrix_copy em a antes.
 */
static matrix_t *matrix_subtract_scalar(matrix_t *a, double scalar)
{
    size_t i;
    for (i = 0; i < a->rows; ++i) {
        size_t j;
        for (j = 0; j < a->columns; ++j) {
            matrix_set_at(a, i, j, matrix_get_at(a, i, j) + scalar);
        }
    }

    return a;
}

/** 
 *  Subtrai b de a.
 *  @return a, após a subtração. 
 *  \warning A matriz retornada não é uma cópia de a. Caso deseje uma cópia, use matrix_copy em a antes.
 */
static matrix_t *matrix_subtract_matrix(matrix_t *a, matrix_t *b)
{
    size_t i;
    
    assert(a->rows == b->rows);
    assert(a->columns == b->columns);
    
    for (i = 0; i < a->rows; ++i) {
        size_t j;
        for (j = 0; j < a->columns; ++j) {
            matrix_set_at(a, i, j, matrix_get_at(a, i, j) - matrix_get_at(b, i, j));
        }
    }

    return a;
}

/** 
 *  Multiplica a por um escalar.
 *  @return a, após a multiplicação
 *  \warning A matrix retornada não é uma cópia de a. Caso deseje uma cópia, use matrix_copy em a antes.
 */
static matrix_t *matrix_mul_scalar(matrix_t *a, double scalar)
{
    size_t i;
    for (i = 0; i < a->rows; ++i) {
        size_t j;
        for (j = 0; j < a->columns; ++j) {
            matrix_set_at(a, i, j, matrix_get_at(a, i, j) * scalar);
        }
    }

    return a;
}

/**
 *  Multiplica a por b, guardando os resultados em uma nova matriz
 *  @return Uma nova matriz, com os resultados da multiplicação.
 */
static matrix_t *matrix_mul_matrix(matrix_t *a, matrix_t *b)
{
    matrix_t *result;
    size_t i;

    assert(a->columns == b->rows);

    result = matrix_new(a->rows, b->columns);
    
    for (i = 0; i < a->rows; ++i) {
        size_t j;
        for (j = 0; j < b->columns; ++j) {
            double sum = 0;
            size_t k;
            for (k = 0; k < a->columns; ++k) {
                sum += matrix_get_at(a, i, k) * matrix_get_at(b, k, j);
            }
            matrix_set_at(result, i, j, sum);
        }
    }

    return result;
}

/** 
 *  Encontra a submatriz de mat formada deletando a linha row e a coluna column
 */
static matrix_t *matrix_get_submatrix(matrix_t *mat, size_t row, size_t column)
{
    matrix_t *sub = matrix_new(mat->rows - 1, mat->columns - 1);
    size_t last_row = 0, last_column = 0;
    size_t i;

    for (i = 0; i < mat->rows; ++i) {
        size_t j;

        if (i == row) {
            continue;
        }
        
        for (j = 0; j < mat->columns; ++j) {
            if (j == column) {
                continue;
            }

            matrix_set_at(sub, last_row, last_column, matrix_get_at(mat, i, j));
            last_column++;
        }

        last_row++;
        last_column = 0;
    }

    return sub;
}

/**
 * Computa o determinante de mat via expansão de Laplace.
 * \warning Computacionalmente ineficiente, O(n!). Tome cuidado.
 */
static double matrix_get_determinant(matrix_t *mat)
{
    size_t  order = mat->rows,
            i;

    double accumulator = 0;

    assert(mat->rows == mat->columns);

    if (order == 2) {
        double primary   = matrix_get_at(mat, 0, 0) * matrix_get_at(mat, 1, 1);
        double secondary = matrix_get_at(mat, 0, 1) * matrix_get_at(mat, 1, 0);

        return primary - secondary;
    }

    for (i = 0; i < order; ++i) {
        matrix_t *sub = matrix_get_submatrix(mat, 0, i);

        accumulator += matrix_get_at(mat, 0, i) * (pow(-1, i) * matrix_get_determinant(sub));

        matrix_free(sub);
    }

    return accumulator;
}

/**
 *  Imprime uma matriz.
 */
static void matrix_pretty_print(matrix_t *mat)
{
    size_t i;
    for (i = 0; i < mat->rows; ++i) {
        size_t j;
        for (j = 0; j < mat->columns; ++j) {
            printf("%f    ", matrix_get_at(mat, i, j));
        }
        printf("\n");
    }
}

#endif /* MATRIX_H */
