/**
 *  @file
 */
/*
 * Copyright (c) 2016 Andrei Parente
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

#ifndef ITERATIVE_METHODS_H
#define ITERATIVE_METHODS_H

#include "matrix.h"
#include "matrix_norms.h"

static size_t __g_iterations;

/**
 *  Informa a quantidade de iterações usadas no último método iterativo chamado.
 *  @return As iterações
 */
static size_t get_iterations(void)
{
    return __g_iterations;
}

/**
 *  Confere o critério das linhas em mat
 *  @author Andrei Parente
 */
static int convergence_rows_check(matrix_t *mat)
{
    size_t i;
    for (i = 0; i < mat->rows; ++i) {
        double sum = 0;
        size_t j;
        for (j = 0; j < mat->columns; ++j) {
            if (i != j) {
                if (matrix_get_at(mat, i, j) < 0) {
                    sum -= matrix_get_at(mat, i, j);
                } else {
                    sum += matrix_get_at(mat, i, j);
                }
            } 
        }
        if (matrix_get_at(mat, i, i) < 0 && sum >= matrix_get_at(mat, i, i) * -1) {
            return 0;
        } else if (sum >= matrix_get_at(mat, i, i)) {
            return 0;
        }
    }
    return 1;
}

/**
 *  Confere o critério das colunas em mat
 *  @author Andrei Parente
 */
static int convergence_columns_check(matrix_t *mat)
{
    size_t j;
    for (j = 0; j < mat->columns; ++j) {
        double sum = 0;
        size_t i;
        for (i = 0; i < mat->rows; ++i) {
            if (i != j) {
                if (matrix_get_at(mat, i, j) < 0) {
                    sum -= matrix_get_at(mat, i, j);
                } else {
                    sum += matrix_get_at(mat, i, j);
                }
            } 
        }
        if (matrix_get_at(mat, j, j) < 0 && sum >= matrix_get_at(mat, j, j) * -1) {
            return 0;
        } else if (sum >= matrix_get_at(mat, j, j)) {
            return 0;
        }
    }

    return 1;
}

/**
 *  Confere o critério de Sassenfeld em mat
 *  @author Andrei Parente
 */
static int convergence_sassenfeld_check(matrix_t *mat)
{
    double *beta = malloc(sizeof(double) * mat->rows);

    size_t i;
    for (i = 0; i < mat->rows; ++i) {
        beta[i] = 1;
    }

    for (i = 0; i < mat->rows; ++i) {
        beta[i] = 0;
        size_t j;
        for (j = 0; j < mat->columns; ++j) {
            if (i != j) {
                if (matrix_get_at(mat, i, j) < 0) {
                    beta[i] -= matrix_get_at(mat, i, j) * beta[j];
                } else {
                    beta[i] += matrix_get_at(mat, i, j) * beta[j];
                }
            } 
        }


        if (matrix_get_at(mat, i, i) < 0) {
            beta[i] /= matrix_get_at(mat, i, i) * -1;
            if (beta[i] >= 1) {
                free(beta);
                return 0;
            }
        } else {
            beta[i] /= matrix_get_at(mat, i, i);
            if (beta[i] >= 1) {
                free(beta);
                return 0;
            }
        }
    }

    free(beta);
    return 1;
}

/**
 *  Confere o critério das normas em mat
 *  @author Andrei Parente
 */
static int convergence_norms_check(matrix_t *mat)
{
    if (row_norm(mat) < 1) {
        return 1;
    }

    if (column_norm(mat) < 1) {
        return 1;
    }

    if (frobenius_norm(mat) < 1) {
        return 1;
    }

    return 0;
}

/**
 *  Confere todos os critérios de convergência em mat
 *  @author Andrei Parente
 */
static int convergence_check_all(matrix_t *mat)
{
    if (convergence_rows_check(mat) == 1 || 
        convergence_columns_check(mat) == 1 || 
        convergence_sassenfeld_check(mat) == 1 || 
        convergence_norms_check(mat) == 1)
    {
        return 1;        
    }

    return 0;
}

/**
 *  Resolve o sistema Ax = b pelo método iterativo de Jacobi.
 *  @return O vetor x.
 *  @author Andrei Parente
 */
static matrix_t *jacobi_solve(matrix_t *A, matrix_t* b, double absolute_error)
{
    matrix_t *x1 = matrix_new(b->rows, 1);

    size_t i;
    for (i = 0; i < x1->rows; ++i) {
        matrix_set_at(x1, i, 0, 0);
    }

    size_t flag = 0;
    __g_iterations = 0;
    while (flag != x1->rows) {
        __g_iterations++;
        matrix_t *x0 = matrix_copy(x1);
        flag = 0;

        for (i = 0; i < A->rows; ++i) {
            double sum = 0;
            size_t j;
            for (j = 0; j < A->columns; ++j) {
                if (i == j) {
                    sum += matrix_get_at(b, i, 0);
                } else {
                    sum -= matrix_get_at(A, i, j)*matrix_get_at(x0, j, 0);
                }
            }
            sum /= matrix_get_at(A, i, i);
            matrix_set_at(x1, i, 0, sum);

            if (fabs(matrix_get_at(x1, i, 0) - matrix_get_at(x0, i, 0)) < absolute_error) {
                flag++;
            }
        }
        matrix_free(x0); 
    }

    return x1;
}

/**
 *  Resolve o sistema Ax = b pelo método iterativo de Sobre-Relaxação Sucessiva.
 *  @return O vetor x.
 */
static matrix_t *sor_solve(matrix_t *A, matrix_t *b, double w, double absolute_error)
{
    matrix_t *x = matrix_new(b->rows, 1);
    size_t done = 0;

    iterations = 0;

    memset(x->elements, 0, sizeof(double) * x->rows);

    while (done != x->rows) {
        iterations++;
        done = 0;
        matrix_t *prevx = matrix_copy(x);
        size_t i;
        for (i = 0; i < x->rows; ++i) {
            double xhat = matrix_get_at(b, i, 0);
            double delta;
            size_t j;

            for (j = 0; j < x->rows; ++j) {
                if (i != j) {
                    xhat -= matrix_get_at(A, i, j) * matrix_get_at(x, j, 0);
                }
            }

            xhat /= matrix_get_at(A, i, i);
            delta = xhat - matrix_get_at(x, i, 0);

            matrix_set_at(x, i, 0, matrix_get_at(x, i, 0) + w * delta);

            if (fabs(matrix_get_at(x, i, 0) - matrix_get_at(prevx, i, 0)) < absolute_error) {
                done++;
            }
        }
        

        matrix_free(prevx);
    }

    return x;
}

/**
 *  Resolve o sistema Ax = b pelo método iterativo de Gauss-Seidel.
 *  @return O vetor x.
 *  @author Márcio Medeiros
 */
static matrix_t *gauss_seidel_solve(matrix_t *A, matrix_t *b, double absolute_error)
{
    return sor_solve(A, b, 1, absolute_error);
}


#endif
