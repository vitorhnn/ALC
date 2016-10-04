/**
 *
 *  @file
 */
/*
 * Copyright (c) 2016 Daniel Juventude
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

#ifndef MATRIX_INVERSE_H
#define MATRIX_INVERSE_H

#include "matrix.h"

//funçao auxiliar para calcular a inversa
static void troca_linha(matrix_t *mat, size_t l, size_t i)
{
    double aux;
    int c;
    for (c = 0; c < mat->columns; ++c) {
        aux = matrix_get_at(mat, l, c);
        matrix_set_at(mat, l, c, matrix_get_at(mat, i, c));
        matrix_set_at(mat, i, c, aux);
    }
}

//funçao auxiliar para calcular a inversa
static void linha_div_scalar(matrix_t *mat, size_t l, double m)
{
    size_t j;
    for (j = 0; j < mat->columns; ++j) {
        matrix_set_at(mat, l, j, matrix_get_at(mat,l, j) / m);
    }
}

//funçao auxiliar para calcular a inversa
static void linha_soma_linha(matrix_t *mat, size_t l, size_t lu, double m)
{
    size_t j;
    for (j = 0; j < mat->columns; ++j) {
        matrix_set_at(mat, l, j, matrix_get_at(mat, l, j) + matrix_get_at(mat, lu, j) * m);
    }
}

/**
 *  Retorna a inversa da Matriz mat
 *  @author Daniel Juventude
 */
static matrix_t *matrix_inverse(matrix_t *mat)
{
    matrix_t *inv;
    inv = matrix_new(mat->rows, mat->columns);
    inv = matrix_load_identity(inv);
    
    double escalar;
    size_t l, c, lu;

    /*
    *   Verifica se nao ha zeros na diagonal principal,
    *e se houver, troca-se as linhas l e i; se possivel.
    */
    for (l = 0; l < mat->rows; ++l) {
        if (matrix_get_at(mat, l, l) == 0) {
            int achou = 1;
            size_t i;
            for (i = 0; (i < mat->rows) && (achou == 1); ++i) {
                if (matrix_get_at(mat, i, l) != 0) {
                    if (matrix_get_at(mat, l, i) != 0) {
                        troca_linha(mat, l, i);
                        troca_linha(inv, l, i);
                        achou = 0;
                    }
                }
            }
            assert(achou == 0);
        }
    }

    for (c = 0; c < mat->columns; ++c) {
        /*
         * transforma o elemento da diagonal principal em 1
         */
        if (matrix_get_at(mat, c, c) != 1) {
            escalar = matrix_get_at(mat, c, c);
            linha_div_scalar(mat, c, escalar);
            linha_div_scalar(inv, c, escalar);
        }
        lu = c;

        /*
         * transforma os elementos em zero
         */
        for (l = 0; l < mat->rows; ++l) {
            if (l != c) {
                if (matrix_get_at(mat, l, c) != 0) {
                    escalar = matrix_get_at(mat, l, c) * -1;
                    linha_soma_linha(mat, l, lu, escalar);
                    linha_soma_linha(inv, l, lu, escalar);
                }
            }
        }
    }
    return inv;
}

#endif
