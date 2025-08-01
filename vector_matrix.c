#include "vector_matrix.h"
#include <stdlib.h>

Vector* create_vector(int size) {
    Vector* v = (Vector*)malloc(sizeof(Vector));
    v->size = size;
    v->data = (double*)malloc(size * sizeof(double));
    return v;
}

Matrix* create_matrix(int rows, int cols) {
    Matrix* m = (Matrix*)malloc(sizeof(Matrix));
    m->rows = rows;
    m->cols = cols;
    m->data = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        m->data[i] = (double*)malloc(cols * sizeof(double));
    }
    return m;
}

void free_vector(Vector* v) {
    if (v) {
        free(v->data);
        free(v);
    }
}

void free_matrix(Matrix* m) {
    if (m) {
        for (int i = 0; i < m->rows; i++) {
            free(m->data[i]);
        }
        free(m->data);
        free(m);
    }
}
