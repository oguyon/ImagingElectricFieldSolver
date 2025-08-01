#ifndef VECTOR_MATRIX_H
#define VECTOR_MATRIX_H

// Represents a vector
typedef struct {
    int size;
    double* data;
} Vector;

// Represents a matrix
typedef struct {
    int rows;
    int cols;
    double** data;
} Matrix;

Vector* create_vector(int size);
Matrix* create_matrix(int rows, int cols);
void free_vector(Vector* v);
void free_matrix(Matrix* m);

#endif // VECTOR_MATRIX_H
