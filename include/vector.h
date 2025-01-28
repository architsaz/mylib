#ifndef VECTOR_H
#define VECTOR_H
    #include <stdlib.h> 
    #include "vector_types.h"
    int checkEIDS(int *elems);
    int compare_int_min(void *a, void *b);
    int compare_int_max(void *a, void *b);
    int compare_double_min(void *a, void *b);
    int compare_double_max(void *a, void *b);
    void *find_extreme(void *array, size_t element_size, size_t num_elements, compare_func comp);
    void crossProduct(double v1[3], double v2[3], double result[3]);
    void normalize(double vector[3]);
#endif
