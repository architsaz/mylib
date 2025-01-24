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
#endif
