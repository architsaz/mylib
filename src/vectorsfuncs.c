#include <stdio.h>
#include <math.h>
#include "vector_types.h"


int compare_int_min(void *a, void *b)
{
	return (*(int *)a < *(int *)b);
}
int compare_int_max(void *a, void *b)
{
	return (*(int *)a > *(int *)b);
}
int compare_double_min(void *a, void *b)
{
	return (*(double *)a < *(double *)b);
}
int compare_double_max(void *a, void *b)
{
	return (*(double *)a > *(double *)b);
}
void *find_extreme(void *array, size_t element_size, size_t num_elements, compare_func comp)
{
	void *extreme = array;

	for (size_t i = 1; i < num_elements; ++i)
	{
		void *current_element = (char *)array + i * element_size;
		if (comp(current_element, extreme))
		{
			extreme = current_element;
		}
	}

	return extreme;
}
// check the start ID of elements in the mesh file
int checkEIDS(int *elems)
{
	size_t int_size = sizeof(*elems) / sizeof(elems[0]);
	// Find min and max for int array
	int *int_min = (int *)find_extreme(elems, sizeof(int), int_size, compare_int_min);
	// printf("--> ID of elements start from %d!\n",*int_min);
	return *int_min;
}
// Function to calculate the cross product of two vectors
void crossProduct(double v1[3], double v2[3], double result[3]) {
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

// Function to normalize a vector
void normalize(double vector[3]) {
    double magnitude = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    if (magnitude > 0.0) {
        vector[0] /= magnitude;
        vector[1] /= magnitude;
        vector[2] /= magnitude;
    }
}