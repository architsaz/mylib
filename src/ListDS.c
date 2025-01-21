#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "listDS_types.h"
int init_LDS (int num_lists,ListData **lists2,int **lists_ptr2,int **lists_size2,int temp_size){
    if ( temp_size == 0 || num_lists == 0){
        fprintf(stderr,"! the temporary size or number of list is not defined.\n");
        return 1;
    }
    ListData *lists = (ListData *)malloc ((size_t)num_lists*(size_t)temp_size*sizeof(ListData));
    int *lists_ptr = (int *)calloc ((size_t)(num_lists+1),sizeof(int));
    int *lists_size = (int *) calloc ((size_t)num_lists,sizeof(int));
    if (lists == NULL || lists_ptr == NULL || lists_ptr == NULL ){
        fprintf(stderr,"! there is a problem in allocation memory\n");
        return 1;
    }
    for (int i=1;i<=num_lists;i++)
        lists_ptr[i]=temp_size*i;
    
    // pass arraies
    *lists2 =lists;
    *lists_ptr2 = lists_ptr;
    *lists_size2 = lists_size;
    return 0; //success signal
}
void LDS_reloc(int listID, int num_lists, ListData **lists_ptr, int *lists_ptr_array, int temp_size) {
    // Extract the current `lists` array
    ListData *lists = *lists_ptr;

    // Calculate the new size for the `lists` array
    int newsize = lists_ptr_array[num_lists] + temp_size;

    // Safely reallocate memory for the `lists` array
    ListData *new_lists = (ListData *)realloc(lists, (size_t)newsize * sizeof(ListData));
    if (new_lists == NULL) {
        fprintf(stderr, "Reallocation failed!\n");
        exit(EXIT_FAILURE);
    }

    // Shift elements in the `lists` array to create space for the new list
    int start = lists_ptr_array[listID + 1];
    int end = lists_ptr_array[num_lists];

    // Move elements in reverse to avoid overwriting
    for (int i = end - 1; i >= start; i--) {
        new_lists[i + temp_size] = new_lists[i];
    }

    // Initialize the newly created space with 0
    for (int i = start; i < start + temp_size; i++) {
        new_lists[i].double_data = 0;
        new_lists[i].int_data = 0;
    }

    // Update `lists_ptr_array` 
    for (int i = listID + 1; i <= num_lists; i++) {
        lists_ptr_array[i] += temp_size;
    }

    // Update the pointer to the reallocated memory
    *lists_ptr = new_lists;
}
int LDS_print(int num_lists, ListData *lists, int *lists_ptr, int *lists_size) {
    if (lists == NULL || lists_ptr == NULL || lists_size == NULL) {
        fprintf(stderr, "! NULL pointer passed to the function\n");
        return 1;
    }

    for (int i = 0; i < num_lists; i++) {
        printf("- List%d size(%d): ", i, lists_size[i]);
        for (int k = 0; k < lists_size[i]; k++) {
            printf(" ");

            // Print data based on the type
            if (lists[lists_ptr[i] + k].int_data != 0) {
                printf("(Int: %d ", lists[lists_ptr[i] + k].int_data);
            }
            if (lists[lists_ptr[i] + k].double_data != 0) {
                printf("Double: %.2f )", lists[lists_ptr[i] + k].double_data);
            }
        }
        printf("\n");
    }
    return 0;
}
void LDS_insert(ListData data, int listID, int num_lists, int temp_size,ListData **lists_ptr, int *lists_ptr_array, int *lists_size) {
    // Extract the current `lists` array
    ListData *lists = *lists_ptr;

    // Check if there's enough space; if not, reallocate
    int lsize = lists_ptr_array[listID + 1] - lists_ptr_array[listID];
    if (lsize <= lists_size[listID]) {
        #ifdef DEBUG
        printf("! Relocate list memory safely\n");
        #endif
        LDS_reloc(listID, num_lists, lists_ptr, lists_ptr_array, temp_size);
        lists = *lists_ptr; // Update the local pointer after reallocation
    }

    // Insert the data
    lists[lists_ptr_array[listID] + lists_size[listID]].double_data = data.double_data;
    lists[lists_ptr_array[listID] + lists_size[listID]].int_data = data.int_data;
    lists_size[listID]++;
}


