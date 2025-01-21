#ifndef LISTDS_H
#define LISTDS_H
    #include "listDS_types.h"
    int init_LDS (int num_lists,ListData **lists2,int **lists_ptr2,int **lists_size2,int temp_size);
    int LDS_print(int num_lists,ListData *lists,int *lists_ptr,int *lists_size);
    void LDS_insert(ListData data, int listID, int num_lists, int temp_size, ListData **lists_ptr, int *lists_ptr_array, int *lists_size);
#endif
