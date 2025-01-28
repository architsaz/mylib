#ifndef MESH_H
#define MESH_H
#include "mesh_types.h"
char *edit_endline_character(char *line, int buffer, FILE *fptr);
int read_zfem(char *path, int *npoin, int *nelem, double **ptxyz, int **elems);
int calc_area_tri3(double *ptxyz, int *elems, int nelem, double **area2);
int save_esurp(int npoin, int nelem, int *elems, int **esurp2, int **esurp_pointer2, int Nredge);
int *find_nei_elem3D(int *esurp_pointer, int *esurp, int *num_nei, int *open, int *elems, int ele, int ele_p1, int ele_p2, int Nredge);
int save_esure(int nelem, int *elems, int *esurp_pointer, int *esurp, int **esue2, int **open2, int Nredge);
int save_fsure(int nelem, int *esure, int **efid2, int *numf, int Nredge);
int save_psurf(int nelem, int numf, int *elems, int *esure, int **psurf2, int Nredge);
int save_esurf(int nelem, int *esure, int numf, int **esurf2, int Nredge);
int save_normele(int nelem, int *elems, double *ptxyz, double **norm);
int save_centri3(int nelem, int *elems, double *ptxyz, double **cen2);
void tri3_to_tri6(mesh *M1, mesh **M2);
void tri3_to_quad4(mesh *M1, mesh **M2);
int ConverMesh(mesh *M1, mesh *M2, ConvertorFunc Func);
void SCA_int_VTK(FILE *fptr, char *name, int col, int num, void *field);
void SCA_double_VTK(FILE *fptr, char *name, int col, int num, void *field);
void VEC_double_VTK(FILE *fptr, char *name, int col, int num, void *field);
void tri3funcVTK(FILE *fptr, int nelem, int *elems);
void tri6funcVTK(FILE *fptr, int nelem, int *elems);
void read_VTK_double(FILE *fptr, int col, int nr, void **field);
void read_VTK_int(FILE *fptr, int col, int nr, void **field);
int ReadVTK(char *dir, char *filenam, int step, FunctionWithArgs2 *prtfield, int countfield);
int SaveVTK(char *dir, char *filenam, int step, mesh *M, elemVTK elemfunc, FunctionWithArgs elefuncs[], size_t nrelefield, FunctionWithArgs pntfuncs[], size_t nrpntfield);
void calculateEdgeNormal(double edge[3], double faceNormal[3], double edgeNormal[3]);
int save_normedge(int nelem, double *ptxyz, int *elems, double *normele, double **normedge2);
int save_cenedgetri3(int nelem,int *elems,double *ptxyz,double **cen2);
#endif