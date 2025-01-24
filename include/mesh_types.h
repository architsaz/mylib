
#ifndef MESH_TYPES_H
#define MESH_TYPES_H
#include <stdio.h>

// mesh structs
typedef struct
{
    char type[10]; // tri or quad
    int nredge;    // 3 or 4
    int nrpts;     //<3,6> or <4,8,9>
    int npoin;
    int numExtraPoints;
    double *ptxyz;
    double *extra_ptxyz;
    int nelem;
    int *elems;
    int *esurp;
    int *esurp_ptr;
    int *esure;
    int numf;
    int *fsure;
    int *psurf;
    int *esurf;
    double *normele;
    int *eledomain;
    int *open; // adjacent elements to the open region (inlet/outlet)
    // elemental and pointal mask
    int *Melem;    // wall charectristics from .wall file
    int *rpts;     // pointal value of regional mask     --> read labels_srf.zfem
    int *relems;   // elemental value of regional mask --> approximate
    double *t;     // thickness shell element
    double *young; // Young modulus for each element [dyne/cm^2]
    int *presmask; // mask of pressure applied for pre-stressed configuration
    int *fixbmask; // mask of fixed element for boundary consition
    int *BCmask;   // mask of Boundary Condition is getting its value from the BCmask.txt file.
} mesh;

// define pointer to fuc for mesh convert
typedef void ConvertorFunc(mesh *, mesh **);

// define typs for writing fields in the VTK format
typedef void (*elefieldVTK)(FILE *, char *, int, int, void *);
typedef struct
{
    char *name;
    int col;
    int nr;
    void *field;
    elefieldVTK function;
} FunctionWithArgs;

// define typs for writing fields in the VTK format
typedef void (*readfieldVTK)(FILE *, int, int, void **);
typedef void(elemVTK)(FILE *, int, int *);
typedef struct
{
    char *name;
    int col;
    int nr;
    void **arr;
    readfieldVTK function;
} FunctionWithArgs2;


#endif // MESH_TYPES_H