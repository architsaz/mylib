#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "mesh_types.h"
#include "vector.h"


char *edit_endline_character(char *line, int buffer, FILE *fptr)
{

	char *str;
	int len;

	str = fgets(line, buffer, fptr);
	len = (int)strlen(str);
	if (str[len - 1] == '\n')
		str[len - 1] = '\0';

	return str;
}
int read_zfem(char *path, int *npoin, int *nelem, double **ptxyz, int **elems)
{
    int e = 0;
    int npoin1 = 0, nelem1 = 0;
    int *elems1 = NULL;
    double *ptxyz1 = NULL;

    /* Open the file */
    FILE *fptr = fopen(path, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open file - %s.\n", path);
        return -1;
    }
    printf("File opened - %s.\n", path);

    /* Read all lines of the file */
    int buffer = 100;
    char *str;
    char line[buffer];
    int endcount = 0, nscan, iline;
    char test[20];

    while (1)
    {
        // Start reading points
        str = edit_endline_character(line, buffer, fptr);
        nscan = sscanf(str, "%s", test);

        if (!strcmp(test, "POINTS"))
        {
            printf("Reading POINTS.\n");

            /* Read Number of Points */
            str = edit_endline_character(line, buffer, fptr);
            nscan = sscanf(str, "%d", &npoin1);
            printf("Number of Points = %d.\n", npoin1);
            if (nscan != 1)
            {
                fprintf(stderr, "ERROR: Incorrect number of entries on POINTS line.\n");
                fclose(fptr);
                return -1;
            }

            /* Allocate and Read Coordinates */
            ptxyz1 = malloc(3 * (size_t)npoin1 * sizeof(*ptxyz1));
            if (!ptxyz1)
            {
                fprintf(stderr, "ERROR: Memory allocation failed for ptxyz array.\n");
                fclose(fptr);
                return -1;
            }
            for (iline = 0; iline < npoin1; iline++)
            {
                str = edit_endline_character(line, buffer, fptr);
                nscan = sscanf(str, "%lf %lf %lf",
                               &ptxyz1[3 * iline],
                               &ptxyz1[3 * iline + 1],
                               &ptxyz1[3 * iline + 2]);
                if (nscan != 3)
                {
                    fprintf(stderr, "ERROR: Incorrect coordinates on line %d of POINTS.\n", iline + 1);
                    free(ptxyz1);
                    fclose(fptr);
                    return -1;
                }
            }
            endcount += 1;
        }
        else if (!strcmp(test, "TRIANGLE"))
        {
            printf("Reading ELEMENTS.\n");

            /* Read Number of Elements */
            str = edit_endline_character(line, buffer, fptr);
            str = edit_endline_character(line, buffer, fptr);
            nscan = sscanf(str, "%d", &nelem1);
            printf("Number of ELEMENTS = %d.\n", nelem1);
            if (nscan != 1)
            {
                fprintf(stderr, "ERROR: Incorrect number of entries for ELEMENTS.\n");
                free(ptxyz1);
                fclose(fptr);
                return -1;
            }

            /* Allocate and Read Connectivity */
            elems1 = malloc(3 * (size_t)nelem1 * sizeof(*elems1));
            if (!elems1)
            {
                fprintf(stderr, "ERROR: Memory allocation failed for elems array.\n");
                free(ptxyz1);
                fclose(fptr);
                return -1;
            }

            for (iline = 0; iline < nelem1; iline++)
            {
                str = edit_endline_character(line, buffer, fptr);
                nscan = sscanf(str, "%d %d %d",
                               &elems1[3 * iline],
                               &elems1[3 * iline + 1],
                               &elems1[3 * iline + 2]);
                if (nscan != 3)
                {
                    fprintf(stderr, "ERROR: Incorrect connectivity on line %d of ELEMENTS.\n", iline + 1);
                    free(ptxyz1);
                    free(elems1);
                    fclose(fptr);
                    return -1;
                }
            }
            endcount += 1;
        }

        // Break loop if both POINTS and ELEMENTS sections are read
        if (endcount == 2)
            break;
    }

    /* Close the file */
    if (fclose(fptr) == EOF)
    {
        fprintf(stderr, "ERROR: Failed to close file %s.\n", path);
        free(ptxyz1);
        free(elems1);
        return -1;
    }

    // Assign results to output pointers
    *npoin = npoin1;
    *nelem = nelem1;
    *ptxyz = ptxyz1;
    *elems = elems1;

    printf("* Exiting function for reading flds.zfem file.\n");
    checkEIDS(elems1); // Call any post-processing/check function as needed

    return e;
}
int calc_area_tri3(double *ptxyz, int *elems, int nelem, double **area2)
{
	double p1[3], p2[3], p3[3];
	double u[3], v[3];
	int np1, np2, np3;
	double *area;

	// Allocate memory for area array
	area = (double *)calloc((size_t)nelem, sizeof(*area));
	if (area == NULL)
	{
		fprintf(stderr, "ERROR: Memory allocation failed\n");
		return 1;
	}

	// Loop through each element (triangle)
	for (int ele = 0; ele < nelem; ele++)
	{
		// Get vertex indices (assuming elems is 1-based indexing, hence the -1)
		np1 = elems[3 * ele] - 1;
		np2 = elems[3 * ele + 1] - 1;
		np3 = elems[3 * ele + 2] - 1;

		// Retrieve the coordinates of the vertices
		for (int i = 0; i < 3; i++)
			p1[i] = ptxyz[3 * np1 + i];
		for (int i = 0; i < 3; i++)
			p2[i] = ptxyz[3 * np2 + i];
		for (int i = 0; i < 3; i++)
			p3[i] = ptxyz[3 * np3 + i];

		// Compute vectors u = p2 - p1 and v = p3 - p1
		for (int i = 0; i < 3; i++)
			u[i] = p2[i] - p1[i];
		for (int i = 0; i < 3; i++)
			v[i] = p3[i] - p1[i];

		// Compute the cross product of u and v and find the magnitude
		area[ele] = 0.5 * sqrt(
							  pow(u[1] * v[2] - u[2] * v[1], 2) +
							  pow(u[2] * v[0] - u[0] * v[2], 2) +
							  pow(u[0] * v[1] - u[1] * v[0], 2));

		// Check for invalid (non-positive) area values
		if (area[ele] <= 0)
		{
			fprintf(stderr, "ERROR: the area of ele[%d] : %lf\n", ele, area[ele]);
			free(area); // Free allocated memory in case of error
			return 1;
		}
	}

	// Assign the calculated areas to the output pointer
	*area2 = area;
	printf("* area of each tri3 calculated!\n");
	return 0;
}
// make data structure for elements surrounding a point
int save_esurp(int npoin, int nelem, int *elems, int **esurp2, int **esurp_pointer2, int Nredge)
{
	int e = 0;
	// check the start ID element
	if (checkEIDS(elems) != 1)
	{
		fprintf(stderr, "ERROR: The element ID should start from 1 for save_esurp function!\n");
		return -1;
	}
	// define parameter
	int *pointer, *pointer2, *esurp;
	// allocate memory for pointer
	pointer = calloc((size_t)npoin + 2, sizeof(*(pointer)));
	pointer2 = calloc((size_t)npoin + 2, sizeof(*(pointer2)));
	// find the nr of Elements surround each point
	for (int ele = 0; ele < nelem; ele++)
	{
		for (int i = 0; i < Nredge; i++)
			pointer[elems[Nredge * ele + i] + 1]++;
	}

	for (int pt = 1; pt <= npoin; pt++)
		pointer[pt + 1] += pointer[pt];
	for (int pt = 0; pt <= npoin + 1; pt++)
		pointer2[pt] = pointer[pt];
	// allocate memory for the esurp
	esurp = malloc((size_t)pointer[npoin + 1] * sizeof(*(esurp)));

	// find elements surround each point
	for (int ele = 0; ele < nelem; ele++)
	{
		for (int i = 0; i < Nredge; i++)
			esurp[pointer[elems[Nredge * ele + i]]++] = ele;
	}

	// done
	free(pointer);
	*esurp2 = esurp;
	*esurp_pointer2 = pointer2;
	return e;
}
// find neighbor element arround each element from esurp data structure
int *find_nei_elem3D(int *esurp_pointer, int *esurp, int *num_nei, int *open, int *elems, int ele, int ele_p1, int ele_p2, int Nredge)
{

	int elemnum, *p, *order;
	static int *nei; // if nei[0] is -9999 it means that there is problem to find a neighbour element ----------->nei[1] indicate to the number of neighbour in the nei[0]
	nei = calloc((size_t)2, sizeof(*(nei)));
	nei[0] = -9999;
	p = calloc((size_t)Nredge, sizeof(*p));
	order = calloc(2 * (size_t)Nredge, sizeof(*order));
	for (int i = 1; i < Nredge; i++)
	{
		order[2 * i - 1] = i;
		order[2 * i] = i;
	}
	/* find neighbour of ele */
	int *lesps; // list of element around ele_p1 and ele_p2
	int j = 0;
	int nr = esurp_pointer[ele_p1 + 1] + esurp_pointer[ele_p2 + 1] - esurp_pointer[ele_p1] - esurp_pointer[ele_p2];
	lesps = calloc((size_t)nr, sizeof(*lesps));
	for (int i = esurp_pointer[ele_p1]; i < esurp_pointer[ele_p1 + 1]; i++)
	{
		lesps[j] = esurp[i];
		j++;
	}
	for (int i = esurp_pointer[ele_p2]; i < esurp_pointer[ele_p2 + 1]; i++)
	{
		lesps[j] = esurp[i];
		j++;
	}

	for (int i = 0; i < nr; i++)
	{

		elemnum = lesps[i];
		if (num_nei[elemnum] < Nredge && open[elemnum] == 0)
		{

			if (elemnum == ele)
				continue; // checking the same element ID

			for (int k = 0; k < Nredge; k++)
				p[k] = elems[Nredge * elemnum + k];
			for (int k = 0; k < Nredge; k++)
			{
				if (p[order[2 * k]] == ele_p1 && p[order[2 * k + 1]] == ele_p2)
				{
					nei[0] = elemnum;
					nei[1] = k;
					break;
				}
				if (p[order[2 * k]] == ele_p2 && p[order[2 * k + 1]] == ele_p1)
				{
					nei[0] = elemnum;
					nei[1] = k;
					break;
				}
			}
		}
	}
	free(lesps);
	free(order);
	free(p);
	return nei;
}
// make data structure for elements surrounding an element
int save_esure(int nelem, int *elems, int *esurp_pointer, int *esurp, int **esue2, int **open2, int Nredge)
{
	int e = 0;
	// check the start ID element
	if (checkEIDS(elems) != 1)
	{
		fprintf(stderr, "ERROR: The element ID should start from 1 for save_esurp function!\n");
		return -1;
	}
	// parameters
	int *p, *order, *nei, *num_nei, *open;

	/* Allocate space to nei pointer */
	p = calloc((size_t)Nredge, sizeof(*p));
	order = calloc(2 * (size_t)Nredge, sizeof(*order));
	nei = calloc((size_t)Nredge * (size_t)nelem, sizeof(*(nei)));
	num_nei = calloc((size_t)nelem, sizeof(*(num_nei)));
	open = calloc((size_t)nelem, sizeof(*(open)));
	// initializing
	for (int ele = 0; ele < nelem; ele++)
	{
		for (int j = 0; j < Nredge; j++)
			nei[Nredge * ele + j] = -1;
	}
	for (int i = 1; i < Nredge; i++)
	{
		order[2 * i - 1] = i;
		order[2 * i] = i;
	}

	for (int ele = 0; ele < nelem; ele++)
	{
		for (int j = 0; j < Nredge; j++)
			p[j] = elems[Nredge * ele + j];

		// controller condition
		if (num_nei[ele] == Nredge)
			continue;
		for (int j = 0; j < Nredge; j++)
		{
			if (nei[Nredge * ele + j] == -1)
			{
				int *out = find_nei_elem3D(esurp_pointer, esurp, num_nei, open, elems, ele, p[order[2 * j]], p[order[2 * j + 1]], Nredge);
				if (out[0] != -9999)
				{
					nei[Nredge * ele + j] = out[0];
					num_nei[ele]++;
					num_nei[out[0]]++;
					nei[Nredge * out[0] + out[1]] = ele;
				}
				else
				{
					nei[Nredge * ele + j] = -2;
				}
				free(out); // free `out` returned by `find_nei_elem3D` after each use
			}
		}
		// find the element adjacent to hole
		if (num_nei[ele] < Nredge)
		{
			open[ele] = 1;
			// printf("the element %d is near holes\n",ele);
		}
	}
	// Done;
	free(num_nei);
	free(p);
	free(order);
	*esue2 = nei;
	*open2 = open;
	return e;
}
// find Nr of eadge in the mesh and make data structure for adges surrounding an element
int save_fsure(int nelem, int *esure, int **efid2, int *numf, int Nredge)
{
	int e = 0;
	static int *efid;
	int nei, ele, f;
	int num = 0;

	// allocate memory
	efid = calloc((size_t)Nredge * (size_t)nelem, sizeof(*efid));
	for (int i = 0; i < (Nredge * nelem); i++)
		efid[i] = -1;

	for (ele = 0; ele < nelem; ele++)
	{
		for (f = 0; f < Nredge; f++)
		{
			if (efid[Nredge * ele + f] < 0)
			{
				nei = esure[Nredge * ele + f];
				if (nei >= 0)
				{ // this is not boundary face
					efid[Nredge * ele + f] = num;
					for (int j = 0; j < Nredge; j++)
					{
						if (esure[Nredge * nei + j] == ele)
							efid[Nredge * nei + j] = num;
					}
				}
				else
				{ // this is on the boundary face
					efid[Nredge * ele + f] = nei;
				}
				num++;
			}
		}
	}
	printf("* nr face : %d\n", num);
	*numf = num;
	*efid2 = efid;
	return e;
}
// make data structure for points surrounding an edge
int save_psurf(int nelem, int numf, int *elems, int *esure, int **psurf2, int Nredge)
{
	int e = 0;
	int *psurf, *order;
	int *efid, nei, *p;
	int num = 0;

	// allocate memory
	p = calloc((size_t)Nredge, sizeof(*p));
	order = calloc(2 * (size_t)Nredge, sizeof(*order));
	efid = calloc((size_t)Nredge * (size_t)nelem, sizeof(*efid));
	for (int i = 0; i < (Nredge * nelem); i++)
		efid[i] = -1;
	psurf = calloc(2 * (size_t)numf, sizeof(*psurf));

	for (int i = 1; i < Nredge; i++)
	{
		order[2 * i - 1] = i;
		order[2 * i] = i;
	}

	for (int ele = 0; ele < nelem; ele++)
	{
		for (int f = 0; f < Nredge; f++)
		{
			if (efid[Nredge * ele + f] < 0)
			{
				for (int j = 0; j < Nredge; j++)
					p[j] = elems[Nredge * ele + j];
				for (int j = 0; j < Nredge; j++)
				{
					if (f == j)
					{
						psurf[2 * num] = p[order[2 * j]];
						psurf[2 * num + 1] = p[order[2 * j + 1]];
					}
				}
				nei = esure[Nredge * ele + f];
				if (nei >= 0)
				{ // this is not boundary face
					efid[Nredge * ele + f] = num;
					for (int j = 0; j < Nredge; j++)
					{
						if (esure[Nredge * nei + j] == ele)
							efid[Nredge * nei + j] = num;
					}
				}
				num++;
			}
		}
	}

	*psurf2 = psurf;
	free(efid);
	free(p);
	free(order);
	printf("* psurf is done.\n");
	return e;
}
// make data structure for elements surrounding an edge
int save_esurf(int nelem, int *esure, int numf, int **esurf2, int Nredge)
{

	int e = 0;
	int *esurf;
	int *efid, nei, ele, f;
	int num = 0;

	// allocate memory
	efid = calloc((size_t)Nredge * (size_t)nelem, sizeof(*efid));
	for (int i = 0; i < (Nredge * nelem); i++)
		efid[i] = -1;
	esurf = calloc(2 * (size_t)numf, sizeof(*esurf));

	for (ele = 0; ele < nelem; ele++)
	{
		for (f = 0; f < Nredge; f++)
		{
			if (efid[Nredge * ele + f] < 0)
			{
				nei = esure[Nredge * ele + f];
				esurf[2 * num] = ele;
				esurf[2 * num + 1] = nei;
				if (nei >= 0)
				{ // this is not boundary face
					efid[Nredge * ele + f] = num;
					if (esure[Nredge * nei] == ele)
						efid[Nredge * nei] = num;
					if (esure[Nredge * nei + 1] == ele)
						efid[Nredge * nei + 1] = num;
					if (esure[Nredge * nei + 2] == ele)
						efid[Nredge * nei + 2] = num;
				}
				num++;
			}
		}
	}
	free(efid);
	*esurf2 = esurf;
	printf("* esurf is done!\n");
	return e;
}
// find the normal of each element in 3D
int save_normele(int nelem, int *elems, double *ptxyz, double **norm)
{
	int e = 0;
	double *norm2 = NULL;

	// Allocate memory for norm2
	norm2 = calloc(3 * (size_t)nelem, sizeof(*norm2));
	if (norm2 == NULL)
	{
		// Memory allocation failed
		fprintf(stderr, "Memory allocation failed!\n");
		return -1; // Return an error code
	}

	double p1[3] = {0, 0, 0};
	double p2[3] = {0, 0, 0};
	double p3[3] = {0, 0, 0};
	double u[3] = {0, 0, 0};
	double v[3] = {0, 0, 0};

	for (int ele = 0; ele < nelem; ele++)
	{
		for (int i = 0; i < 3; i++)
			p1[i] = ptxyz[3 * (elems[3 * ele] - 1) + i];
		for (int i = 0; i < 3; i++)
			p2[i] = ptxyz[3 * (elems[3 * ele + 1] - 1) + i];
		for (int i = 0; i < 3; i++)
			p3[i] = ptxyz[3 * (elems[3 * ele + 2] - 1) + i];

		for (int i = 0; i < 3; i++)
			u[i] = p2[i] - p1[i];
		for (int i = 0; i < 3; i++)
			v[i] = p3[i] - p2[i];

		norm2[3 * ele] = u[1] * v[2] - u[2] * v[1];
		norm2[3 * ele + 1] = u[2] * v[0] - u[0] * v[2];
		norm2[3 * ele + 2] = u[0] * v[1] - u[1] * v[0];
		double mag = 0;
		for (int i = 0; i < 3; i++)
			mag += norm2[3 * ele + i] * norm2[3 * ele + i];
		if (mag == 0)
		{
			// zero magnitude normal vetor
			fprintf(stderr, "zero magnitude normal vetor at ele : %d!\n", ele);
			return -1; // Return an error code
		}
		for (int i = 0; i < 3; i++)
			norm2[3 * ele + i] = norm2[3 * ele + i] / sqrt(mag);
	}

	*norm = norm2; // Assign the calculated normals to the output pointer
	printf("* normele is done!\n");

	return e;
}
// calculate the center of surface element
int save_centri3(int nelem, int *elems, double *ptxyz, double **cen2)
{
	// Determine if elems starts from 0 or 1
	int min_elems = elems[0];
	for (int i = 0; i < 3 * nelem; i++)
	{
		if (min_elems > elems[i])
			min_elems = elems[i];
	}
	if (min_elems != 1)
	{
		fprintf(stderr, "Problem in elems array: it starts from %d\n", min_elems);
		return 1;
	}
	// Allocate memory for center points
	double *cen = calloc(3 * (size_t)nelem, sizeof(double));
	if (cen == NULL)
	{
		fprintf(stderr, "Memory allocation failed.\n");
		return 1;
	}

	// Calculate the center points
	for (int ele = 0; ele < nelem; ele++)
	{
		int p1 = elems[3 * ele] - 1;
		int p2 = elems[3 * ele + 1] - 1;
		int p3 = elems[3 * ele + 2] - 1;

		for (int i = 0; i < 3; i++)
		{
			cen[3 * ele + i] = (ptxyz[3 * p1 + i] + ptxyz[3 * p2 + i] + ptxyz[3 * p3 + i]) / 3.0;
		}
	}

	*cen2 = cen;
	return 0;
}
// conver mesh from tri3 to other type of mesh
void tri3_to_tri6(mesh *M1, mesh **M2)
{

	// define type of mesh for M2
	strcpy((*M2)->type, "tri");
	(*M2)->nredge = 3;
	(*M2)->nrpts = 6;
	//  define coordinate and elems for M2
	static double *ptxyz2;
	static int npoin2, *elems2, nelem2;
	// allocate memmory
	nelem2 = M1->nelem;
	npoin2 = M1->numf + M1->npoin;
	ptxyz2 = calloc(3 * (size_t)npoin2, sizeof(*ptxyz2));
	elems2 = calloc(6 * (size_t)M1->nelem, sizeof(*elems2));

	// coordinate of all(new+old) points
	for (int i = 0; i < (3 * M1->npoin); i++)
		ptxyz2[i] = M1->ptxyz[i];
	for (int i = 0; i < M1->numf; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ptxyz2[3 * M1->npoin + 3 * i + j] = (M1->ptxyz[3 * (M1->psurf[2 * i] - 1) + j] + M1->ptxyz[3 * (M1->psurf[2 * i + 1] - 1) + j]) / 2;
		}
	}
	// 	new connectivity
	for (int i = 0; i < M1->nelem; i++)
	{
		elems2[6 * i + 0] = M1->elems[3 * i + 0];
		elems2[6 * i + 1] = M1->elems[3 * i + 1];
		elems2[6 * i + 2] = M1->elems[3 * i + 2];
		elems2[6 * i + 3] = M1->fsure[3 * i + 0] + M1->npoin + 1;
		elems2[6 * i + 4] = M1->fsure[3 * i + 1] + M1->npoin + 1;
		elems2[6 * i + 5] = M1->fsure[3 * i + 2] + M1->npoin + 1;
	}
	// return:
	(*M2)->npoin = npoin2;
	(*M2)->elems = elems2;
	(*M2)->ptxyz = ptxyz2;
	(*M2)->nelem = nelem2;
	// all other data structure same as M1
	(*M2)->Melem = (int *)malloc((size_t)M1->nelem * sizeof(int));
	memcpy((*M2)->Melem, M1->Melem, (size_t)M1->nelem * sizeof(int));
	(*M2)->rpts = (int *)malloc((size_t)M1->npoin * sizeof(int));
	memcpy((*M2)->rpts, M1->rpts, (size_t)M1->npoin * sizeof(int));
	(*M2)->relems = (int *)malloc((size_t)M1->nelem * sizeof(int));
	memcpy((*M2)->relems, M1->relems, (size_t)M1->nelem * sizeof(int));
	//(*M2)->Melem = M1->Melem;	// wall charectristics from .wall file
	//(*M2)->rpts = M1->rpts;		// pointal value of regional mask     --> read labels_srf.zfem
	//(*M2)->relems = M1->relems; // elemental value of regional mask --> approximate
	printf("* the tri3 mesh converted to the tri6 mesh.\n- new npoin: %d\n- new nelem: %d\n", npoin2, nelem2);
}
void tri3_to_quad4(mesh *M1, mesh **M2)
{
	// define type of mesh for M2
	strcpy((*M2)->type, "quad");
	(*M2)->nredge = 4;
	(*M2)->nrpts = 4;
	//  define coordinate and elems for M2
	double *ptxyz2;
	int npoin2, *elems2, nelem2, *Melem2, *relems2;

	// allocate memmory
	npoin2 = M1->npoin + M1->nelem + M1->numf;
	nelem2 = 3 * M1->nelem;
	ptxyz2 = calloc(3 * (size_t)npoin2, sizeof(*ptxyz2));
	elems2 = calloc(4 * (size_t)nelem2, sizeof(*elems2));
	Melem2 = calloc((size_t)nelem2, sizeof(*Melem2));
	relems2 = calloc((size_t)nelem2, sizeof(*relems2));

	// coordinate of all(new+old) points
	for (int i = 0; i < (3 * M1->npoin); i++)
		ptxyz2[i] = M1->ptxyz[i];
	for (int i = 0; i < M1->nelem; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ptxyz2[3 * M1->npoin + 3 * i + j] = (M1->ptxyz[3 * (M1->elems[3 * i] - 1) + j] + M1->ptxyz[3 * (M1->elems[3 * i + 1] - 1) + j] + M1->ptxyz[3 * (M1->elems[3 * i + 2] - 1) + j]) / 3;
		}
	}
	for (int i = 0; i < M1->numf; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ptxyz2[3 * M1->npoin + 3 * M1->nelem + 3 * i + j] = (M1->ptxyz[3 * (M1->psurf[2 * i] - 1) + j] + M1->ptxyz[3 * (M1->psurf[2 * i + 1] - 1) + j]) / 2;
		}
	}
	// 	new connectivity
	for (int i = 0; i < M1->nelem; i++)
	{
		// first quadrilateral
		elems2[12 * i + 0 * 4 + 0] = M1->elems[3 * i + 0];
		elems2[12 * i + 0 * 4 + 1] = M1->fsure[3 * i + 0] + M1->nelem + M1->npoin + 1;
		elems2[12 * i + 0 * 4 + 2] = M1->npoin + i + 1;
		elems2[12 * i + 0 * 4 + 3] = M1->fsure[3 * i + 2] + M1->nelem + M1->npoin + 1;
		Melem2[3 * i] = M1->Melem[i];
		relems2[3 * i] = M1->relems[i];

		// second quadrilateral
		elems2[12 * i + 1 * 4 + 0] = M1->elems[3 * i + 1];
		elems2[12 * i + 1 * 4 + 1] = M1->fsure[3 * i + 1] + M1->nelem + M1->npoin + 1;
		elems2[12 * i + 1 * 4 + 2] = M1->npoin + i + 1;
		elems2[12 * i + 1 * 4 + 3] = M1->fsure[3 * i + 0] + M1->nelem + M1->npoin + 1;
		Melem2[3 * i + 1] = M1->Melem[i];
		relems2[3 * i + 1] = M1->relems[i];

		// third quadrilateral
		elems2[12 * i + 2 * 4 + 0] = M1->elems[3 * i + 2];
		elems2[12 * i + 2 * 4 + 1] = M1->fsure[3 * i + 2] + M1->nelem + M1->npoin + 1;
		elems2[12 * i + 2 * 4 + 2] = M1->npoin + i + 1;
		elems2[12 * i + 2 * 4 + 3] = M1->fsure[3 * i + 1] + M1->nelem + M1->npoin + 1;
		Melem2[3 * i + 2] = M1->Melem[i];
		relems2[3 * i + 2] = M1->relems[i];
	}
	// return:
	(*M2)->npoin = npoin2;
	(*M2)->elems = elems2;
	(*M2)->ptxyz = ptxyz2;
	(*M2)->nelem = nelem2;
	// other data structure :
	(*M2)->Melem = Melem2;
	(*M2)->relems = relems2;

	printf("the tri3 mesh converted to the quad4 mesh.\n- new npoin: %d\n- new nelem: %d\n", npoin2, nelem2);
}
int ConverMesh(mesh *M1, mesh *M2, ConvertorFunc Func)
{
	int e = 0;
	// find element surround a point
	CHECK_ERROR(save_esurp(M1->npoin, M1->nelem, M1->elems, &M1->esurp, &M1->esurp_ptr, M1->nredge));
	// find element surround an element
	CHECK_ERROR(save_esure(M1->nelem, M1->elems, M1->esurp_ptr, M1->esurp, &M1->esure, &M1->open, M1->nredge));
	// find Nr of eadge and given id to adges*/
	CHECK_ERROR(save_fsure(M1->nelem, M1->esure, &M1->fsure, &M1->numf, M1->nredge));
	printf(" the number of face : %d \n", M1->numf);
	// find point surround a face*/
	CHECK_ERROR(save_psurf(M1->nelem, M1->numf, M1->elems, M1->esure, &M1->psurf, M1->nredge));
	// convert M1 mesh to M2
	Func(M1, &M2);
	return e;
}
void SCA_int_VTK(FILE *fptr, char *name, int col, int num, void *field)
{
	int *int_field = (int *)field;
	fprintf(fptr, "SCALARS %s int %d\nLOOKUP_TABLE default\n\n", name, col);
	for (int ie = 0; ie < num; ie++)
	{
		fprintf(fptr, "%d\n", int_field[ie]);
	}
}
void SCA_double_VTK(FILE *fptr, char *name, int col, int num, void *field)
{
	double *double_field = (double *)field;
	fprintf(fptr, "SCALARS %s double %d\nLOOKUP_TABLE default\n\n", name, col);
	for (int ie = 0; ie < num; ie++)
	{
		fprintf(fptr, "%lf\n", double_field[ie]);
	}
}
void VEC_double_VTK(FILE *fptr, char *name, int col, int num, void *field)
{
	double *double_field = (double *)field;
	fprintf(fptr, "VECTORS %s double\n", name);
	for (int ie = 0; ie < num; ie++)
	{
		for (int j = 0; j < col; j++)
			fprintf(fptr, "%lf ", double_field[col * ie + j]);
		fprintf(fptr, "\n");
	}
}
void tri3funcVTK(FILE *fptr, int nelem, int *elems)
{
	fprintf(fptr, "CELLS %d %d\n", nelem, 4 * nelem);
	for (int ie = 0; ie < nelem; ie++)
	{
		fprintf(fptr, "3 %d %d %d\n", elems[3 * ie] - 1, elems[3 * ie + 1] - 1, elems[3 * ie + 2] - 1);
	}
	fprintf(fptr, "\n");

	fprintf(fptr, "CELL_TYPES %d\n", nelem);
	for (int ie = 0; ie < nelem; ie++)
	{
		fprintf(fptr, "5\n");
	}
	fprintf(fptr, "\n");
}
void tri6funcVTK(FILE *fptr, int nelem, int *elems)
{
	fprintf(fptr, "CELLS %d %d\n", nelem, 7 * nelem);
	for (int ie = 0; ie < nelem; ie++)
	{
		fprintf(fptr, "6 %d %d %d %d %d %d\n", elems[6 * ie] - 1, elems[6 * ie + 1] - 1, elems[6 * ie + 2] - 1, elems[6 * ie + 3] - 1, elems[6 * ie + 4] - 1, elems[6 * ie + 5] - 1);
	}
	fprintf(fptr, "\n");

	fprintf(fptr, "CELL_TYPES %d\n", nelem);
	for (int ie = 0; ie < nelem; ie++)
	{
		fprintf(fptr, "22\n");
	}
	fprintf(fptr, "\n");
}
void read_VTK_double(FILE *fptr, int col, int nr, void **field)
{
	void *arr;
	int buffer = 100;
	char line[buffer];
	char *token;
	const char delimiters[] = " \t\n"; // Delimiters: space, tab, and newline
	int nscan;
	arr = malloc((size_t)col * (size_t)nr * sizeof(double));
	for (int iline = 0; iline < nr; iline++)
	{
		if (fgets(line, buffer, fptr) == NULL)
		{
			if (feof(fptr))
			{
				break;
			}
			else
			{
				exit(EXIT_FAILURE);
			}
		}
		nscan = 0;
		// Get the first token
		token = strtok(line, delimiters);

		// Continue getting tokens until NULL is returned
		while (token != NULL)
		{
			// printf("Token: %s\n", token);
			((double *)arr)[iline] = atof(token);
			// sscanf(token, "%lf", &arr[iline]);
			token = strtok(NULL, delimiters);
			nscan++;
		}
		if (nscan != col)
		{
			fprintf(stderr, "ERROR: Incorrect number of coordinates on line %d of POINTS.\n", iline + 1);
			exit(EXIT_FAILURE);
		}
	}
	*field = arr;
}
void read_VTK_int(FILE *fptr, int col, int nr, void **field)
{
	void *arr;
	int buffer = 100;
	char line[buffer];
	char *token;
	const char delimiters[] = " \t\n"; // Delimiters: space, tab, and newline
	int nscan;
	arr = malloc((size_t)col * (size_t)nr * sizeof(int));
	for (int iline = 0; iline < nr; iline++)
	{
		if (fgets(line, buffer, fptr) == NULL)
		{
			if (feof(fptr))
			{
				break;
			}
			else
			{
				exit(EXIT_FAILURE);
			}
		}
		nscan = 0;
		// Get the first token
		token = strtok(line, delimiters);

		// Continue getting tokens until NULL is returned
		while (token != NULL)
		{
			// printf("Token: %s\n", token);
			// sscanf(token, "%d", &arr[iline]);
			((int *)arr)[iline] = atoi(token);
			token = strtok(NULL, delimiters);
			nscan++;
		}
		if (nscan != col)
		{
			fprintf(stderr, "ERROR: Incorrect number of coordinates on line %d of POINTS.\n", iline + 1);
			exit(EXIT_FAILURE);
		}
	}
	*field = arr;
}
int ReadVTK(char *dir, char *filenam, int step, FunctionWithArgs2 *prtfield, int countfield)
{
	int e = 0;
	char num[10];
	sprintf(num, "%d", step);
	char path[500];
	strcpy(path, dir);
	strcat(path, filenam);
	strcat(path, "_");
	strcat(path, num);
	strcat(path, ".vtk");
	/* define File pointer:*/
	FILE *fptr;
	fptr = calloc(1, sizeof(*fptr));
	printf("open file - %s.\n", path);
	/* Opening File */
	fptr = fopen(path, "r");
	if (fptr == NULL)
	{
		fprintf(stderr, "ERROR: Cannot open file - %s.\n", path);
		return -1;
	}
	/* Read all lines of the file */
	int buffer = 100;
	char *str;
	char line[buffer];
	int endcount = 0;

	char test1[20], test[20];

	while (1)
	{
		// start reading points:
		str = edit_endline_character(line, buffer, fptr);
		sscanf(str, "%s %s ", test1, test);
		for (int ifield = 0; ifield < countfield; ifield++)
		{
			if (!strcmp(test, prtfield[ifield].name))
			{
				printf("    Reading %s.\n", prtfield[ifield].name);
				/* Read header of field */
				str = edit_endline_character(line, buffer, fptr);
				sscanf(str, "%s", test1);
				if (!strcmp(test1, "LOOKUP_TABLE"))
				{
					str = edit_endline_character(line, buffer, fptr);
					sscanf(str, "%s", test1);
				}
				if (!strcmp(test1, ""))
				{
					str = edit_endline_character(line, buffer, fptr);
				}

				/* Read value of field */
				prtfield[ifield].function(fptr, prtfield[ifield].col, prtfield[ifield].nr, prtfield[ifield].arr);
				endcount += 1;
			}
		}
		if (endcount == countfield)
		{
			printf("  Done Reading all %d fields.\n", countfield);
			break;
		}
	}
	if (fclose(fptr) == EOF)
	{
		// If fclose returns EOF, it means there was an error closing the file
		printf("Error closing %s\n", path);
		return -1;
	}

	return e;
}
int SaveVTK(char *dir, char *filenam, int step, mesh *M, elemVTK elemfunc, FunctionWithArgs elefuncs[], size_t nrelefield, FunctionWithArgs pntfuncs[], size_t nrpntfield)
{
	int e = 0;
	// check the start ID element
	if (checkEIDS(M->elems) != 1)
	{
		fprintf(stderr, "ERROR: The element ID should start from 1 for SaveVTK function!\n");
		return -1;
	}
	char num[10];
	sprintf(num, "%d", step);
	char path[500];
	strcpy(path, dir);
	strcat(path, filenam);
	strcat(path, "_");
	strcat(path, num);
	strcat(path, ".vtk");
	char command[500];
	strcpy(command, "rm ");
	strcat(command, path);
	/* define File pointer:*/
	FILE *fptr;
	fptr = calloc(1, sizeof(*fptr));
	/* Opening File */
	fptr = fopen(path, "w");
	if (fptr == NULL)
	{
		fprintf(stderr, "ERROR: Cannot open file - %s.\n", path);
		return -1;
	}
	/*write the header of file : */
	fprintf(fptr, "# vtk DataFile Version 3.0\n");
	fprintf(fptr, "3D Unstructured Surface Grid  with %s%d mesh type\n", M->type, M->nrpts);
	fprintf(fptr, "ASCII\n\n");
	/*write the position of file : */
	fprintf(fptr, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fptr, "POINTS %d float\n", M->npoin + M->numExtraPoints);
	for (int ip = 0; ip < M->npoin; ip++)
	{
		fprintf(fptr, "%lf %lf %lf\n", M->ptxyz[3 * ip], M->ptxyz[3 * ip + 1], M->ptxyz[3 * ip + 2]);
	}
	for (int ip = 0; ip < M->numExtraPoints; ip++)
	{
		fprintf(fptr, "%lf %lf %lf\n", M->extra_ptxyz[3 * ip], M->extra_ptxyz[3 * ip + 1], M->extra_ptxyz[3 * ip + 2]);
	}
	fprintf(fptr, "\n");
	/*write the elems and cell type : */
	elemfunc(fptr, M->nelem, M->elems);

	// write SCALER pointal fields in the file:
	if (nrpntfield != 0)
		fprintf(fptr, "POINT_DATA %d\n", M->npoin + M->numExtraPoints);
	for (size_t i = 0; i < nrpntfield; ++i)
	{
		pntfuncs[i].function(fptr, pntfuncs[i].name, pntfuncs[i].col, pntfuncs[i].nr, pntfuncs[i].field); // Call each function with its array and size
	}
	// write SCALER elemental fields in the file:
	if (nrelefield != 0)
		fprintf(fptr, "CELL_DATA %d\n", M->nelem);
	for (size_t i = 0; i < nrelefield; ++i)
	{
		elefuncs[i].function(fptr, elefuncs[i].name, elefuncs[i].col, elefuncs[i].nr, elefuncs[i].field); // Call each function with its array and size
																										  // printf("field %ld done.\n",i);
	}
	if (fclose(fptr) == EOF)
	{
		// If fclose returns EOF, it means there was an error closing the file
		printf("Error closing %s\n", path);
		return -1;
	}
	printf("* wrote %s in the VTK format!\n", path);
	return e;
}