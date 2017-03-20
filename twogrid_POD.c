/*************************************************************************
	> File Name: APODheat.c
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2016年10月21日 星期五 21时26分17秒
 ************************************************************************/

#include "phg.h"
#include "fun.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <math.h>
#include <stdlib.h>

int
main(int argc, char *argv[])
{
    int i, j, nT, k, kk, rank, numprocs, TotalN, TotalN1, numt1, tmpt;
	GRID *g;
	MAT *A, *G1, *G2, *CosG1, *CosG2, *C;
	VEC *b, *V, *V1, *f_h;
	MAP *map, *map1;
	SOLVER *solver;
	ELEMENT *e;
    DOF *u_h;  /* numerical solution at time t_n */
    DOF  *B1, *B2, *B3, *C1, *C2, *C3, *grad_u, *u_p, **uhs, **FEMuhs, **coarseuhs;
    FLOAT tt[3], tt1[3], tt2[3], *svdM;
	FLOAT *U, *S, value, *U1, *S1;
    int  context, flag = 0, cflag;
    long long   Numdof = 0; 
    long long Nelem = 0; /*  elements */
    char sn[40];

    FILE *snD, *snD2;
	int maxnlocal, mapnlocal, M, N, info, ZERO, ONE , N0, N1, M1, N2;
	static BOOLEAN export_mesh = FALSE;
	static BOOLEAN export_D11 = TRUE;
	FLOAT localnorm, globalnorm, D11, err;
   
	ZERO =0; ONE =1;
	snD = fopen("APODD11_0.5.text", "w");
	snD2 = fopen("APODerr.text", "w");

    phgInit(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	g = phgNewGrid(-1);

	phgSetPeriodicity(g, X_MASK | Y_MASK | Z_MASK);
    if (!phgImport(g, fn, TRUE))
        phgError(1, "can't read file \"%s\".\n", fn);
	phgRefineAllElements(g, RefineN);
	
	u_p = phgDofNew(g, DOF_P1, 1, "u_p", func_u0);
    u_h = phgDofNew(g, DOF_P1, 1, "u_h", DofInterpolation);
    f_h = phgDofNew(g, DOF_P1, 1, "f_h", DofInterpolation);
	B1=phgDofNew(g, DOF_P1, 1, "B1", func_B1);
	B2=phgDofNew(g, DOF_P1, 1, "B2", func_B2);
	B3=phgDofNew(g, DOF_P1, 1, "B3", func_B3);
	C1=phgDofNew(g, DOF_P1, 1, "C1", func_C1);
	C2=phgDofNew(g, DOF_P1, 1, "C2", func_C2);
	C3=phgDofNew(g, DOF_P1, 1, "C3", func_C3);
    
	map = phgMapCreate(u_h, NULL);
	G1 = phgMapCreateMat(map, map);
	G2 = phgMapCreateMat(map, map);
	CosG1 = phgMapCreateMat(map, map);
	CosG2 = phgMapCreateMat(map, map);
	C = phgMapCreateMat(map, map);

	build_rhs_Mat(C, map, u_h);
    build_stiff_Mat1(stime, D0, G1, map, B1, B2, B3, u_h);
    build_stiff_Mat2(stime, G2, map, C1, C2, C3, u_h);
    build_stiff_Mat1(stime1, D0, CosG1, map, B1, B2, B3, u_h);
    build_stiff_Mat2(stime1, CosG2, map, C1, C2, C3, u_h);


    Nelem = g->nleaf_global;
    Numdof = DofGetDataCountGlobal(u_h);
     
    nT = PODT0 / stime + 1;
    TotalN = T / stime + 1;
	TotalN1 = T / stime1 + 1;
    numt1 = numt * stime1/stime;
    
    uhs = phgAlloc(TotalN * sizeof(*uhs));
    coarseuhs = phgAlloc(TotalN1 * sizeof(*uhs));

    uhs[0] = phgDofCopy(u_p, NULL, DOF_P1, "u_p1");
    uhs[0]->userfunc = DofInterpolation;
    coarseuhs[0] = phgDofCopy(u_p, NULL, DOF_P1, "u_p1");
    coarseuhs[0]->userfunc = DofInterpolation;
    /*******************computing coarse grid DOF*************************/
	kk = 0;  k = 0;
   phgGetTime(tt);
  while(crtime < T - 1e-8)
  {
   /********************************************************************/	
   	phgGetTime(tt1);
   	ptime = crtime;
    crtime += stime1;
   	phgPrintf("\n/********* start new time layer(coarse grid) *************/\n");
    phgPrintf("current time layer: [%lf, %lf]\n", (FLOAT)ptime, (FLOAT)crtime);
       
    flag++;
     if (flag > 1)
   	{   /* update u_p */
         phgDofFree(&u_p);
         u_p = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
         u_p->userfunc = DofInterpolation;
     }
    b = phgMapCreateVec(map, 1);
   	phgVecDisassemble(b);
   	phgDofSetDataByFunction(f_h, func_f);
   	build_rhs_Vec(stime1, b, f_h, map, u_h);
    phgFEMLoopTime3(crtime, u_h, u_p, C, CosG1, CosG2, b);
    kk++;
   	coarseuhs[kk] = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
   	coarseuhs[kk]->userfunc = DofInterpolation;
   	phgPrintf("\ntime step: %lg\n\n", (FLOAT)stime);
  }
   phgDofFree(&u_p);
   u_p = phgDofNew(g, DOF_P1, 1, "u_p", func_u0);
   /**********************************************************************/
  /* init parameters */
   crtime = 0; flag =0;
   while(crtime < PODT0 - 1e-8)
   {
	/********************************************************************/	
    	phgGetTime(tt1);
		ptime = crtime;
        crtime += stime;
		phgPrintf("\n/********* start new time layer *************/\n");
        phgPrintf("current time layer: [%lf, %lf]\n", (FLOAT)ptime, (FLOAT)crtime);
        
        flag++;
        if (flag > 1)
		{   /* update u_p */
          phgDofFree(&u_p);
          u_p = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
          u_p->userfunc = DofInterpolation;
        }
        b = phgMapCreateVec(map, 1);
		phgVecDisassemble(b);
		phgDofSetDataByFunction(f_h, func_f);
		build_rhs_Vec(stime, b, f_h, map, u_h);
	    phgFEMLoopTime3(crtime, u_h, u_p, C, G1, G2, b);
	
        k++;
		uhs[k] = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
		uhs[k]->userfunc = DofInterpolation;
		phgPrintf("\ntime step: %lg\n\n", (FLOAT)stime);
		if (export_mesh) 
		 {
		    *sn = '\0';
			 sprintf(sn, "../exportdata/T100t0.2Refine4_2pi/a4_%d.vtk", k);
			 phgPrintf("output mesh to %s:\n", sn);
			 phgExportVTK(g, sn, u_h, NULL);
		}
    }

/*************************get Vector from Dofs**************************/
	V = phgMapCreateVec(map, nT);
	phgMapDofArraysToVec(map, nT, FALSE, &V, uhs, NULL);
/******************************SVD decomposition***********************/
	MPI_Reduce(&map->nlocal, &maxnlocal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxnlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
	M = numprocs * maxnlocal;
	N = nT;

    U = (FLOAT *)malloc((maxnlocal * M ) * sizeof(FLOAT));
    U1 = (FLOAT *)malloc((maxnlocal * M ) * sizeof(FLOAT));

    S = (FLOAT *)malloc(N*sizeof(FLOAT));
    phgMatSvd(V->data, map, nT, U, S, numprocs, rank);
	for(i=0; i< N; i++)
		phgPrintf("%f\n", S[i]);
///*************************************************************/
	N0 = phgGetPodNumbasis(gamma1, S, nT);
	free(S);
/*******************************POD********************************************/
/**************assemble POD basis matrix*******************/ 
    FLOAT *PODmat, *PODC, *PODu_p;
    VEC *PODu0, *VL, *VR, *Vtmp, *Vtf;
	FLOAT *PODG1, *PODG2, *PODb, *y;
    int  PODflag = 0, *IPIV;
	double  done, dzero;
    done = 1.0;  dzero = 0.0;

	Vtmp = phgMapCreateVec(map, 1);
	PODu0 = phgMapCreateVec(map, 1);
    VL = phgMapCreateVec(map, 1);
	VR = phgMapCreateVec(map, 1);
/**********start  loop time********/
   tmpt =0 ; PODflag = 0;
   phgPrintf("crtime: %f\n", crtime);
   while(crtime < T - 1e-8)
   {
    if(tmpt == 0 || tmpt == numt1)
	{
		if(tmpt == numt1)
	    {
			free(PODmat);  free(PODC); free(PODu_p);free(y);  free(IPIV);
			free(PODb);  free(PODG1); free(PODG2);
		}
//*********************SVD1*************************/
		if((T - crtime) / stime1 < numt)
			numt = (T - crtime)/stime1;
	    V1 = phgMapCreateVec(map, numt);
		cflag = crtime/stime1 + 1;
		phgPrintf("cflag: %d\n", cflag);
		phgMapDofArraysToVec(map, numt, FALSE, &V1, &coarseuhs[cflag], NULL);
		S1 = (FLOAT *)malloc(numt * sizeof(FLOAT));
        phgMatSvd(V1->data, map, numt, U1, S1, numprocs, rank);
		N1 = phgGetPodNumbasis(gamma2, S1, numt);
/*******************SVD2*************************/
		svdM =(FLOAT *) malloc(map->nlocal*(N0 + N1)*sizeof(FLOAT));
        for(i=0; i< N0 + N1; i++)
		{
			if(i< N0)
			memcpy(&svdM[i*map->nlocal], &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
	        if(i >= N0)
			memcpy(&svdM[i*map->nlocal], &U1[(i-N0)*maxnlocal], map->nlocal*sizeof(FLOAT));
		}
    	MPI_Barrier(MPI_COMM_WORLD);
		S = (FLOAT *)malloc((N0 + N1) * sizeof(FLOAT));
		phgMatSvd(svdM, map, N0+N1, U, S, numprocs, rank);
		N0 = phgGetPodNumbasis(gamma3, S, N0 + N1);
		free(S1);  free(S); free(svdM);
/*****************************************************/
	    PODmat = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		PODC = (FLOAT *)malloc(N0*N0 * sizeof(FLOAT));
		PODu_p = (FLOAT *)malloc(N0 * sizeof(FLOAT));
    
		PODG1 = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		PODG2 = (FLOAT *)malloc(N0*N0*sizeof(FLOAT));
		PODb = (FLOAT *)malloc(N0 * sizeof(FLOAT));
		y = (FLOAT *)malloc(N0 * sizeof(FLOAT));
		IPIV = (int *)malloc(N0 * sizeof(int));
/************build POD linear system***********************/	
    	phgMapDofArraysToVec(map, 1, FALSE, &PODu0, &uhs[k], NULL);
    	for(i=0; i< N0; i++)
    	{
	    	memcpy(VL->data, &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
	    	for(j=0; j< N0; j++)
	        {
		       memcpy(VR->data, &U[j*maxnlocal], map->nlocal*sizeof(FLOAT));
			   phgMatVec(0, 1.0, G1, VR, 0.0, &Vtmp);
	 	       PODG1[j*N0 + i] = phgVecDot(Vtmp, 0, VL, 0, NULL);
		       phgMatVec(0, 1.0, G2, VR, 0.0, &Vtmp);
		       PODG2[j*N0 + i] = phgVecDot(Vtmp, 0, VL, 0, NULL);
		       phgMatVec(0, 1.0, C, VR, 0.0, &Vtmp);
		       PODC[j*N0 + i]  = phgVecDot(Vtmp, 0, VL, 0, NULL);
		    }
           PODu_p[i] = phgVecDot(VL, 0, PODu0, 0, NULL);
  	   }
		tmpt = 1;
		PODflag = 0;
	}
	ptime = crtime;
    crtime += stime;
    phgPrintf("current time layer: [%lf, %lf]\n", (FLOAT)ptime, (FLOAT)crtime);
	if(PODflag > 0)
		memcpy(PODu_p, PODb, N0*sizeof(FLOAT));

    b = phgMapCreateVec(map, 1);
	phgVecDisassemble(b);
	phgDofSetDataByFunction(f_h, func_f);
	build_rhs_Vec(stime, b, f_h, map, u_h);
/**********************************************************/
    for(i=0; i < N0; i++)
    {    for(j=0 ;j< N0; j++)
	    	PODmat[i*N0 + j] = PODG1[i*N0 + j] + eta_t(crtime)*PODG2[i*N0 + j];
	    	memcpy(VL->data, &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
            PODb[i] = phgVecDot(VL, 0, b, 0, NULL);
	} 

	dgemv_("T", &N0, &N0, &done, PODC, &N0, PODu_p, &ONE, &dzero, y, &ONE);
    
	for(i=0; i< N0; i++)
		PODb[i] += y[i];
  
    dgesv_(&N0, &ONE, PODmat, &N0, IPIV, PODb, &N0, &info);

	phgPrintf("info :%d\n", info);
    /***************transform PODb to DOF****************/
	phgPODDofToFEMDof(map, N0, U, maxnlocal, PODb, u_h);
	if (export_mesh) 
	{
		*sn = '\0';
		sprintf(sn, "../exportdata/T100t0.2Refine4_2pi/a4_%d.vtk", k);
		phgPrintf("output mesh to %s:\n", sn);
		phgExportVTK(g, sn, u_h, NULL);
	}
	k++;
    uhs[k] = phgDofCopy(u_h, NULL, DOF_P1, "u_p1");
    uhs[k]->userfunc = DofInterpolation;
    tmpt++;
	PODflag++;
 }
 /*****************computing D11*******************************/
    phgPrintf("N0:%d\n", N0);
    phgGetTime(tt1);
	phgPrintf("%d\t%d\n", TotalN, k);
	phgPrintf("\n Total time: %f \n ", (FLOAT)(tt1[2]-tt[2]));
 //  for(i = 0; i< TotalN; i++)
	{
       D11 = D0*(1 + GetD11(uhs[i]));
	   if(rank==0 && export_D11)
	     fprintf(snD,"%f\t%f\n",  i*stime, D11);
	}
/*******************************computing err***********************************/
//    phgRefineAllElements(g,2);
	phgRefineFEMLoopTime(stime, g, TotalN, uhs, rank, snD2);
/**************************************************************/
  phgPrintf("\n DOF:%"dFMT", elements:%"dFMT"\n", 
  		DofGetDataCountGlobal(u_h), g->nleaf_global);

  phgPrintf("\n Total time: %f \n ", (FLOAT)(tt1[2]-tt[2]));
	free(PODmat);  free(PODC); free(PODu_p); 
    free(U); free(U1); free(PODb);
  phgMatDestroy(&G1);
  phgMatDestroy(&G2);
	for(i=0; i< TotalN; i++)
	   phgDofFree(uhs+i);
  for(i=0; i< TotalN1; i++)
     phgDofFree(coarseuhs+i);
   phgFree(uhs);
	phgFree(coarseuhs);
	phgDofFree(&B1);
    phgDofFree(&B2);
    phgDofFree(&B3);
    phgDofFree(&C1);
    phgDofFree(&C2);
    phgDofFree(&C3);
    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&u_p);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
