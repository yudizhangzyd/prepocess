//
//#include <stdio.h>              /* I/O lib         ISOC  */
//#include <stdlib.h>             /* Standard Lib    ISOC  */
//#include <Rembedded.h>
//#include <Rinternals.h>
//
//
//void source(const char *name)
//{
//	SEXP e;
//	PROTECT(e = lang2(install("source"), mkString(name)));
//	R_tryEval(e, R_GlobalEnv, NULL);
//	UNPROTECT(1);
//}
//
//void R_qr(int **mat, int n, int p, int *val, int *length) {
//	SEXP arg;
//	PROTECT(arg = allocMatrix(INTSXP, n, p));
//	for (int i = 0; i < n; ++i) {
//		for (int k = 0; k < p; ++k) {
//			INTEGER(arg)[k + i * p] = mat[i][k];
//		}
//	}
//	
//	SEXP qr_call;
//	PROTECT(qr_call = lang2(install("getNullSpaceCols"), arg));
//	int errorOccurred;
//	SEXP ret = R_tryEval(qr_call, R_GlobalEnv, &errorOccurred);
//	
//	if (!errorOccurred)
//	{
//		int *val = INTEGER(ret);
//		length[0] = LENGTH(ret);
//		printf("R returned: ");
//		for (int i = 0; i < LENGTH(ret); i++)
//		printf("%d", val[i]);
//		printf("\n");
//	}
//	else
//	{
//		printf("Error occurred calling R:getNullSpaceCols\n");
//	}
//	// Unprotect add1_call and arg
//	UNPROTECT(2);
//}
//
//int main(int argc, char **argv) {
//	
//	// Intialize the embedded R environment.
//	int r_argc = 2;
//	char *r_argv[] = { "R", "--silent" };
//	Rf_initEmbeddedR(r_argc, r_argv);
//	
//	int **mat = NULL;
//	mat = malloc(sizeof(*mat) * 2);
//	for (int i = 0; i < 2; i++) {
//		mat[i] = malloc(sizeof(**mat) * 3);
//	}
//	
//	mat[0][0] = 1;
//	mat[0][1] = 2;
//	mat[0][2] = 3;
//	mat[1][0] = 1;
//	mat[1][1] = 2;
//	mat[1][2] = 2;
//	
//	source("qr.R");
//	int **v;
//	int *length;
//	length = malloc(sizeof(*length) * 1);
//	v = malloc(sizeof(*v));
//	R_qr(mat, 2, 3, &v, length);
//	
//	for (int i = 0; i < length[0]; ++i) {
//		printf("%d", &v([i]));
//	}
//	// Release R environment
//	Rf_endEmbeddedR(0);
//	
//	return 0;
//} /* end func main */
