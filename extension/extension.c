#include <stdlib.h>
#include <string.h>
#include "lapack.h"
#include "libguile.h"

int
eigendecomp_(char jobz, char range, char uplo, int n,
	     double *a, int lda, double vl, double vu, int il,
	     int iu, double abstol, int m, double *w, double *z, int ldz)
{
	int isuppz[n * 2];
	int info;
	int lwork = -1, liwork = -1;

	double work_temp[1];
	int iwork_temp[1];
	dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu,
		&abstol, &m, w, z, &ldz, isuppz, work_temp, &lwork,
		iwork_temp, &liwork, &info);
	if (info)
		return info;
	lwork = *work_temp;
	liwork = *iwork_temp;
	double work[lwork];
	int iwork[liwork];
	dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu,
		&abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork,
		&liwork, &info);
	if (info)
		return info;
	return 0;
}

SCM
eigendecomp(SCM jobz, SCM range, SCM uplo, SCM n,
	    SCM a, SCM vl, SCM vu, SCM il, SCM iu, SCM abstol)
{
	int N = scm_to_int(n);
	double *values = malloc(N * sizeof(double)),
	    *vectors = malloc(N * N * sizeof(double));
	int info = eigendecomp_(scm_to_char(jobz), scm_to_char(range),
				scm_to_char(uplo),
				N, scm_to_pointer(a), N,
				scm_to_double(vl), scm_to_double(vu),
				scm_to_int(il), scm_to_int(iu),
				scm_to_double(abstol), 0,
				values, vectors, N);
	SCM vals[3] = { scm_from_pointer(values, NULL),
		scm_from_pointer(vectors, NULL),
		scm_from_int(info)
	};

	return
	    scm_c_values(memcpy(malloc(sizeof(vals)), vals, sizeof(vals)),
			 3);
}

void
init_eigendecomp(void)
{
	scm_c_define_gsubr("eigendecomp", 10, 0, 0, eigendecomp);
}
