/*        Generated by TAPENADE     (INRIA, Tropics team)
    Tapenade 3.10 (r5498) - 20 Jan 2015 09:48
*/
#include "gmm_dv.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
/*  Hint: NBDirsMax should be the maximum number of differentiation directions
*/

/*
  Differentiation of arr_max in forward (tangent) mode:
   variations   of useful results: arr_max
   with respect to varying inputs: *x
   Plus diff mem management of: x:in
*/
void arr_max_dv(int n, double *x, double (*xd)[NBDirsMax], double *arr_max, 
        double arr_maxd[NBDirsMax], int nbdirs) {
    double m;
    double md[NBDirsMax];
    int i;
    int nd;
    for (nd = 0; nd < nbdirs; ++nd)
        md[nd] = xd[0][nd];
    m = x[0];
    for (i = 1; i < n; ++i)
        if (x[i] > m) {
            for (nd = 0; nd < nbdirs; ++nd)
                md[nd] = xd[i][nd];
            m = x[i];
        }
    *arr_max = m;
    for (nd = 0; nd < nbdirs; ++nd)
        arr_maxd[nd] = md[nd];
}
/*  Hint: NBDirsMax should be the maximum number of differentiation directions
*/


/*
  Differentiation of logsumexp in forward (tangent) mode:
   variations   of useful results: logsumexp
   with respect to varying inputs: *x
   Plus diff mem management of: x:in
*/
void logsumexp_dv(int n, double *x, double (*xd)[NBDirsMax], double *logsumexp
        , double logsumexpd[NBDirsMax], int nbdirs) {
    int i;
    double mx, semx;
    double mxd[NBDirsMax], semxd[NBDirsMax];
    double arg1;
    double arg1d[NBDirsMax];
    int nd;
    arr_max_dv(n, x, xd, &mx, mxd, nbdirs);
    semx = 0.;
    for (nd = 0; nd < nbdirs; ++nd)
        semxd[nd] = 0.0;
    for (i = 0; i < n; ++i) {
        arg1 = x[i] - mx;
        for (nd = 0; nd < nbdirs; ++nd) {
            arg1d[nd] = xd[i][nd] - mxd[nd];
            semxd[nd] = semxd[nd] + arg1d[nd]*exp(arg1);
        }
        semx += exp(x[i] - mx);
    }
    *logsumexp = log(semx) + mx;
    for (nd = 0; nd < nbdirs; ++nd)
        logsumexpd[nd] = semxd[nd]/semx + mxd[nd];
}
/*  Hint: NBDirsMax should be the maximum number of differentiation directions
*/


/*
  Differentiation of log_wishart_prior in forward (tangent) mode:
   variations   of useful results: log_wishart_prior
   with respect to varying inputs: *icf
   Plus diff mem management of: icf:in
*/
void log_wishart_prior_dv(int p, int k, Wishart wishart, double *icf, double (
        *icfd)[NBDirsMax], double *log_wishart_prior, double 
        log_wishart_priord[NBDirsMax], int nbdirs) {
    int n, ik, i, icf_sz;
    double out, C, frobenius, sum_log_diag, tmp;
    double outd[NBDirsMax], frobeniusd[NBDirsMax], sum_log_diagd[NBDirsMax], 
    tmpd[NBDirsMax];
    float arg1;
    double result1;
    int nd;
    n = p + wishart.m + 1;
    icf_sz = p*(p+1)/2;
    arg1 = 0.5*n;
    result1 = log_gamma_distrib(arg1, p);
    C = n*p*(log(wishart.gamma)-0.5*log(2)) - result1;
    out = 0;
    for (nd = 0; nd < nbdirs; ++nd)
        outd[nd] = 0.0;
    for (ik = 0; ik < k; ++ik) {
        frobenius = 0;
        sum_log_diag = 0;
        for (nd = 0; nd < nbdirs; ++nd) {
            frobeniusd[nd] = 0.0;
            sum_log_diagd[nd] = 0.0;
        }
        for (i = 0; i < p; ++i) {
            tmp = icf[icf_sz*ik + i];
            for (nd = 0; nd < nbdirs; ++nd) {
                tmpd[nd] = icfd[icf_sz*ik + i][nd];
                sum_log_diagd[nd] = sum_log_diagd[nd] + tmpd[nd];
                tmpd[nd] = tmpd[nd]*exp(tmp);
            }
            sum_log_diag = sum_log_diag + tmp;
            tmp = exp(tmp);
            for (nd = 0; nd < nbdirs; ++nd)
                frobeniusd[nd] = frobeniusd[nd] + tmpd[nd]*tmp + tmp*tmpd[nd];
            frobenius = frobenius + tmp*tmp;
        }
        for (i = p; i < icf_sz; ++i) {
            tmp = icf[icf_sz*ik + i];
            for (nd = 0; nd < nbdirs; ++nd) {
                tmpd[nd] = icfd[icf_sz*ik + i][nd];
                frobeniusd[nd] = frobeniusd[nd] + tmpd[nd]*tmp + tmp*tmpd[nd];
            }
            frobenius = frobenius + tmp*tmp;
        }
        for (nd = 0; nd < nbdirs; ++nd)
            outd[nd] = outd[nd] + 0.5*(wishart.gamma*wishart.gamma)*frobeniusd
                [nd] - wishart.m*sum_log_diagd[nd];
        out = out + 0.5*wishart.gamma*wishart.gamma*frobenius - wishart.m*
            sum_log_diag;
    }
    *log_wishart_prior = out - k*C;
    for (nd = 0; nd < nbdirs; ++nd)
        log_wishart_priord[nd] = outd[nd];
}
/*  Hint: NBDirsMax should be the maximum number of differentiation directions
*/


/*
  Differentiation of gmm_objective in forward (tangent) mode:
   variations   of useful results: *err
   with respect to varying inputs: *means *icf *alphas
   RW status of diff variables: *err:out *means:in *icf:in *alphas:in
   Plus diff mem management of: err:in means:in icf:in alphas:in
*/
void gmm_objective_dv(int d, int k, int n, double *alphas, double (*alphasd)[
        NBDirsMax], double *means, double (*meansd)[NBDirsMax], double *icf, 
        double (*icfd)[NBDirsMax], double *x, Wishart wishart, double *err, 
        double (*errd)[NBDirsMax], int nbdirs) {
    int ik, ix, id, i, j, icf_sz, icf_off, Lparamsidx;
    double *lse, *Ldiag, *xcentered, *mahal;
	double(*lsed)[NBDirsMax], (*Ldiagd)[NBDirsMax], 
		(*xcenteredd)[NBDirsMax], (*mahald)[NBDirsMax];
    double sumlog_Ldiag, sqsum_mahal, slse, lse_alphas, CONSTANT;
    double sumlog_Ldiagd[NBDirsMax], sqsum_mahald[NBDirsMax], slsed[NBDirsMax]
    , lse_alphasd[NBDirsMax];
    double result1;
    double result1d[NBDirsMax];
    double result2;
    int nd;
    result1 = sqrt(2*3.14159265359);
    result2 = pow(result1, d);
    CONSTANT = 1/result2;
    icf_sz = d*(d+1)/2;
	lsed = (double(*)[NBDirsMax])malloc(k*sizeof(double)*NBDirsMax);
    lse = (double *)malloc(k*sizeof(double));
	Ldiagd = (double(*)[NBDirsMax])malloc(d*sizeof(double)*NBDirsMax);
    Ldiag = (double *)malloc(d*sizeof(double));
	xcenteredd = (double(*)[NBDirsMax])malloc(d*sizeof(double)*NBDirsMax);
    xcentered = (double *)malloc(d*sizeof(double));
	mahald = (double(*)[NBDirsMax])malloc(d*sizeof(double)*NBDirsMax);
    mahal = (double *)malloc(d*sizeof(double));
    slse = 0.;
	for (id = 0; id < d; id++)
	{
		memset(mahald[id], 0, NBDirsMax * sizeof(double));
		memset(xcenteredd[id], 0, NBDirsMax * sizeof(double));
		memset(Ldiagd[id], 0, NBDirsMax * sizeof(double));
	}
	for (ik = 0; ik < k; ik++)
	{
		memset(lsed[ik], 0, NBDirsMax * sizeof(double));
	}
	memset(slsed, 0, NBDirsMax * sizeof(double));
    for (ix = 0; ix < n; ++ix) {
        for (ik = 0; ik < k; ++ik) {
            icf_off = ik*icf_sz;
            sumlog_Ldiag = 0.;
            for (nd = 0; nd < nbdirs; ++nd)
                sumlog_Ldiagd[nd] = 0.0;
            for (id = 0; id < d; ++id) {
                for (nd = 0; nd < nbdirs; ++nd) {
                    sumlog_Ldiagd[nd] = sumlog_Ldiagd[nd] + icfd[icf_off + id]
                        [nd];
                    Ldiagd[id][nd] = icfd[icf_off+id][nd]*exp(icf[icf_off+id])
                    ;
                }
                sumlog_Ldiag = sumlog_Ldiag + icf[icf_off + id];
                Ldiag[id] = exp(icf[icf_off + id]);
            }
            for (id = 0; id < d; ++id) {
                xcentered[id] = x[ix*d + id] - means[ik*d + id];
                for (nd = 0; nd < nbdirs; ++nd) {
                    xcenteredd[id][nd] = -meansd[ik*d+id][nd];
                    mahald[id][nd] = Ldiagd[id][nd]*xcentered[id] + Ldiag[id]*
                        xcenteredd[id][nd];
                }
                mahal[id] = Ldiag[id]*xcentered[id];
            }
            Lparamsidx = d;
            for (i = 0; i < d; ++i)
                for (j = i+1; j < d; ++j) {
                    for (nd = 0; nd < nbdirs; ++nd)
                        mahald[j][nd] = mahald[j][nd] + icfd[icf_off+
                            Lparamsidx][nd]*xcentered[i] + icf[icf_off+
                            Lparamsidx]*xcenteredd[i][nd];
                    mahal[j] = mahal[j] + icf[icf_off+Lparamsidx]*xcentered[i]
                    ;
                    Lparamsidx = Lparamsidx + 1;
                }
            sqsum_mahal = 0.;
            for (nd = 0; nd < nbdirs; ++nd)
                sqsum_mahald[nd] = 0.0;
            for (id = 0; id < d; ++id) {
                for (nd = 0; nd < nbdirs; ++nd)
                    sqsum_mahald[nd] = sqsum_mahald[nd] + mahald[id][nd]*mahal
                        [id] + mahal[id]*mahald[id][nd];
                sqsum_mahal = sqsum_mahal + mahal[id]*mahal[id];
            }
            for (nd = 0; nd < nbdirs; ++nd)
                lsed[ik][nd] = alphasd[ik][nd] + sumlog_Ldiagd[nd] - 0.5*
                    sqsum_mahald[nd];
            lse[ik] = alphas[ik] + sumlog_Ldiag - 0.5*sqsum_mahal;
        }
        logsumexp_dv(k, lse, lsed, &result1, result1d, nbdirs);
        for (nd = 0; nd < nbdirs; ++nd)
            slsed[nd] = slsed[nd] + result1d[nd];
        slse = slse + result1;
    }
    free(mahald);
    free(mahal);
    free(xcenteredd);
    free(xcentered);
    free(Ldiagd);
    free(Ldiag);
    free(lsed);
    free(lse);
    logsumexp_dv(k, alphas, alphasd, &lse_alphas, lse_alphasd, nbdirs);
    for (nd = 0; nd < nbdirs; ++nd)
        (*errd)[nd] = slsed[nd] - n*lse_alphasd[nd];
    *err = n*log(CONSTANT) + slse - n*lse_alphas;
    log_wishart_prior_dv(d, k, wishart, icf, icfd, &result1, result1d, nbdirs)
    ;
    for (nd = 0; nd < nbdirs; ++nd)
		(*errd)[nd] = (*errd)[nd] + result1d[nd];
    *err = *err + result1;
    // this is here so that tapenade would recognize that means and inv_cov_factors are variables
    *err = *err + (means[0] - means[0] + (icf[0] - icf[0]));
}
