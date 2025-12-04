#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#define SWAP(a,b)   { int t; t=a; a=b; b=t; }

/* ------------------------------------------------------------------
   Forward declarations
   ------------------------------------------------------------------ */

static double *tmplike = NULL;

/* cost functions (defined in cost_general_functions.c) */
double mll_var(double x, double x2, double x3, int n, double shape);
double mll_mean(double x, double x2, double x3, int n, double shape);
double mll_meanvar(double x, double x2, double x3, int n, double shape);
double mll_meanvar_exp(double x, double x2, double x3, int n, double shape);
double mll_meanvar_gamma(double x, double x2, double x3, int n, double shape);
double mll_meanvar_poisson(double x, double x2, double x3, int n, double shape);

double mbic_var(double x, double x2, double x3, int n, double shape);
double mbic_mean(double x, double x2, double x3, int n, double shape);
double mbic_meanvar(double x, double x2, double x3, int n, double shape);
double mbic_meanvar_exp(double x, double x2, double x3, int n, double shape);
double mbic_meanvar_gamma(double x, double x2, double x3, int n, double shape);
double mbic_meanvar_poisson(double x, double x2, double x3, int n, double shape);

/* utility to find minimum and its index */
void min_which(double *x, int n, double *minout, int *whichout);

/* ------------------------------------------------------------------
   FreePELT: free temporary memory, called from R if desired
   ------------------------------------------------------------------ */

void FreePELT(int *err)
{
    /* Only free if there was no error and memory was allocated */
    if (err != NULL && *err == 0 && tmplike != NULL) {
        free((void *) tmplike);
        tmplike = NULL;
    }
}

/* ------------------------------------------------------------------
   PELT_online: core online PELT recursion
   ------------------------------------------------------------------ */

void PELT_online(char   **cost_func,
                 double *sumstat,        /* summary statistics for series (n+1 rows, 3 cols) */
                 int    *ndone,          /* length already processed */
                 int    *nupdate,        /* length of new data in this update */
                 double *pen,            /* penalty */
                 int    *cptsout,        /* changepoint locations up to n */
                 int    *err,            /* 0 by default, non-zero for error */
                 double *shape,          /* gamma shape parameter (if used) */
                 int    *minseglen,      /* minimum segment length */
                 double *lastchangelike, /* DP costs up to n (to be extended) */
                 int    *lastchangecpts, /* DP backpointers */
                 int    *checklist,      /* candidate last-changepoints */
                 int    *nchecklist)     /* number of candidates in checklist */
{
    /* choose cost function based on string */
    double (*costfunction)(double, double, double, int, double) = NULL;

    if (strcmp(*cost_func, "var.norm") == 0) {
        costfunction = &mll_var;
    } else if (strcmp(*cost_func, "mean.norm") == 0) {
        costfunction = &mll_mean;
    } else if (strcmp(*cost_func, "meanvar.norm") == 0) {
        costfunction = &mll_meanvar;
    } else if (strcmp(*cost_func, "meanvar.exp") == 0) {
        costfunction = &mll_meanvar_exp;
    } else if (strcmp(*cost_func, "meanvar.gamma") == 0) {
        costfunction = &mll_meanvar_gamma;
    } else if (strcmp(*cost_func, "meanvar.poisson") == 0) {
        costfunction = &mll_meanvar_poisson;
    } else if (strcmp(*cost_func, "mean.norm.mbic") == 0) {
        costfunction = &mbic_mean;
    } else if (strcmp(*cost_func, "var.norm.mbic") == 0) {
        costfunction = &mbic_var;
    } else if (strcmp(*cost_func, "meanvar.norm.mbic") == 0) {
        costfunction = &mbic_meanvar;
    } else if (strcmp(*cost_func, "meanvar.exp.mbic") == 0) {
        costfunction = &mbic_meanvar_exp;
    } else if (strcmp(*cost_func, "meanvar.gamma.mbic") == 0) {
        costfunction = &mbic_meanvar_gamma;
    } else if (strcmp(*cost_func, "meanvar.poisson.mbic") == 0) {
        costfunction = &mbic_meanvar_poisson;
    } else {
        /* unknown cost; mark error and return */
        if (err) *err = 5;
        return;
    }

    /* allocate temporary like vector (global pointer so FreePELT can free it) */
    int n      = *ndone + *nupdate;
    int maxlen = *nupdate + *nchecklist + 1;

    tmplike = (double *) calloc((size_t) maxlen, sizeof(double));
    if (tmplike == NULL) {
        if (err) *err = 4;
        return;
    }

    double minout;
    int    tstar;
    int    i, whichout, nchecktmp;

    int min = 2 * (*minseglen);

    /* Initialisation block: if nothing processed yet */
    if (*ndone == 0) {
        /* convention: F(0) = -pen */
        lastchangelike[0] = -(*pen);
        lastchangecpts[0] = 0;

        if (min > *nupdate) {
            min = *nupdate;    /* can't add a change if not enough data */
        }

        /* first-pass costs for segments [1..i], length i */
        for (i = *minseglen; i < min; i++) {
            /* sumstat layout: (n+1) x 3 stored column-major:
               col1: index 0..n
               col2: offset (n+1)
               col3: offset 2*(n+1)
             */
            double x  = *(sumstat + i);
            double x2 = *(sumstat + n + 1 + i);
            double x3 = *(sumstat + n + n + 2 + i);

            lastchangelike[i] = costfunction(x, x2, x3, i, *shape);
        }

        for (i = *minseglen; i < min; i++) {
            lastchangecpts[i] = 0;
        }

        *ndone    = min;
        *nupdate -= min;

        if (*nupdate == 0) {
            /* nothing left to update: free and return */
            free(tmplike);
            tmplike = NULL;
            return;
        }

        *nchecklist      = 2;
        checklist[0]     = 0;
        checklist[1]     = *minseglen;
    }

    /* Main online PELT recursion */
    for (tstar = *ndone; tstar < (n + 1); tstar++) {
        R_CheckUserInterrupt();

        /* compute candidates' costs for changepoint at each checklist entry */
        for (i = 0; i < *nchecklist; i++) {
            int tau = checklist[i];

            double x  = (*(sumstat + tstar)                - *(sumstat + tau));
            double x2 = (*(sumstat + n + 1 + tstar)        - *(sumstat + n + 1 + tau));
            double x3 = (*(sumstat + n + n + 2 + tstar)    - *(sumstat + n + n + 2 + tau));
            int    len = tstar - tau;

            tmplike[i] =
                lastchangelike[tau] +
                costfunction(x, x2, x3, len, *shape) +
                *pen;
        }

        /* choose best candidate */
        min_which(tmplike, *nchecklist, &minout, &whichout);
        lastchangelike[tstar] = minout;
        lastchangecpts[tstar] = checklist[whichout];

        /* update checklist for next iteration */
        nchecktmp = 0;
        for (i = 0; i < *nchecklist; i++) {
            if (tmplike[i] <= lastchangelike[tstar] + *pen) {
                checklist[nchecktmp] = checklist[i];
                nchecktmp++;
            }
        }
        *nchecklist = nchecktmp;

        /* new candidate: tstar - (minseglen - 1) ensures minseglen on next seg */
        checklist[nchecktmp] = tstar - (*minseglen - 1);
        (*nchecklist)++;
    }

    /* reconstruct changepoints: walk back from n to 0 */
    int ncpts = 0;
    int last  = n;
    while (last != 0) {
        cptsout[ncpts] = last;
        last = lastchangecpts[last];
        ncpts++;
    }

    free(tmplike);
    tmplike = NULL;
}
