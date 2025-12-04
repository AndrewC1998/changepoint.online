#include <Rcpp.h>
#include <math.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

// user includes
#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <numeric>   // for std::inner_product

using namespace std;
using namespace Rcpp;

// declarations
extern "C" {
SEXP eFastC_delta_online(SEXP Z_, SEXP K_, SEXP delta_, SEXP alpha_, SEXP verbose_, SEXP oldlength_, SEXP newlength_, SEXP cSum_, SEXP DLL_, SEXP DRR_, SEXP DLR_, SEXP Left_, SEXP Right_);
}

double dst(const NumericVector& X, const NumericVector& Y, double alpha);
std::vector<double> MEAN_VAR(const std::vector<double>& X);
double delta_sum(NumericMatrix& X, int a, int b, double alpha);
double dist_X(NumericMatrix& X, double alpha);
double dist_XY(NumericMatrix& X, NumericMatrix& Y, double alpha);
std::vector<std::vector<int> > find_locations(const NumericMatrix& A);

// definitions

SEXP eFastC_delta_online(SEXP Z_, SEXP K_, SEXP delta_, SEXP alpha_, SEXP verbose_, SEXP oldlength_, SEXP newlength_, SEXP cSum_, SEXP DLL_, SEXP DRR_, SEXP DLR_, SEXP Left_, SEXP Right_) {
BEGIN_RCPP

    //Convert SEXP variables to types Rcpp/C++ knows
    int K         = as<int>(K_);
    int delta     = as<int>(delta_);
    double alpha  = as<double>(alpha_);
    bool verbose  = as<bool>(verbose_);
    int oldlength = as<int>(oldlength_);
    int newlength = as<int>(newlength_);
    double deltlengthplus = oldlength + delta;
    NumericMatrix Z(Z_);

    int N = newlength + oldlength; // number of observations

    // Reset K to be maximum number of possible segments
    if (N / (delta + 1) < K)
        K = N / (delta + 1);

    NumericMatrix FF(2, N), A(K, N); // GOF and changepoint location matrices
    NumericVector GOF(K);

    std::fill(FF.begin(), FF.end(), R_NegInf);
    std::fill(A.begin(),  A.end(),  -1);

    // Sum of within sample for left and right segments as well as between segment distances
    double dll = 0.0, drr = 0.0, dlr = 0.0;
    // Counters for number of terms in each sum
    int cll = 0, crr = 0, clr = 0;

    if (verbose)
        Rcout << "Starting pre-processing" << std::endl;

    int minsize   = delta + 1;   // minseglen in PELT
    int oldlength1 = oldlength + 1;

    // Calculations for complete portion of statistics in delta neighborhoods
    NumericVector DLL(DLL_), DRR(DRR_), DLR(DLR_);
    // DRR[s] = within sum for Z[s+1], ..., Z[s+delta]
    // DLL[s] = within sum for Z[s-delta+1], ..., Z[s]
    // DLR[s] = between sum for the sets used to calculate DLL[s] and DRR[s]

    for (int s = delta; s < newlength; ++s) {
        DLL[s + oldlength] = delta_sum(Z, s - delta + 1, s, alpha);
        if (s >= delta) // avoid array out of bounds error
            DRR[s + oldlength - delta] = DLL[s + oldlength];
    }
    if (oldlength > 0) {
        for (int i = oldlength - delta; i < oldlength; ++i) {
            DRR[i] = DRR[oldlength - delta - 1];
        }
        for (int i = oldlength; i < oldlength + delta; ++i) {
            DLL[i] = DLL[oldlength + delta];
        }
    }

    // Calculate DLR array in O(delta*N) time
    NumericMatrix Left(Left_), Right(Right_);
    // Left(i,0)  = sum of dst(Z[i], Z[i-1..i-delta])
    // Right(i,0) = sum of dst(Z[i], Z[i+1..i+delta])
    // Left(i,1)  = sum of dst(Z[i], Z[i-delta..i-2*delta+1])
    // Right(i,1) = sum of dst(Z[i], Z[i+delta..i+2*delta-1])

    for (int i = deltlengthplus; i < N - delta; ++i) {
        int ithpointlr = i - oldlength;
        for (int j1 = ithpointlr - delta; j1 < ithpointlr; ++j1)
            Left(i, 0) += dst(Z(ithpointlr, _), Z(j1, _), alpha);
        for (int j2 = ithpointlr + 1; j2 < ithpointlr + delta + 1; ++j2)
            Right(i, 0) += dst(Z(ithpointlr, _), Z(j2, _), alpha);

        if (i >= 2 * delta - 1 + oldlength)
            for (int j3 = i - oldlength - 2 * delta + 1; j3 < i - oldlength - delta + 1; ++j3)
                Left(i, 1) += dst(Z(ithpointlr, _), Z(j3, _), alpha);
        if (i + 2 * delta - 1 < N)
            for (int j4 = ithpointlr + delta; j4 < ithpointlr + 2 * delta; ++j4)
                Right(i, 1) += dst(Z(ithpointlr, _), Z(j4, _), alpha);
    }

    // cover zero cells in Left and Right
    if (oldlength > 0) {
        // update end of previous
        for (int i = oldlength - delta; i < oldlength; ++i) {
            Left(i, 0)  = Left(oldlength - delta - 1, 0);
            Right(i, 0) = Right(oldlength - delta - 1, 0);
        }
        for (int i = oldlength - 2 * delta - 1; i < oldlength; ++i) {
            Left(i, 1)  = Left(oldlength - 2 * delta - 1, 1);
            Right(i, 1) = Right(oldlength - 2 * delta - 1, 1);
        }
        for (int i = oldlength; i < oldlength + delta; ++i) {
            // beginning of new
            Left(i, 0)  = Left(oldlength + delta, 0);
            Right(i, 0) = Right(oldlength + delta, 0);
            Right(i, 1) = Right(oldlength + delta, 1);
        }
        for (int i = oldlength; i < 2 * delta + oldlength; ++i) {
            Left(i, 1) = Left(2 * delta + oldlength, 1);
        }
    }

    // Update DLR
    if (oldlength == 0) {
        for (int i = 1; i < minsize; ++i)
            for (int j = minsize; j < minsize + delta; ++j)
                DLR[minsize - 1] += dst(Z(i, _), Z(j, _), alpha);

        for (int s = minsize; s < N - delta; ++s) {
            double r1 = Left(s, 0);
            double r2 = Right(s - delta, 1);
            double r3 = dst(Z(s, _), Z(s + delta, _), alpha);
            double a1 = Left(s + delta, 1);
            double a2 = Right(s, 0);
            double a3 = dst(Z(s, _), Z(s - delta, _), alpha);
            DLR[s] = DLR[s - 1] - r1 - r2 - r3 + a1 + a2 + a3;
        }
    } else {
        int minold = oldlength + minsize;
        for (int i = oldlength1; i < minold; ++i) {
            for (int j = minsize; j < minsize + delta; ++j) {
                int ithpoint = i - oldlength;
                DLR[minold - 1] += dst(Z(ithpoint, _), Z(j, _), alpha);
            }
        }
        for (int s = minold; s < N - delta; ++s) {
            double r1 = Left(s, 0);
            double r2 = Right(s - delta, 1);
            double r3 = dst(Z(s - oldlength, _), Z(s - oldlength + delta, _), alpha);
            double a1 = Left(s + delta, 1);
            double a2 = Right(s, 0);
            double a3 = dst(Z(s - oldlength, _), Z(s - oldlength - delta, _), alpha);
            DLR[s] = DLR[s - 1] - r1 - r2 - r3 + a1 + a2 + a3;
        }
        for (int i = oldlength - delta; i < oldlength + delta; ++i) {
            double r1 = Left(i, 0);
            double r2 = Right(i - delta, 1);
            double a1 = Left(i + delta, 1);
            double a2 = Right(i, 0);
            DLR[i] = DLR[i - 1] - r1 - r2 + a1 + a2;
        }
    }

    // cumulative sum of distances for adjacent observations
    NumericVector cSum(cSum_);
    if (oldlength > 0) {
        cSum[oldlength] = cSum[oldlength - 1];
    }
    for (int i = oldlength1; i < N; ++i) {
        int ithpoint = i - oldlength;
        cSum[i] = cSum[i - 1] + dst(Z(ithpoint, _), Z(ithpoint - 1, _), alpha);
    }

    if (verbose)
        Rcout << "Pre-processing complete. Starting optimization." << std::endl;

    // Solve for K = 1
    std::vector<std::set<int> > testPoints(N);

    if (K == 1) {
        int t = N - 1;
        std::set<int> cpSet;
        for (int a = delta; a <= t - minsize; ++a)
            cpSet.insert(a);

        for (std::set<int>::iterator si = cpSet.begin(); si != cpSet.end(); ++si) {
            int s = *si;
            int u = -1;

            cll = crr = delta * (delta - 1) / 2;
            clr = delta * delta;

            dll = DLL[s];
            drr = DRR[s];
            dlr = DLR[s];

            dll += cSum[s - delta];
            drr += (cSum[t] - cSum[s + delta]);
            cll += (s - delta - u);
            crr += (t - s - delta);

            double stat = 2.0 * dlr / (clr + 0.0) - drr / (crr + 0.0) - dll / (cll + 0.0);
            double num  = (s - u) * (t - s) * 1.0;
            double dnom = std::pow(t - u + 0.0, 2.0);
            stat *= (num / dnom);

            if (stat > FF(0, t)) {
                FF(0, t) = stat;
                A(0, t)  = s;
            }
        }
    } else {
        for (int t = 2 * minsize - 1; t < N; ++t) {
            std::set<int> cpSet;
            for (int a = delta; a <= t - minsize; ++a)
                cpSet.insert(a);

            for (std::set<int>::iterator si = cpSet.begin(); si != cpSet.end(); ++si) {
                int s = *si;
                int u = -1;

                cll = crr = delta * (delta - 1) / 2;
                clr = delta * delta;

                dll = DLL[s];
                drr = DRR[s];
                dlr = DLR[s];

                dll += cSum[s - delta];
                drr += (cSum[t] - cSum[s + delta]);
                cll += (s - delta - u);
                crr += (t - s - delta);

                double stat = 2.0 * dlr / (clr + 0.0) - drr / (crr + 0.0) - dll / (cll + 0.0);
                double num  = (s - u) * (t - s) * 1.0;
                double dnom = std::pow(t - u + 0.0, 2.0);
                stat *= (num / dnom);

                if (stat > FF(0, t)) {
                    FF(0, t) = stat;
                    A(0, t)  = s;
                }
            }

            // pruning
            std::set<int> cpSetCopy = cpSet;
            int a2 = t - minsize;
            int b2 = A(0, a2);

            double V2 = (b2 <= 0) ? 0.0 : cSum[b2 - 1];
            double cst2  = (a2 - b2) * (t - a2) / std::pow(t - b2 + 0.0, 2.0);
            double stat2 = 2 * DLR[a2] / (delta * delta);
            double t1_2  = (DRR[a2] + cSum[t] - cSum[a2 - delta]) /
                           (delta * (delta - 1) / 2 + t - a2 - delta);
            double t2_2  = (DLL[a2] + cSum[a2 - delta] - V2) /
                           (delta * (delta - 1) / 2 + a2 - b2 - delta);
            stat2 -= (t1_2 + t2_2);
            stat2 *= cst2;

            for (std::set<int>::iterator si = cpSetCopy.begin(); si != cpSetCopy.end(); ++si) {
                int a = *si;
                int b = A(0, a);
                double V = (b <= 0) ? 0.0 : cSum[b - 1];

                double cst  = (a - b) * (t - a) / std::pow(t - b + 0.0, 2.0);
                double stat = 2 * DLR[a] / (delta * delta);
                double t1   = (DRR[a] + cSum[t] - cSum[a - delta]) /
                              (delta * (delta - 1) / 2 + t - a - delta);
                double t2   = (DLL[a] + cSum[a - delta] - V) /
                              (delta * (delta - 1) / 2 + a - b - delta);
                stat -= (t1 + t2);
                stat *= cst;

                if (FF(0, a) + stat < FF(0, a2) + stat2)
                    cpSet.erase(a);
            }
            testPoints[t] = cpSet;
        }
    }

    GOF[0] = FF(0, N - 1);

    if (K == 1) {
        if (verbose)
            Rcout << "Finished optimization." << std::endl;
        return wrap(List::create(
            _["number"]   = 1,
            _["estimates"]= A(0, N - 1) + 2,
            _["gofM"]     = GOF,
            _["delta"]    = delta,
            _["alpha"]    = alpha,
            _["verbose"]  = verbose,
            _["csum"]     = cSum,
            _["dll"]      = DLL,
            _["dlr"]      = DLR,
            _["drr"]      = DRR,
            _["left"]     = Left,
            _["right"]    = Right
        ));
    }

    // Solve for K > 1
    bool flag = true;

    for (int k = 1; k < K; ++k) {
        if ((k + 1) >= K) {
            int t = N - 1;
            FF(flag, t) = R_NegInf;
            std::set<int> cpSet = testPoints[t];

            for (std::set<int>::iterator si = cpSet.begin(); si != cpSet.end(); ++si) {
                int s = *si;
                int u = A(k - 1, s);

                cll = crr = delta * (delta - 1) / 2;
                clr = delta * delta;

                dll = DLL[s];
                drr = DRR[s];
                dlr = DLR[s];

                dll += cSum[s - delta];
                if (u > 0)
                    dll -= cSum[u - 1];
                drr += (cSum[t] - cSum[s + delta]);
                cll += (s - delta - u);
                crr += (t - s - delta);

                double stat = 2 * dlr / (clr + 0.0) - drr / (crr + 0.0) - dll / (cll + 0.0);
                double num  = (s - u) * (t - s) * 1.0;
                double dnom = std::pow(t - u + 0.0, 2.0);
                stat *= (num / dnom);
                if (u > 0)
                    stat += FF(1 - flag, s);

                if (stat > FF(flag, t)) {
                    FF(flag, t) = stat;
                    A(k, t)     = s;
                }
            }
        } else {
            for (int t = 2 * minsize - 1; t < N; ++t) {
                FF(flag, t) = R_NegInf;
                std::set<int> cpSet = testPoints[t];

                for (std::set<int>::iterator si = cpSet.begin(); si != cpSet.end(); ++si) {
                    int s = *si;
                    int u = A(k - 1, s);

                    cll = crr = delta * (delta - 1) / 2;
                    clr = delta * delta;

                    dll = DLL[s];
                    drr = DRR[s];
                    dlr = DLR[s];

                    dll += cSum[s - delta];
                    if (u > 0)
                        dll -= cSum[u - 1];
                    drr += (cSum[t] - cSum[s + delta]);
                    cll += (s - delta - u);
                    crr += (t - s - delta);

                    double stat = 2 * dlr / (clr + 0.0) - drr / (crr + 0.0) - dll / (cll + 0.0);
                    double num  = (s - u) * (t - s) * 1.0;
                    double dnom = std::pow(t - u + 0.0, 2.0);
                    stat *= (num / dnom);
                    if (u > 0)
                        stat += FF(1 - flag, s);

                    if (stat > FF(flag, t)) {
                        FF(flag, t) = stat;
                        A(k, t)     = s;
                    }
                }

                std::set<int> cpSetCopy = cpSet;
                int a2 = t - minsize;
                int b2 = A(k, a2);

                double V2   = (b2 <= 0) ? 0.0 : cSum[b2 - 1];
                double cst2 = (a2 - b2) * (t - a2) * 1.0 / ((t - b2) * (t - b2));
                double stat2 = 2 * DLR[a2] / (delta * delta);
                double t1_2  = (DRR[a2] + cSum[t] - cSum[a2 - delta]) /
                               (delta * (delta - 1) / 2 + t - a2 - delta);
                double t2_2  = (DLL[a2] + cSum[a2 - delta] - V2) /
                               (delta * (delta - 1) / 2 + a2 - b2 - delta);
                stat2 -= (t1_2 + t2_2);
                stat2 *= cst2;

                for (std::set<int>::iterator si = cpSetCopy.begin(); si != cpSetCopy.end(); ++si) {
                    int a = *si;
                    int b = A(k, a);

                    double V = (b <= 0) ? 0.0 : cSum[b - 1];
                    double cst  = (a - b) * (t - a) * 1.0 / ((t - b) * (t - b));
                    double stat = 2 * DLR[a] / (delta * delta);
                    double t1   = (DRR[a] + cSum[t] - cSum[a - delta]) /
                                  (delta * (delta - 1) / 2 + t - a - delta);
                    double t2   = (DLL[a] + cSum[a - delta] - V) /
                                  (delta * (delta - 1) / 2 + a - b - delta);
                    stat -= (t1 + t2);
                    stat *= cst;

                    if (FF(flag, a) + stat < FF(flag, a2) + stat2)
                        cpSet.erase(a);
                }
                testPoints[t] = cpSet;
            }
        }

        GOF[k] = FF(flag, N - 1);
        flag   = 1 - flag;
    }

    if (verbose)
        Rcout << "Finished Optimization." << std::endl;

    // point of inflection in GOF
    std::vector<double> POI;
    for (int i = 1; i < (GOF.size() - 1); ++i) {
        double x1_mean = 0.0;
        for (int j = 0; j <= i; ++j) x1_mean += j;
        x1_mean /= (i + 1);

        double y1_mean = 0.0;
        for (int j = 0; j <= i; ++j) y1_mean += GOF[j];
        y1_mean /= (i + 1);

        double beta1_num = 0.0, beta1_denom = 0.0;
        for (int j = 0; j <= i; ++j) {
            beta1_num   += (j - x1_mean) * (GOF[j] - y1_mean);
            beta1_denom += std::pow(j - x1_mean, 2);
        }
        double beta1  = beta1_num / beta1_denom;
        double alpha1 = y1_mean - beta1 * x1_mean;

        double x2_mean = 0.0;
        for (int j = i; j < (int)GOF.size(); ++j) x2_mean += j;
        x2_mean /= (GOF.size() - i);

        double y2_mean = 0.0;
        for (int j = i; j < (int)GOF.size(); ++j) y2_mean += GOF[j];
        y2_mean /= (GOF.size() - i);

        double beta2_num = 0.0, beta2_denom = 0.0;
        for (int j = i; j < (int)GOF.size(); ++j) {
            beta2_num   += (j - x2_mean) * (GOF[j] - y2_mean);
            beta2_denom += std::pow(j - x2_mean, 2);
        }
        double beta2  = beta2_num / beta2_denom;
        double alpha2 = y2_mean - beta2 * x2_mean;

        double SSE = 0.0;
        for (int j = 0; j <= i; ++j)
            SSE += std::pow(alpha1 + beta1 * j - GOF[j], 2);
        for (int j = i; j < (int)GOF.size(); ++j)
            SSE += std::pow(alpha2 + beta2 * j - GOF[j], 2);

        POI.push_back(SSE);
    }

    int k = -1;
    double SSE_inf = R_PosInf;
    for (int i = 0; i < (int)POI.size(); ++i) {
        if (POI[i] < SSE_inf) {
            SSE_inf = POI[i];
            k       = i;
        }
    }

    std::vector<std::vector<int> > cpLocs = find_locations(A);

    // FIXED: use lambda instead of deprecated bind2nd
    for (int i = 0; i < K; ++i) {
        std::transform(cpLocs[i].begin(), cpLocs[i].end(),
                       cpLocs[i].begin(),
                       [](int x) { return x + 2; });
    }

    std::vector<int> cps = cpLocs[k + 1];

    return wrap(List::create(
        _["number"]   = 1,
        _["estimates"]= A(0, N - 1) + 2,
        _["gofM"]     = GOF,
        _["cpLoc"]    = cpLocs,
        _["delta"]    = delta,
        _["alpha"]    = alpha,
        _["verbose"]  = verbose,
        _["csum"]     = cSum,
        _["dll"]      = DLL,
        _["dlr"]      = DLR,
        _["drr"]      = DRR,
        _["left"]     = Left,
        _["right"]    = Right
    ));

END_RCPP
}

// sub-functions used within the ecp3o functions.

std::vector<double> MEAN_VAR(const std::vector<double>& X) {
    if (X.empty()) {
        std::vector<double> ret(2, 0.0);
        return ret;
    }

    double mean = *(X.begin()), var = 0.0;

    std::vector<double>::const_iterator i = X.begin();
    ++i;
    int cnt = 2;
    for (; i != X.end(); ++i, ++cnt) {
        double dif = *i - mean;
        var  = var * (cnt - 2) / ((double)cnt - 1) + dif * dif / cnt;
        mean = mean + dif / cnt;
    }

    std::vector<double> res;
    res.push_back(mean);
    res.push_back(var);
    return res;
}

double dst(const NumericVector& X, const NumericVector& Y, double alpha) {
    NumericVector res = X - Y;
    double ip = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
    return std::pow(ip, alpha / 2.0);
}

double delta_sum(NumericMatrix& X, int a, int b, double alpha) {
    double ret = 0.0;
    for (int i = a; i < b; ++i)
        for (int j = i + 1; j <= b; ++j)
            ret += dst(X(i, _), X(j, _), alpha);
    return ret;
}

double dist_X(NumericMatrix& X, double alpha) {
    int n = X.nrow();
    double ret = 0.0;
    for (int i = 0; i < n - 1; ++i)
        for (int j = i + 1; j < n; ++j)
            ret += dst(X(i, _), X(j, _), alpha);
    return ret;
}

double dist_XY(NumericMatrix& X, NumericMatrix& Y, double alpha) {
    int n = X.nrow();
    int m = Y.nrow();
    double ret = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            ret += dst(X(i, _), Y(j, _), alpha);
    return ret;
}

std::vector<std::vector<int> > find_locations(const NumericMatrix& A) {
    std::vector<std::vector<int> > res;
    int K = A.nrow();
    int N = A.ncol();

    for (int k = 0; k < K; ++k) {
        int cp = A(k, N - 1);
        int k1 = k;
        std::vector<int> cps;
        do {
            cps.push_back(cp);
            --k1;
            if (k1 >= 0)
                cp = A(k1, cp);
        } while (k1 >= 0);
        std::sort(cps.begin(), cps.end());
        res.push_back(cps);
    }
    return res;
}
