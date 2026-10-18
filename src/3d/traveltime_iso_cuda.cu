//
// © 2026. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by
// Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
// Security Administration. All rights in the program are reserved by
// Triad National Security, LLC, and the U.S. Department of Energy/National
// Nuclear Security Administration. The Government is granted for itself and
// others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works,
// distribute copies to the public, perform publicly and display publicly,
// and to permit others to do so.
//
// Author:
//    Kai Gao, kaigao@lanl.gov
//

//
// CUDA fast sweeping for the 3D factorized eikonal equation and the
// adjoint-state equation, called from module_traveltime_iso_cuda.f90.
//
// The GPU runs the same iterations as fast_sweep_forward and
// fast_sweep_adjoint in module_traveltime_iso.f90: the same eight sweep
// orders, the same level-by-level parallel sweeping (Detrixhe et al., 2013),
// the same local solvers, and the same stopping rule. One kernel updates
// all the points of one level. The points of a level are not neighbors of
// each other, so the result does not depend on the order in which the GPU
// updates them, and it equals that of the sequential sweep.
//
// The local solvers repeat the Fortran expressions operation by operation,
// in double precision. Compiled with --fmad=false, like the Fortran code with
// -fp-model=consistent, the GPU and the CPU give the same traveltime and
// adjoint fields to the last bit.
//

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

// Threads per block of the sweep kernels: 32 consecutive j at one k. Small
// blocks spread the few points of a level over more multiprocessors.
const int block_j = 32;
const int block_k = 1;

// Blocks of the reduction kernel; the host adds their partial sums in a fixed order
const int reduce_threads = 256;
const int reduce_blocks = 512;

//
// Device buffers, kept between calls and enlarged when needed
//
struct buffers {
    size_t n = 0;
    double *a[7] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    char *flag = nullptr;
    double *partial = nullptr;
};

buffers buf;
int device = -1;
char message[512] = "";

// On a failure, keep the message, clear the error and return 1
#define check(call) \
    do { \
        cudaError_t e = (call); \
        if (e != cudaSuccess) { \
            snprintf(message, sizeof(message), "%s (%s:%d)", cudaGetErrorString(e), __FILE__, __LINE__); \
            cudaGetLastError(); \
            return 1; \
        } \
    } while (0)

void release()
{
    for (int l = 0; l < 7; l++) {
        cudaFree(buf.a[l]);
        buf.a[l] = nullptr;
    }
    cudaFree(buf.flag);
    cudaFree(buf.partial);
    buf.flag = nullptr;
    buf.partial = nullptr;
    buf.n = 0;
}

int allocate(size_t n)
{
    for (int l = 0; l < 7; l++) {
        check(cudaMalloc(&buf.a[l], n*sizeof(double)));
    }
    check(cudaMalloc(&buf.flag, n));
    check(cudaMalloc(&buf.partial, reduce_blocks*sizeof(double)));
    buf.n = n;
    return 0;
}

// Make room for n grid points; on a failure, e.g., too little memory, free all
int reserve(size_t n)
{
    if (n <= buf.n) {
        return 0;
    }
    release();
    if (allocate(n) != 0) {
        release();
        return 1;
    }
    return 0;
}

//
// The point of a level that a thread updates. A sweep in the direction
// (sx, sy, sz) visits the points with io + jo + ko = level, where io = i or
// nx - 1 - i for sx = 1 or -1 (0-based indices), and so on. The threads
// cover jo = jlo, ..., jlo + nj - 1 and ko = klo, ..., klo + nk - 1; those
// whose io falls outside the grid have no point.
//
struct level_point {
    int i, j, k;
    bool valid;
};

__device__ level_point point_of_thread(int nx, int ny, int nz, int sx, int sy, int sz,
    int level, int jlo, int klo, int nj, int nk)
{
    level_point p;
    int jo = blockIdx.x*blockDim.x + threadIdx.x;
    int ko = blockIdx.y*blockDim.y + threadIdx.y;
    p.valid = jo < nj && ko < nk;
    jo += jlo;
    ko += klo;
    int io = level - jo - ko;
    p.valid = p.valid && io >= 0 && io < nx;
    p.i = sx > 0 ? io : nx - 1 - io;
    p.j = sy > 0 ? jo : ny - 1 - jo;
    p.k = sz > 0 ? ko : nz - 1 - ko;
    return p;
}

//
// Update the points of one level; see fast_sweep_forward in module_traveltime_iso.f90.
//
// A point takes the minimum of its eight octant solutions (local_solver_3d).
// An octant uses, along each axis, the neighbor on one side of the point:
// side 0 is the neighbor before the point (sign 1 in the Fortran code), and
// side 1 the one after it (sign -1). When the 3D solution of an octant fails,
// the octant falls back to three 2D solutions (local_solver_2d), and a failed
// 2D solution falls back to two 1D solutions. The octants of a point share
// these: each of the twelve 2D solutions serves two octants, and each of the
// six 1D solutions serves four 2D solutions. Here each one runs at most once,
// and only when needed, instead of once for every octant that falls back to
// it. Every number comes from the same operations as in the Fortran code, so
// the result is the same.
//
__global__ void __launch_bounds__(block_j*block_k) sweep_forward_level(
    int nx, int ny, int nz, int sx, int sy, int sz, int level, int jlo, int klo, int nj, int nk,
    double dx, double dy, double dz, double huge_value,
    const double *__restrict__ vp, const double *__restrict__ t0, const double *__restrict__ px0,
    const double *__restrict__ py0, const double *__restrict__ pz0, double *__restrict__ tau)
{
    level_point p = point_of_thread(nx, ny, nz, sx, sy, sz, level, jlo, klo, nj, nk);
    if (!p.valid) {
        return;
    }

    const long long stride[3] = {1, nx, (long long)nx*ny};
    const long long id = p.i + stride[1]*p.j + stride[2]*p.k;
    const int index[3] = {p.i, p.j, p.k};
    const int n[3] = {nx, ny, nz};
    const double h[3] = {dx, dy, dz};

    const double t0c = t0[id];
    const double v = vp[id];
    const double p0[3] = {px0[id], py0[id], pz0[id]};

    // The neighbors: whether each exists, its tau and its time tau*t0, and
    // the products of the coefficients a and b of the local solvers
    bool has[3][2];
    double taun[3][2], tn[3][2], aa[3][2], ab[3][2], bb[3][2];
#pragma unroll
    for (int m = 0; m < 3; m++) {
        has[m][0] = index[m] > 0;
        has[m][1] = index[m] < n[m] - 1;
#pragma unroll
        for (int s = 0; s < 2; s++) {
            taun[m][s] = 0;
            tn[m][s] = 0;
            aa[m][s] = 0;
            ab[m][s] = 0;
            bb[m][s] = 0;
            if (has[m][s]) {
                const double sign = s == 0 ? 1 : -1;
                const long long idn = s == 0 ? id - stride[m] : id + stride[m];
                taun[m][s] = tau[idn];
                tn[m][s] = taun[m][s]*t0[idn];
                const double a = p0[m] + t0c/h[m]*sign;
                const double b = -(taun[m][s]*t0c/h[m]*sign);
                aa[m][s] = a*a;
                ab[m][s] = a*b;
                bb[m][s] = b*b;
            }
        }
    }

    const double vv = 1.0/(v*v);

    // The 3D solutions of the octants, and whether each falls back to 2D
    double r3[2][2][2];
    bool f3[2][2][2];
#pragma unroll
    for (int l = 0; l < 2; l++) {
#pragma unroll
        for (int m = 0; m < 2; m++) {
#pragma unroll
            for (int o = 0; o < 2; o++) {
                r3[l][m][o] = huge_value;
                f3[l][m][o] = false;
                if (has[0][l] && has[1][m] && has[2][o]) {
                    const double a = aa[0][l] + aa[1][m] + aa[2][o];
                    const double b = 2*(ab[0][l] + ab[1][m] + ab[2][o]);
                    const double c = bb[0][l] + bb[1][m] + bb[2][o] - vv;
                    f3[l][m][o] = true;
                    if (b*b - 4*a*c > 0) {
                        const double tau1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
                        const double tau2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
                        const double tmax = fmax(fmax(tn[0][l], tn[1][m]), tn[2][o]);
                        const bool causality1 = tau1*t0c >= tmax;
                        const bool causality2 = tau2*t0c >= tmax;
                        if (causality1 && causality2) {
                            r3[l][m][o] = fmin(tau1, tau2);
                            f3[l][m][o] = false;
                        } else if (causality1) {
                            r3[l][m][o] = tau1;
                            f3[l][m][o] = false;
                        } else if (causality2) {
                            r3[l][m][o] = tau2;
                            f3[l][m][o] = false;
                        }
                    }
                }
            }
        }
    }

    // The 2D solutions on the planes xy, xz and yz, for the sides of the
    // first and the second axis of the plane, and whether each falls back to 1D
    const int plane[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    double r2[3][2][2];
    bool f2[3][2][2];
#pragma unroll
    for (int q = 0; q < 3; q++) {
#pragma unroll
        for (int l = 0; l < 2; l++) {
#pragma unroll
            for (int m = 0; m < 2; m++) {
                const bool need = q == 0 ? f3[l][m][0] || f3[l][m][1]
                    : q == 1 ? f3[l][0][m] || f3[l][1][m] : f3[0][l][m] || f3[1][l][m];
                const int m1 = plane[q][0], m2 = plane[q][1];
                r2[q][l][m] = huge_value;
                f2[q][l][m] = false;
                if (need) {
                    const double a = aa[m1][l] + aa[m2][m];
                    const double b = 2*(ab[m1][l] + ab[m2][m]);
                    const double c = bb[m1][l] + bb[m2][m] - vv;
                    f2[q][l][m] = true;
                    if (b*b - 4*a*c > 0) {
                        const double tau1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
                        const double tau2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
                        const bool causality1 = tau1*t0c >= tn[m1][l] && tau1*t0c >= tn[m2][m];
                        const bool causality2 = tau2*t0c >= tn[m1][l] && tau2*t0c >= tn[m2][m];
                        if (causality1 && causality2) {
                            r2[q][l][m] = fmin(tau1, tau2);
                            f2[q][l][m] = false;
                        } else if (causality1) {
                            r2[q][l][m] = tau1;
                            f2[q][l][m] = false;
                        } else if (causality2) {
                            r2[q][l][m] = tau2;
                            f2[q][l][m] = false;
                        }
                    }
                }
            }
        }
    }

    // The 1D solutions along each axis from the neighbor on each side
    double r1[3][2];
#pragma unroll
    for (int m = 0; m < 3; m++) {
#pragma unroll
        for (int s = 0; s < 2; s++) {
            const bool need = m == 0 ? f2[0][s][0] || f2[0][s][1] || f2[1][s][0] || f2[1][s][1]
                : m == 1 ? f2[0][0][s] || f2[0][1][s] || f2[2][s][0] || f2[2][s][1]
                : f2[1][0][s] || f2[1][1][s] || f2[2][0][s] || f2[2][1][s];
            r1[m][s] = huge_value;
            if (need) {
                const double sign = s == 0 ? 1 : -1;
                r1[m][s] = fmax((t0c*taun[m][s] + h[m]/v)/(t0c + p0[m]*h[m]*sign), tn[m][s]/t0c);
            }
        }
    }

    // The fallbacks
#pragma unroll
    for (int q = 0; q < 3; q++) {
#pragma unroll
        for (int l = 0; l < 2; l++) {
#pragma unroll
            for (int m = 0; m < 2; m++) {
                if (f2[q][l][m]) {
                    r2[q][l][m] = fmin(r1[plane[q][0]][l], r1[plane[q][1]][m]);
                }
            }
        }
    }

    double t = tau[id];
#pragma unroll
    for (int l = 0; l < 2; l++) {
#pragma unroll
        for (int m = 0; m < 2; m++) {
#pragma unroll
            for (int o = 0; o < 2; o++) {
                if (f3[l][m][o]) {
                    r3[l][m][o] = fmin(fmin(r2[0][l][m], r2[1][l][o]), r2[2][m][o]);
                }
                t = fmin(t, r3[l][m][o]);
            }
        }
    }
    tau[id] = t;
}

//
// Update the points of one level; see fast_sweep_adjoint in module_traveltime_iso.f90
//
__global__ void __launch_bounds__(block_j*block_k) sweep_adjoint_level(
    int nx, int ny, int nz, int sx, int sy, int sz, int level, int jlo, int klo, int nj, int nk,
    double d1, double d2, double d3,
    const double *__restrict__ tt, const char *__restrict__ recrflag, const double *__restrict__ src,
    double *__restrict__ lambda)
{
    level_point p = point_of_thread(nx, ny, nz, sx, sy, sz, level, jlo, klo, nj, nk);
    if (!p.valid) {
        return;
    }

    const long long sj = nx;
    const long long sk = (long long)nx*ny;
    const long long c = p.i + sj*p.j + sk*p.k;

    if (recrflag[c]) {
        return;
    }

    // At a boundary, the missing neighbor is the point itself
    const long long i1 = p.i > 0 ? c - 1 : c;
    const long long i2 = p.i < nx - 1 ? c + 1 : c;
    const long long j1 = p.j > 0 ? c - sj : c;
    const long long j2 = p.j < ny - 1 ? c + sj : c;
    const long long k1 = p.k > 0 ? c - sk : c;
    const long long k2 = p.k < nz - 1 ? c + sk : c;

    // Solve equation (A-9) in Taillandier et al. (2009)
    const double ttc = tt[c];
    double ap = (tt[i2] - ttc)/d1;
    double am = (ttc - tt[i1])/d1;

    double bp = (tt[j2] - ttc)/d2;
    double bm = (ttc - tt[j1])/d2;

    double cp = (tt[k2] - ttc)/d3;
    double cm = (ttc - tt[k1])/d3;

    double app = (ap + fabs(ap))/2.0;
    double apm = (ap - fabs(ap))/2.0;

    double amp = (am + fabs(am))/2.0;
    double amm = (am - fabs(am))/2.0;

    double bpp = (bp + fabs(bp))/2.0;
    double bpm = (bp - fabs(bp))/2.0;

    double bmp = (bm + fabs(bm))/2.0;
    double bmm = (bm - fabs(bm))/2.0;

    double cpp = (cp + fabs(cp))/2.0;
    double cpm = (cp - fabs(cp))/2.0;

    double cmp = (cm + fabs(cm))/2.0;
    double cmm = (cm - fabs(cm))/2.0;

    double lhs = (apm - amp)/d1 + (bpm - bmp)/d2 + (cpm - cmp)/d3;
    double rhs = (amm*lambda[i1] - app*lambda[i2])/d1
        + (bmm*lambda[j1] - bpp*lambda[j2])/d2
        + (cmm*lambda[k1] - cpp*lambda[k2])/d3;

    double t;
    if (lhs == 0) {
        t = 0;
    } else if (src) {
        t = (rhs - src[c])/lhs;
    } else {
        t = rhs/lhs;
    }

    lambda[c] = fmin(lambda[c], t);
}

//
// Partial sums of |a - b|, one per block
//
__global__ void __launch_bounds__(reduce_threads) sum_abs_diff(size_t n, const double *__restrict__ a,
    const double *__restrict__ b, double *__restrict__ partial)
{
    __shared__ double s[reduce_threads];

    double sum = 0;
    for (size_t l = blockIdx.x*size_t(blockDim.x) + threadIdx.x; l < n; l += size_t(gridDim.x)*blockDim.x) {
        sum += fabs(a[l] - b[l]);
    }
    s[threadIdx.x] = sum;
    __syncthreads();

    for (int h = blockDim.x/2; h > 0; h /= 2) {
        if (threadIdx.x < h) {
            s[threadIdx.x] += s[threadIdx.x + h];
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        partial[blockIdx.x] = s[0];
    }
}

//
// Mean of |a - b|, as mean(abs(a - b)) in the Fortran code
//
int mean_abs_diff(size_t n, const double *a, const double *b, double *mean)
{
    double partial[reduce_blocks];

    sum_abs_diff<<<reduce_blocks, reduce_threads>>>(n, a, b, buf.partial);
    check(cudaGetLastError());
    check(cudaMemcpy(partial, buf.partial, reduce_blocks*sizeof(double), cudaMemcpyDeviceToHost));

    double sum = 0;
    for (int l = 0; l < reduce_blocks; l++) {
        sum += partial[l];
    }
    *mean = sum/n;
    return 0;
}

//
// The eight sweep directions, in the order of the Fortran code
//
const int sweep_direction[8][3] = {
    {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1},
    {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};

//
// The range of jo and ko that holds the points of a level
//
struct level_range {
    int jlo, klo, nj, nk;
};

level_range range_of_level(int nx, int ny, int nz, int level)
{
    level_range r;
    int klo = level - (nx - 1) - (ny - 1);
    r.klo = klo > 0 ? klo : 0;
    int khi = level < nz - 1 ? level : nz - 1;
    int jlo = level - khi - (nx - 1);
    r.jlo = jlo > 0 ? jlo : 0;
    int jhi = level - r.klo < ny - 1 ? level - r.klo : ny - 1;
    r.nj = jhi - r.jlo + 1;
    r.nk = khi - r.klo + 1;
    return r;
}

} // namespace

//
// Select a GPU, by the rank of the process on its node when several share
// it, and return its index, or -1 when there is none
//
extern "C" int latte_cuda_init(char *name, int len)
{
    if (device >= 0) {
        return device;
    }

    int ndevice = 0;
    if (cudaGetDeviceCount(&ndevice) != cudaSuccess || ndevice == 0) {
        cudaGetLastError();
        return -1;
    }

    int local = 0;
    const char *vars[] = {"MPI_LOCALRANKID", "OMPI_COMM_WORLD_LOCAL_RANK", "MV2_COMM_WORLD_LOCAL_RANK",
        "PMI_LOCAL_RANK", "PALS_LOCAL_RANKID", "SLURM_LOCALID"};
    for (const char *var : vars) {
        const char *value = getenv(var);
        if (value != nullptr) {
            local = atoi(value);
            break;
        }
    }

    int id = local % ndevice;
    cudaDeviceProp prop;
    if (cudaSetDevice(id) != cudaSuccess || cudaGetDeviceProperties(&prop, id) != cudaSuccess) {
        cudaGetLastError();
        return -1;
    }

    // A space-padded Fortran string
    int m = int(strlen(prop.name));
    for (int l = 0; l < len; l++) {
        name[l] = l < m ? prop.name[l] : ' ';
    }

    device = id;
    return device;
}

//
// The message of the last failure, as a space-padded Fortran string
//
extern "C" void latte_cuda_message(char *text, int len)
{
    int m = int(strlen(message));
    for (int l = 0; l < len; l++) {
        text[l] = l < m ? message[l] : ' ';
    }
}

//
// Fast sweeping iterations for the factorized eikonal equation; see forward_iso.
// The arrays are (nx, ny, nz) in Fortran order. On return, tau is the
// multiplicative time field, niter the number of iterations and diff the
// mean change of tau in the last iteration. Returns 0 on success.
//
extern "C" int latte_cuda_forward_3d(int nx, int ny, int nz, double dx, double dy, double dz,
    const double *vp, const double *t0, const double *px0, const double *py0, const double *pz0,
    double *tau, double huge_value, double threshold, int niter_max, int *niter, double *diff)
{
    const size_t n = size_t(nx)*ny*nz;
    if (reserve(n) != 0) {
        return 1;
    }

    double *d_vp = buf.a[0], *d_t0 = buf.a[1], *d_px0 = buf.a[2], *d_py0 = buf.a[3], *d_pz0 = buf.a[4];
    double *d_tau = buf.a[5], *d_prev = buf.a[6];

    check(cudaMemcpy(d_vp, vp, n*sizeof(double), cudaMemcpyHostToDevice));
    check(cudaMemcpy(d_t0, t0, n*sizeof(double), cudaMemcpyHostToDevice));
    check(cudaMemcpy(d_px0, px0, n*sizeof(double), cudaMemcpyHostToDevice));
    check(cudaMemcpy(d_py0, py0, n*sizeof(double), cudaMemcpyHostToDevice));
    check(cudaMemcpy(d_pz0, pz0, n*sizeof(double), cudaMemcpyHostToDevice));
    check(cudaMemcpy(d_tau, tau, n*sizeof(double), cudaMemcpyHostToDevice));

    int iter = 0;
    double ttdiff = huge_value;
    const dim3 block(block_j, block_k);

    while (ttdiff >= threshold && iter < niter_max) {

        check(cudaMemcpyAsync(d_prev, d_tau, n*sizeof(double), cudaMemcpyDeviceToDevice));

        for (int s = 0; s < 8; s++) {
            const int *sd = sweep_direction[s];
            for (int level = 0; level <= nx + ny + nz - 3; level++) {
                level_range r = range_of_level(nx, ny, nz, level);
                const dim3 grid((r.nj + block_j - 1)/block_j, (r.nk + block_k - 1)/block_k);
                sweep_forward_level<<<grid, block>>>(nx, ny, nz, sd[0], sd[1], sd[2], level,
                    r.jlo, r.klo, r.nj, r.nk, dx, dy, dz, huge_value, d_vp, d_t0, d_px0, d_py0, d_pz0, d_tau);
            }
        }
        check(cudaGetLastError());

        if (mean_abs_diff(n, d_tau, d_prev, &ttdiff) != 0) {
            return 1;
        }
        iter++;

    }

    check(cudaMemcpy(tau, d_tau, n*sizeof(double), cudaMemcpyDeviceToHost));

    *niter = iter;
    *diff = ttdiff;
    return 0;
}

//
// Fast sweeping iterations for the adjoint-state equation; see adjoint_iso.
// The arrays are (nx, ny, nz) in Fortran order; recrflag is 1 at the receiver
// points and 0 elsewhere. When use_src is 1, src holds the point sources of
// the adjoint-state equation; otherwise src is not read. On return, lambda is
// the adjoint field, niter the number of iterations and diff the mean change
// of lambda in the last iteration. Returns 0 on success.
//
extern "C" int latte_cuda_adjoint_3d(int nx, int ny, int nz, double dx, double dy, double dz,
    const double *tt, const char *recrflag, const double *src, int use_src, double *lambda,
    double huge_value, double threshold, int niter_max, int *niter, double *diff)
{
    const size_t n = size_t(nx)*ny*nz;
    if (reserve(n) != 0) {
        return 1;
    }

    double *d_tt = buf.a[0], *d_lambda = buf.a[1], *d_prev = buf.a[2];
    double *d_src = use_src ? buf.a[3] : nullptr;
    char *d_flag = buf.flag;

    check(cudaMemcpy(d_tt, tt, n*sizeof(double), cudaMemcpyHostToDevice));
    check(cudaMemcpy(d_flag, recrflag, n, cudaMemcpyHostToDevice));
    check(cudaMemcpy(d_lambda, lambda, n*sizeof(double), cudaMemcpyHostToDevice));
    if (d_src) {
        check(cudaMemcpy(d_src, src, n*sizeof(double), cudaMemcpyHostToDevice));
    }

    int iter = 0;
    double lambda_diff = huge_value;
    const dim3 block(block_j, block_k);

    while (lambda_diff >= threshold && iter < niter_max) {

        check(cudaMemcpyAsync(d_prev, d_lambda, n*sizeof(double), cudaMemcpyDeviceToDevice));

        for (int s = 0; s < 8; s++) {
            const int *sd = sweep_direction[s];
            for (int level = 0; level <= nx + ny + nz - 3; level++) {
                level_range r = range_of_level(nx, ny, nz, level);
                const dim3 grid((r.nj + block_j - 1)/block_j, (r.nk + block_k - 1)/block_k);
                sweep_adjoint_level<<<grid, block>>>(nx, ny, nz, sd[0], sd[1], sd[2], level,
                    r.jlo, r.klo, r.nj, r.nk, dx, dy, dz, d_tt, d_flag, d_src, d_lambda);
            }
        }
        check(cudaGetLastError());

        if (mean_abs_diff(n, d_lambda, d_prev, &lambda_diff) != 0) {
            return 1;
        }
        iter++;

    }

    check(cudaMemcpy(lambda, d_lambda, n*sizeof(double), cudaMemcpyDeviceToHost));

    *niter = iter;
    *diff = lambda_diff;
    return 0;
}
