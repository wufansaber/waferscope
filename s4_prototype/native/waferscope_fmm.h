/* Project adapter for the pinned GPL S4 runtime. See generated build patch. */
#ifndef WAFERSCOPE_FMM_H
#define WAFERSCOPE_FMM_H
#include <config.h>
#include <complex>
#include <vector>
#include <cmath>
#include "S4.h"
#include "RNP/TBLAS.h"
#ifdef HAVE_BLAS
#include "RNP/TBLAS_ext.h"
#endif
#include "RNP/LinearSolve.h"

// Return the interface-normal axis, or -1 for genuinely 2D patterns.
static int WS4StripeAxis(const Simulation *S, const Layer *L){
    if(S->Lr[0] <= 0 || S->Lr[3] <= 0 || S->Lr[1] != 0 || S->Lr[2] != 0 ||
       L->pattern.nshapes == 0){ return -1; }
    const Material *background = Simulation_GetMaterialByName(S, L->material, NULL);
    if(background->type != 0){ return -1; }
    bool full[2] = {true, true};
    const double period[2] = {S->Lr[0], S->Lr[3]};
    for(int i = 0; i < L->pattern.nshapes; ++i){
        const shape &a = L->pattern.shapes[i];
        if(a.type != RECTANGLE || a.angle != 0 || L->pattern.parent[i] >= 0 ||
           Simulation_GetMaterialByIndex(S, a.tag)->type != 0){ return -1; }
        for(int axis = 0; axis < 2; ++axis){
            full[axis] = full[axis] &&
                std::abs(2*a.vtab.rectangle.halfwidth[axis]-period[axis]) < 1e-12*period[axis];
        }
    }
    const int normal = full[0] ? 1 : (full[1] ? 0 : -1);
    if(normal < 0){ return -1; }
    const double p = period[normal];
    for(int i = 0; i < L->pattern.nshapes; ++i){
        const shape &a = L->pattern.shapes[i];
        const double wa = a.vtab.rectangle.halfwidth[normal];
        if(wa <= 0 || 2*wa > p*(1+1e-12)){ return -1; }
        for(int j = 0; j < i; ++j){
            const shape &b = L->pattern.shapes[j];
            double distance = std::fmod(std::abs(a.center[normal]-b.center[normal]), p);
            distance = std::fmin(distance, p-distance);
            if(distance < wa+b.vtab.rectangle.halfwidth[normal]-1e-12*p){ return -1; }
        }
    }
    return normal;
}

// Apply the 1D Li inverse rule to a stripe layer in a fully 2D modal basis.
// No symmetry reduction is applied to other layers or the scattering matrix.
static int WS4StripeMaterial(const Simulation *S, const Layer *L, int n,
        std::complex<double> *epsilon2, std::complex<double> *epsilon_inv){
    const int normal = WS4StripeAxis(S, L);
    if(normal < 0){ return 0; }
    typedef std::complex<double> Complex;
    const int nn = n*n, n2 = 2*n;
    std::vector<double> values(2*(1+L->pattern.nshapes)), inverse(values.size());
    for(int i = -1; i < L->pattern.nshapes; ++i){
        const Material *m = i < 0 ? Simulation_GetMaterialByName(S, L->material, NULL) :
            Simulation_GetMaterialByIndex(S, L->pattern.shapes[i].tag);
        const Complex e(m->eps.s[0], m->eps.s[1]);
        if(std::abs(e) == 0){ return -1; }
        values[2*(i+1)] = e.real(); values[2*(i+1)+1] = e.imag();
        inverse[2*(i+1)] = (1./e).real(); inverse[2*(i+1)+1] = (1./e).imag();
    }
    std::vector<Complex> epsilon(nn), eta(nn), scratch(nn);
    const int *g = S->solution->G;
    for(int j = 0; j < n; ++j){
        for(int i = 0; i < n; ++i){
            if(g[2*i+1-normal] != g[2*j+1-normal]){ continue; }
            const int dx = g[2*i]-g[2*j], dy = g[2*i+1]-g[2*j+1];
            const double f[2] = {dx*S->Lk[0]+dy*S->Lk[2], dx*S->Lk[1]+dy*S->Lk[3]};
            double ft[2];
            Pattern_GetFourierTransform(&L->pattern, &values[0], f, 2, Simulation_GetUnitCellSize(S), ft);
            epsilon[i+j*n] = Complex(ft[0], ft[1]);
            Pattern_GetFourierTransform(&L->pattern, &inverse[0], f, 2, Simulation_GetUnitCellSize(S), ft);
            eta[i+j*n] = Complex(ft[0], ft[1]);
        }
    }
    RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., epsilon2,n2);
    // S4 uses (Hx,Hy), so Ey's inverse rule occupies the Hx block.
    const int inverse_offset = normal == 1 ? 0 : n+n*n2;
    const int direct_offset = normal == 1 ? n+n*n2 : 0;
    RNP::TBLAS::CopyMatrix<'A'>(n,n, &epsilon[0],n, epsilon2+direct_offset,n2);
    RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., epsilon2+inverse_offset,n2);
    int info = 0;
    RNP::LinearSolve<'N'>(n,n, &eta[0],n, epsilon2+inverse_offset,n2, &info,NULL);
    if(info != 0){ return -1; }
    scratch = epsilon;
    RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., epsilon_inv,n);
    RNP::LinearSolve<'N'>(n,n, &scratch[0],n, epsilon_inv,n, &info,NULL);
    return info == 0 ? 1 : -1;
}
#endif
