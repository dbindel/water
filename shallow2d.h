#ifndef SHALLOW2D_H
#define SHALLOW2D_H

#include <cmath>
#include <array>
#include <algorithm>

#ifndef RESTRICT
#define RESTRICT
#endif

//ldoc on
/**
 * # Shallow water equations
 * 
 * ## Physics picture
 * 
 * The shallow water equations treat water as incompressible and
 * inviscid, and assume that the horizontal velocity remains constant
 * in any vertical column of water.  The unknowns at each point are
 * the water height and the total horizontal momentum in a water
 * column; the equations describe conservation of mass (fluid is
 * neither created nor destroyed) and conservation of linear momentum.
 * We will solve these equations with a numerical method that also
 * exactly conserves mass and momentum (up to rounding error), though
 * it only approximately conserves energy.
 * 
 * The basic variables are water height ($h$), and the velocity components
 * ($u, v$).  We write the governing equations in the form
 * $$
 *   U_t = F(U)_x + G(U)_y
 * $$
 * where
 * $$
 *   U = \begin{bmatrix} h \\ hu \\ hv \end{bmatrix},
 *   F = \begin{bmatrix} hu \\ h^2 u + gh^2/2 \\ huv \end{bmatrix}
 *   G = \begin{bmatrix} hv \\ huv \\ h^2 v + gh^2/2 \end{bmatrix}
 * $$
 * The functions $F$ and $G$ are called *fluxes*, and describe how the
 * conserved quantities (volume and momentum) enter and exit a region
 * of space.
 * 
 * Note that we also need a bound on the characteristic wave speeds
 * for the problem in order to ensure that our method doesn't explode;
 * we use this to control the Courant-Friedrichs-Levy (CFL) number
 * relating wave speeds, time steps, and space steps.  For the shallow
 * water equations, the characteristic wave speed is $\sqrt{g h}$
 * where $g$ is the gravitational constant and $h$ is the height of the
 * water; in addition, we have to take into account the velocity of
 * the underlying flow.
 * 
 * ## Implementation
 * 
 * Our solver takes advantage of C++ templates to get (potentially)
 * good performance while keeping a clean abstraction between the
 * solver code and the details of the physics.  The `Shallow2D`
 * class specifies the precision of the comptutation (single precision),
 * the data type used to represent vectors of unknowns and fluxes
 * (the C++ `std::array`).  We are really only using the class as 
 * name space; we never create an instance of type `Shallow2D`,
 * and the `flux` and `wave_speed` functions needed by the solver are
 * declared as static (and inline, in the hopes of getting the compiler
 * to optimize for us).
 */

struct Shallow2D {

    // Type parameters for solver
    static constexpr int nfield = 3;
    typedef float real;
    typedef std::array<real,nfield> vec;

    // Gravitational force (compile time constant)
    static constexpr real g = 9.8;

    // Compute shallow water fluxes F(U), G(U)
    static void flux(real* FU, real* GU, const real* U,
                     int ncell, int field_stride) {

        real* RESTRICT fh  = FU;
        real* RESTRICT fhu = FU +   field_stride;
        real* RESTRICT fhv = FU + 2*field_stride;

        real* RESTRICT gh  = GU;
        real* RESTRICT ghu = GU +   field_stride;
        real* RESTRICT ghv = GU + 2*field_stride;

        const real* RESTRICT h  = U;
        const real* RESTRICT hu = U +   field_stride;
        const real* RESTRICT hv = U + 2*field_stride;

        std::copy(hu, hu+ncell, fh);
        std::copy(hv, hv+ncell, gh);

        for (int i = 0; i < ncell; ++i) {
            float hi = h[i], hui = hu[i], hvi = hv[i];
            float inv_h = 1/hi;
            fhu[i] = hui*hui*inv_h + (0.5*g)*hi*hi;
            fhv[i] = hui*hvi*inv_h;
            ghu[i] = hui*hvi*inv_h;
            ghv[i] = hvi*hvi*inv_h + (0.5*g)*hi*hi;
        }
    }

    // Compute shallow water wave speed
    static void wave_speed(real* cxy, const real* U,
                           int ncell, int field_stride) {
        using namespace std;

        const real* RESTRICT h  = U;
        const real* RESTRICT hu = U +   field_stride;
        const real* RESTRICT hv = U + 2*field_stride;

        real cx = cxy[0];
        real cy = cxy[1];
        for (int i = 0; i < ncell; ++i) {
            real hi = h[i], hui = hu[i], hvi = hv[i];
            real root_gh = sqrt(g * hi);
            cx = max(cx, abs(hui/hi) + root_gh);
            cy = max(cy, abs(hvi/hi) + root_gh);
        }
        cxy[0] = cx;
        cxy[1] = cy;
    }
};

//ldoc off
#undef RESTRICT
#endif /* SHALLOW2D_H */
