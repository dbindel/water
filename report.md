# Shallow water code performance tuning

In these notes (timed around the same time as the mid-project report
is due), I sketch how I got to where I am with the shallow water code.

## General process

My OS X laptop has the same general processor architecture (Haswell)
as the main chips on the Totient cluster.  As the initial issues that
I discovered in the code performance mostly had to do with taking full
advantage of the Haswell vector units (as opposed to doing something
clever with cache blocking), I decided to tune on my laptop initially,
sanity-checking on the cluster periodically.  In general, I alternated
between the following steps:

- Use `amplxe-cl` (on the cluster) and Instruments.app (on my laptop)
  to gather timing information about which parts of the computation
  took the most time.
- Use vectorization reports from ICC and CLang to determine which
  functions and loops were not vectorized (and were important).
- Fix minor performance issues based on measurements.
- More infrequently, make structural changes to computation order to
  improve locality, memory access, or the amount of computation.

I did try MAQAO as well, at least initially.  It mostly told me that I
ought to vectorize, which I already knew.

I was surprised at how good the Clang compiler is at vectorization
(though it is not as good as the Intel compiler).  I was also
surprised at how difficult it was to convince the various compilers
that I used about lack of aliasing in the C++ code.

My code went through several stages on the way to the current
implementation.

## Baseline

The baseline code uses a Lua interface, but otherwise is very similar
to the C++ reference code distributed in class.  For the tiny
dam-break example, the computations for each frame take around 50-60
ms (with more of the latter than the former).  The Intel compiler
refuses to vectorize most loops, with a few exceptions.  In many
cases, the message from the compiler is something like

    remark #15523: loop was not vectorized: loop control variable ix was found, but loop iteration count cannot be computed before executing the loop

Invariably, the loop bounds in such cases turned out to be simple
functions of `const` members of the class.  From this, I can only
infer that the compiler is being conservative about the possibility
that an insane programmer can force a `const` member to be modified.
It is far easier to establish that there is no aliasing of local
variables declared inside a function or passed into a function as
arguments, which suggested one set of modifications to attempt.  I
also noticed several places where the compiler assumed a flow
(read-after-write) or anti (write-after-read) dependency.
Unfortunately, the compiler doesn't always give easy-to-understand
pointers to where these dependencies are, but looking around the
lines of code mentioned showed me that the compiler was assuming
aliasing unless otherwise stated.

In an ideal world, I would be able to tell the C++ compiler that
there is no aliasing across my `std::vector` objects.  Either that, or
`std::valarray` would have evolved with the C++ standard instead of
being created and almost immediately abandoned.  As it is, the
simplest way to convince the compiler that there is no aliasing seemed
to be to move the performance critical parts into C-style code, where
the size parameters used to compute loop bounds are passed in as
explicit arguments, along with (`restrict`) pointers to any data.

## array4vec

My initial attempt at vectorizing the code was too ambitious: I went
directly from the original C++ code to a C code.  I did not take
incremental steps, and so was absolutely stumped about where to look
when I saw my code fail.  Lesson: change the code incrementally, even
if you think you know where you're headed, if only so that you can
test as you go.

## Devec

One of the first things in the Intel performance optimization
reference discussing vectorization is a discussion of
"array-of-struct" versus "struct-of-array" data structures.  The C++
code follows the "array-of-struct" style, which is great for memory
locality but lousy for vectorization.  The "struct-of-array" (aka
parallel array) style is in some ways less natural for high-level C++
programming, since the natural mechanism for data abstraction (an object)
typically involves storing everything together in one place.

The `devec` branch of the code switches to a "struct-of-array" style.
It doesn't really make things much faster on its own; worse, I
introduced a bug around this time that broke conservation of momentum.
Nonetheless, it was a useful development stage.

In order to use the parallel array data structure and maintain an
abstraction barrier between the physics and the solver, I forced the
arrays for the different fields (`h`, `hu`, and `hv` in the shallow
water case) to be separated by a fixed amount (`field_stride`).  I
also changed the physics interface so that flux computation and wave
speed functions, rather than dealing with a single cell, dealt with
arrays of cells.  This allowed me to move the flux and wave speed
computations into C code.  With a judicious use of `restrict`
pointers, I started to see some compiler vectorization.
Unfortunately, the compiler (both Intel and Clang) still told me that
it could not figure out how to vectorize the wave speed computation,
which involves independent per-cell computations and a max reduction.
I assumed the culprit to be the max function (and fixed the issue, but
after the `devec` tag).

The `devec` branch still has many of the same issues with conservative
aliasing analysis by the compiler, but it functions and it uses a
struct-of-arrays style.  However, moving routines into code that the
compiler could understand how to vectorize was increasingly pushing me
to writing C-style routines, even in the cases where I used the C++
compiler.

## C99 version

In the next round, I converted the code completely to C99.  This
involved incrementally moving functionality out of C++ code and into
C99 code that was called by the C++ code, until there was almost no
C++ code left.  Then I converted the small amount of data management
that was left in the C++ side.

As I converted from C++ to C, I made minor adjustments to try to get
better vectorization.  In many cases, this succeeded.  This version of
the code takes 10-15 ms per step, about a factor of four improvement
over the baseline code.

## Vectorized

The next (and current) stage of the code got me down to 7-10 ms per
frame for the 200-by-200 case -- not quite an order of magnitude
improvement over the original, but close.  The performance improvement
is more dramatic for larger grids, as the current version uses memory
more efficiently than the original as well.  There are several
individual improvements that contributed to this performance gain:

### Limiter improvements

The MinMod limiter was a nuisance for a long time.  I did various
things to try to improve performance, such as modestly reducing the
amount of arithmetic per cell by rewriting the minmod computation to
be a scaled version of minmod.  But the greatest improvement came from
replacing a call to `fminf` with a conditional statement (I used the
ternary operator in C, but I suspect that using an `if` statement
would have worked as well).

In addition to changing the limiter code, I also changed the code so
that the limited derivatives were computed just before they were used.
The derivatives of $F$ and $G$ only appear in the half-step predictor;
the derivatives of $U$ only appear in the corrector.  Moreover, if we
sweep from one side of the board to the other, we only need to be
concerned with the limited derivatives in two successive columns.
Computing these values a column at a time helps with vectorization;
computing only one column, just before it is needed, lets us save
memory and reduce cache pressure.

### Sums and differences

In the corrector step, we compute differences of $x$ and $y$
derivatives in two directions based on values at four points.
For the $x$ derivatives, we combine quantities `zkl` at distances
`k` and `l` from a base point through the stencil

    +z01 -z11
    +z00 -z10

We think of computing a sum for each row, i.e.

    s1 = z01 - z11
    s0 = z00 - z10

Similarly, we combine quantities for the $y$ derivative of a $w$ quantity
through the stencil

    +w01 +w11
    -w00 -w10

and we can think of adding up the rows (with a sign change) here as
well:

    d1 = w01 + w11
    d0 = w00 + w01

In this case, we have that the total contribution looks like `s0 +
s1 - (d0-d1)`; and when we shift up to the next row, there is no need
to re-compute the row sums.  This cuts down somewhat on the arithmetic
of the update.  It also simplifies the process of tracking which
limited derivatives we need at what time: when we compute the `s` and
`d` sums for a row, we incorporate any limited derivatives of `u` on
that row, and then don't have to touch those derivatives again.

### Alternating storage

In the original code, we compute each time step by writing into a
temporary array `v`, and then copying the data back into the reference
array `u`.  An alternate strategy is to alternate the role of the two
arrays, first moving from `u` to `v` and then back.  This reduces the
cost of the tuned code by another 10 or 15 percent.
