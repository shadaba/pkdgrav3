To: Douglas Potter <douglas.potter@uzh.ch>, Joachim Stadel <stadel@physik.uzh.ch>
Subject: pkdgrav3: wrong sign of the Ewald potential for the leftover SIMD images (one-line fix)

Dear Doug and Joachim,

We found a bug in the Ewald correction in pkdgrav3. It is present in 3.5.1 and in the current master of
bitbucket.org/dpotter/pkdgrav3 (gravity/ewald.cxx, line 252, master dcce745f14). It affects only the
**potential**; accelerations are correct.

## The bug
`pkdParticleEwald()` (gravity/ewald.cxx, `USE_SIMD_EWALD`) evaluates the real-space Ewald images in batches of
SIMD width. `evalEwaldSIMD()` -> `evalEwald()` *subtracts* each image's potential from its accumulator
(`fPot -= g0*m - ...`).

- Full batches are accumulated straight into `dPot`, which is correct.
- The last, partial batch (0 to width-1 images) is evaluated into a zeroed `tpot`, masked, and then
  `dPot -= tpot;`. That subtracts a second time, so the potential of those images is added with the wrong sign.
- `dax += tax` (and y, z) is correct, which is why forces are unaffected.

Fix (patch attached, `ewald_simd_remainder_sign.patch`):

```diff
-        dPot -= tpot;
+        dPot += tpot;
```

How many images land in the last partial batch depends on the particle's position (how many images lie within
`dEwCut`), on `nReplicas` and on the SIMD width: 2 on ARM/NEON, 4 with AVX, 8 with AVX-512. So the error
looks like position-dependent offsets. The CUDA Ewald path has no such remainder and looks fine.

## Symptoms we saw
Settings: production-like runs, `bEwald = 1`, `dEwCut = 2.6`, `classic_theta_switch()` and
`classic_replicas_switch()` (so `nReplicas = 2` above z = 20). Potentials are shown in peculiar (km/s)^2
(phi_code * U^2 / a).

1. **x86-64 with AVX-512 (64 threads), 128^3, L = 12.5 Mpc/h, z = 49 -> 13 in 10 steps**, run with the binary
   before and after the fix, all parameters the same:
   - Positions and velocities are **bit-identical** at z = 19.3 and z = 13.
   - The potential differs for 20.6-20.7% of particles.
   - Before the fix, 5.1% of particles have potentials 10 to 1000 times deeper than the median (minimum
     -5.3e6 (km/s)^2 at z = 13). After the fix the minimum is -9.3e4 and the median is -7036 (km/s)^2.
   - The same median, -7036, comes from an independent run on an Apple M4 (NEON) with the same parameters.
2. **512^3, L = 50 Mpc/h, 64 threads, x86-64 with AVX-512, before the fix:**
   - About 5% of particles have potentials of -5e6 to -8.5e7 (km/s)^2, mixed with normal particles in the same
     2 Mpc/h cells.
   - Restarting from a checkpoint changed which particles are affected.
3. **Apple M4 (NEON, width 2), 64^3, L = 6.25 Mpc/h, one step to z = 45, 8 threads:**
   - Before the fix, `nReplicas = 2` (or 3) minus `nReplicas = 1` is -0.17 in code units (-3.7e5 (km/s)^2) for
     interior particles and 0 for particles within an eighth of the box of a face.
   - The difference does not depend on dTheta (0.4 or 0.7), so it is not the tree walk.
   - After the fix, the difference is -7 +- 11 (km/s)^2 for 2 replicas and -9 +- 14 for 3. The median with 1
     replica is -13381 (km/s)^2, which is the softened self-interaction -(105/64) m/eps.

## Consequences
Potentials written by unfixed builds (the HDF5 `Potential` column, `bDoPotOutput` `.pot` files, and
`pkd->ga[].minPot` used for group centres) are wrong for some particles. The energy diagnostic `dEnergyU` is
affected too. Dynamics are not affected.

Best regards,
Shadab Alam

(Found while adding a compressed snapshot writer to pkdgrav3 3.5.1: https://github.com/shadaba/pkdgrav3,
branch main, commit 2219e678. The tests are bench/phase5/ewfix.pbs and section 19 of
docs/compressed_io_plan.md there.)
