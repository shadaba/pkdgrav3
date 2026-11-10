/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MSR::WriteCompressed: snapshot output in the snapcomp format (bCompressedOutput). See
 *  core/writecompressed.h, io/snapcomp/ and docs/snapcomp_format.md.
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include "master.h"

#include <cmath>
#include <ctime>
#include <unistd.h>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/writecompressed.h"
#include "snapcomp.h"

using namespace fmt::literals;

void MSR::WriteCompressed(std::string_view FileName, double dTime) {
    using S = ServiceWriteCompressed;
    S::input in;
    std::memset(&in, 0, sizeof(in));
    if (FileName.size() >= sizeof(in.achPrefix)) throw std::runtime_error("compressed output: file name too long");
    std::memcpy(in.achPrefix, FileName.data(), FileName.size());
    if (!csm->val.bComove) throw std::runtime_error("compressed output needs a cosmological (comoving) run");

    const double a = csmTime2Exp(csm, dTime);
    const double L = units.dKpcUnit * 1e-3 * csm->val.h; // Mpc/h
    in.dvFac = 1.0 / (a * a);                            // code velocity -> dx/dt, as in AllNodeWrite
    in.dVelToKms = a * L * 100.0 / csm->val.dHubble0;    // dx/dt (box/code time) -> peculiar km/s
    auto period = parameters.get_dPeriod();
    for (int d = 0; d < 3; ++d) {
        in.dCenter[d] = fCenter[d];
        in.dPeriod[d] = period[d];
    }
    in.dBoxMpch = L;
    in.dRedshift = 1.0 / a - 1.0;
    in.dTime = dTime;
    in.dOmega0 = csm->val.dOmega0;
    in.dLambda = csm->val.dLambda;
    in.dHubble = csm->val.h;
    in.dVelStep = parameters.get_dCompVelStep();
    // Potential: comoving code units (G = 1, box = 1); peculiar physical potential = phi * U^2 / a with U the
    // code velocity unit in km/s (checked on Mac snapshots: <v^2>/well depth is then independent of a).
    in.dPotStep = parameters.get_dCompPotStep();
    in.dPotToKms2 = std::pow(L * 100.0 / csm->val.dHubble0, 2) / a;
    if (in.dPotStep < 0) throw std::runtime_error("compressed output: dCompPotStep must be >= 0");
    // Only some snapshots (low redshift) with the potential: dCompPotMaxRedshift.
    const bool bPotHere = in.dPotStep > 0 && 1.0 / a - 1.0 <= parameters.get_dCompPotMaxRedshift() + 1e-9;
    if (!bPotHere) in.dPotStep = 0;
    in.nTotal = N;
    in.iSeed = uint64_t(parameters.get_iCompSubsampleSeed());
    in.nBlock = uint32_t(parameters.get_iCompBlockSize());
    in.iLevel = int32_t(parameters.get_iCompLevel());
    in.iWriteThrough = int32_t(parameters.get_iCompWriteThrough());
    in.dBytesPerProcess = 10.0 * double(N) / double(mdl->Procs()); // generous: measured 5-9 B/particle
    auto n = uint64_t(std::llround(std::cbrt(double(N))));
    in.nLattice = (n * n * n == N) ? n : 0; // IDs from pkdgrav3 ICs are (iz*N+iy)*N+ix; blocks fall back if not

    std::vector<double> pc = parameters.get_lstCompSubsamples();
    if (pc.empty()) pc = {1, 3, 5, 9, 10, 12, 15, 20, 25}; // Sahyadri
    if (pc.size() > size_t(S::MAX_SUB)) throw std::runtime_error("compressed output: at most 16 subsamples");
    double sum = 0;
    for (size_t k = 0; k < pc.size(); ++k) sum += in.dPercents[k] = pc[k];
    if (std::fabs(sum - 100.0) > 1e-9) throw std::runtime_error("compressed output: lstCompSubsamples must sum to 100");
    in.nPercents = uint32_t(pc.size());
    if (parameters.get_dCompPosStep() <= 0 || in.dVelStep <= 0) throw std::runtime_error("compressed output: steps must be positive");
    if (in.nBlock < 1) throw std::runtime_error("compressed output: iCompBlockSize must be >= 1");
    snapcomp::choose_grid(L * 1000.0, parameters.get_dCompPosStep(), N, in.nCell, in.nSub);

    // Self-description in the info text (JSON; the binary header keeps only Omega_m, Omega_lambda and h).
    // Groups follow the simulation part of an AbacusSummit header: sim (provenance), box, ic, cosmology (with
    // derived quantities at this snapshot, from pkdgrav3's own csm, so w0-wa/radiation/CLASS runs are covered),
    // numerics, steps.
    auto js = [](std::string s) {
        std::string o;
        for (char c : s) {
            if (c == '"' || c == '\\') o += '\\';
            o += c;
        }
        return o;
    };
    const auto& cv = csm->val;
    const double H0 = cv.dHubble0, h = cv.h, z = 1.0 / a - 1.0;
    // step number from the snapshot name (<name>.<step>), as BuildName makes it
    {
        auto dot = FileName.rfind('.');
        auto s = (dot == std::string_view::npos) ? std::string_view() : FileName.substr(dot + 1);
        if (!s.empty() && s.find_first_not_of("0123456789") == std::string_view::npos) in.iStep = std::stoull(std::string(s));
    }
    // sim
    std::string simName(parameters.get_achSimName());
    if (simName.empty()) {
        char cwd[4096];
        if (getcwd(cwd, sizeof(cwd))) {
            std::string c(cwd);
            simName = c.substr(c.find_last_of('/') + 1);
        }
    }
    char host[256] = "";
    gethostname(host, sizeof(host) - 1);
    char when[32] = "";
    std::time_t now = std::time(nullptr);
    std::strftime(when, sizeof(when), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
#ifndef PKD_GIT_DESCRIBE
    #define PKD_GIT_DESCRIBE ""
#endif
    std::string sim = fmt::format("{{\"name\":\"{}\",\"set\":\"{}\",\"comment\":\"{}\",\"code\":\"pkdgrav3 {}\",\"code_git\":\"{}\",\"written_utc\":\"{}\","
                                  "\"host\":\"{}\",\"mpi_processes\":{},\"threads\":{},\"gpu\":{}}}",
                                  js(simName), js(std::string(parameters.get_achSimSet())), js(std::string(parameters.get_achSimComment())), PACKAGE_VERSION,
                                  PKD_GIT_DESCRIBE, when, js(host), mdl->Procs(), mdl->Threads(),
#ifdef USE_CUDA
                                  parameters.get_bGPU() ? "true" : "false"
#else
                                  "false"
#endif
    );
    // box
    const double mpart = 2.77536627e11 * cv.dOmega0 * L * L * L / double(N); // Msun/h
    const bool bSoftGlobal = parameters.has_dSoft();
    const double epsCode = bSoftGlobal ? parameters.get_dSoft() : (parameters.get_nGrid() > 0 ? 1.0 / (50.0 * parameters.get_nGrid()) : 0.0);
    std::string box = fmt::format("{{\"L_mpch\":{:.10g},\"L_mpc\":{:.10g},\"n_particles\":{},\"n_grid\":{},\"particle_mass_msunh\":{:.10g},"
                                  "\"particle_mass_msun\":{:.10g},\"softening_kpch\":{:.6g},\"softening_source\":\"{}\",\"softening_max\":{:.6g},"
                                  "\"softening_max_is_multiplier\":{},\"max_physical_softening\":{:.6g}}}",
                                  L, L / h, N, parameters.get_nGrid(), mpart, mpart / h, epsCode * L * 1e3,
                                  bSoftGlobal ? "dSoft" : (parameters.get_nGrid() > 0 ? "default 1/(50 nGrid)" : "per particle (input file)"),
                                  parameters.get_dSoftMax(), parameters.get_bSoftMaxMul() ? "true" : "false", parameters.get_dMaxPhysicalSoft());
    // ic
    std::string ic = fmt::format("{{\"z_initial\":{:.10g},\"lpt_order\":{},\"seed\":{},\"fixed_amplitude\":{},\"phase_pi\":{:.6g},\"input_file\":\"{}\","
                                 "\"tf_file\":\"{}\",\"class\":{},\"class_file\":\"{}\"}}",
                                 parameters.get_dRedFrom(), parameters.get_iLPT(), parameters.get_iSeed(), parameters.get_bFixedAmpIC() ? "true" : "false",
                                 parameters.get_dFixedAmpPhasePI(), js(std::string(parameters.get_achInFile())), js(std::string(parameters.get_achTfFile())),
                                 cv.classData.bClass ? "true" : "false", js(std::string(parameters.get_achClassFilename())));
    // cosmology: parameters and derived quantities at this snapshot
    auto E = [&](double x) { return csmExp2Hub(csm, x) / H0; };
    const double Ez = E(a);
    const double Ok = 1.0 - cv.dOmega0 - cv.dLambda - cv.dOmegaDE - cv.dOmegaRad;
    const double rhoDE = cv.dLambda + cv.dOmegaDE * std::pow(a, -3.0 * (1.0 + cv.w0 + cv.wa)) * std::exp(-3.0 * cv.wa * (1.0 - a));
    double D1, D2, D3a, D3b, D3c, f1, f2, f3a, f3b, f3c, D10, f10;
    csmComoveGrowth(csm, a, &D1, &D2, &D3a, &D3b, &D3c, &f1, &f2, &f3a, &f3b, &f3c);
    csmComoveGrowth(csm, 1.0, &D10, &D2, &D3a, &D3b, &D3c, &f10, &f2, &f3a, &f3b, &f3c);
    double chi = 0; // comoving distance, Mpc/h: (c/100) int_a^1 da / (a^2 E), Simpson in a
    {
        const int n = 2000;
        const double da = (1.0 - a) / n;
        for (int i = 0; i <= n; ++i) {
            double x = a + i * da, w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
            chi += w / (x * x * E(x));
        }
        chi *= da / 3.0 * 2997.92458;
    }
    const double tH = 9.777922216807891 / h; // 1/H0 in Gyr
    std::string cosmo =
        fmt::format("{{\"h\":{:.10g},\"H0\":{:.10g},\"Omega_m\":{:.10g},\"Omega_lambda\":{:.10g},\"Omega_b\":{:.10g},\"Omega_rad\":{:.10g},"
                    "\"Omega_DE\":{:.10g},\"Omega_K\":{:.10g},\"w0\":{:.10g},\"wa\":{:.10g},\"omega_b\":{:.10g},\"omega_cdm\":{:.10g},\"sigma8\":{:.10g},"
                    "\"normalization\":{:.10g},\"n_s\":{:.10g},\"running\":{:.10g},\"pivot\":{:.10g},\"class\":{},\"class_file\":\"{}\",\"tf_file\":\"{}\","
                    "\"redshift\":{:.10g},\"scale_factor\":{:.10g},\"E_z\":{:.10g},\"H_z_kms_mpc\":{:.10g},\"Omega_m_z\":{:.10g},\"Omega_DE_z\":{:.10g},"
                    "\"growth_D\":{:.10g},\"growth_f\":{:.10g},\"age_gyr\":{:.10g},\"hubble_time_gyr\":{:.10g},\"comoving_distance_mpch\":{:.10g}}}",
                    h, 100 * h, cv.dOmega0, cv.dLambda, cv.dOmegab, cv.dOmegaRad, cv.dOmegaDE, Ok, cv.w0, cv.wa, cv.dOmegab * h * h,
                    (cv.dOmega0 - cv.dOmegab) * h * h, cv.dSigma8, cv.dNormalization, cv.dSpectral, cv.dRunning, cv.dPivot,
                    cv.classData.bClass ? "true" : "false", js(std::string(parameters.get_achClassFilename())), js(std::string(parameters.get_achTfFile())), z,
                    a, Ez, 100 * h * Ez, cv.dOmega0 / (a * a * a * Ez * Ez), rhoDE / (Ez * Ez), D1 / D10, f1, csmExp2Time(csm, a) * H0 * tH, tH, chi);
    // numerics
    std::string num = fmt::format("{{\"theta\":{:.6g},\"replicas\":{},\"ewald\":{},\"ewald_cut\":{:.6g},\"eta\":{:.6g},\"eps_acc_step\":{},\"new_kdk\":{}}}",
                                  parameters.get_dTheta(), parameters.get_nReplicas(), parameters.get_bEwald() ? "true" : "false", parameters.get_dEwCut(),
                                  parameters.get_dEta(), parameters.get_bEpsAccStep() ? "true" : "false", parameters.get_bNewKDK() ? "true" : "false");
    // steps: output redshifts (dRedTo) and steps per segment
    auto nSteps = parameters.get_nSteps();
    auto zTo = parameters.get_dRedTo();
    int64_t nTot = 0;
    std::string ns, zs;
    bool same = true;
    for (size_t k = 0; k < nSteps.size(); ++k) {
        nTot += nSteps[k];
        same &= nSteps[k] == nSteps[0];
        ns += (k ? "," : "") + std::to_string(nSteps[k]);
    }
    for (size_t k = 0; k < zTo.size(); ++k) zs += fmt::format("{}{:.6g}", k ? "," : "", zTo[k]);
    std::string steps = fmt::format("{{\"step\":{},\"total_steps\":{},\"steps_per_segment\":{},\"segment_end_redshifts\":[{}],\"output_interval\":{}}}",
                                    in.iStep, nTot, same && !nSteps.empty() ? std::to_string(nSteps[0]) : "[" + ns + "]", zs, parameters.get_iOutInterval());
    std::string info =
        fmt::format("{{\"code\":\"pkdgrav3 {}\",\"pos_step_kpch_requested\":{},\"vel_step_kms\":{},\"pot_step_kms2\":{},\"name\":\"{}\","
                    "\"sim\":{},\"box\":{},\"ic\":{},\"cosmology\":{},\"numerics\":{},\"steps\":{}}}",
                    PACKAGE_VERSION, parameters.get_dCompPosStep(), in.dVelStep, in.dPotStep, js(std::string(FileName)), sim, box, ic, cosmo, num, steps);
    if (info.size() >= sizeof(in.achInfo)) throw std::runtime_error(fmt::format("compressed output: info text too long ({} bytes)", info.size()));
    std::strncpy(in.achInfo, info.c_str(), sizeof(in.achInfo) - 1);

    print_detail("Writing {name} in compressed format in parallel ({files} files per process, dx {dx:.4g} kpc/h, dv {dv:g} km/s{pot}) ...\n",
                 "name"_a = FileName, "files"_a = pc.size(), "dx"_a = L * 1e3 / (double(in.nCell) * in.nSub), "dv"_a = in.dVelStep,
                 "pot"_a = in.dPotStep > 0 ? fmt::format(", potential dphi {:g} (km/s)^2", in.dPotStep)
                           : parameters.get_dCompPotStep() > 0
                               ? fmt::format(", no potential (z > dCompPotMaxRedshift = {:g})", parameters.get_dCompPotMaxRedshift())
                               : std::string());
    print_detail("Time:{time:g} Redshift:{z:g}\n", "time"_a = dTime, "z"_a = in.dRedshift);

    TimerStart(TIMER_IO);
    S::output out;
    mdl->RunService<S>(sizeof(in), &in, &out, sizeof(out));
    TimerStop(TIMER_IO);
    double dsec = TimerGet(TIMER_IO);

    if (out.nWritten + out.nSkipped != N) throw std::runtime_error("compressed output: particle count mismatch");
    if (out.nSkipped) print_detail("WARNING: {n} non-dark-matter particles were not written (compressed output is dark matter only)\n", "n"_a = out.nSkipped);
    print_detail("Compressed {n} particles into {mb:.1f} MB ({bpp:.3f} B/particle); encode {enc:.2f} s, append {app:.2f} s summed over threads, slowest "
                 "thread {mx:.2f} s\n",
                 "n"_a = out.nWritten, "mb"_a = out.nBytes / 1e6, "bpp"_a = out.nWritten ? double(out.nBytes) / out.nWritten : 0.0, "enc"_a = out.dEncode,
                 "app"_a = out.dWrite, "mx"_a = out.dMaxThread);
    if (in.dPotStep > 0)
        print_detail("  of which potential: {mb:.1f} MB ({bpp:.3f} B/particle)\n", "mb"_a = out.nPotBytes / 1e6,
                     "bpp"_a = out.nWritten ? double(out.nPotBytes) / out.nWritten : 0.0);
    print_detail("Compressed output scratch memory: {tot:.1f} MB over all threads, {mx:.1f} MB per thread (freed after the write)\n",
                 "tot"_a = out.nScratch / 1e6, "mx"_a = out.nScratchMax / 1e6);
    print_detail("Compressed output write-through (page cache bypassed): {on} of {n} threads; node MemAvailable {gb:.1f} GiB, iCompWriteThrough {m}\n",
                 "on"_a = out.nWriteThrough, "n"_a = out.nThreads, "gb"_a = out.dMemAvailGB < 1e29 ? out.dMemAvailGB : -1.0,
                 "m"_a = parameters.get_iCompWriteThrough());
    print_detail("Output file has been successfully written, Wallclock: {seconds:.5} secs.\n", "seconds"_a = dsec);
}
