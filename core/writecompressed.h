/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Compressed parallel snapshot output (snapcomp, docs/snapcomp_format.md).
 *  Every leaf thread compresses its own particles and appends blocks to one file per subsample that is
 *  shared by all threads of its MPI process. No gathering, no reordering.
 */
#ifndef SERVICE_WRITECOMPRESSED_H
#define SERVICE_WRITECOMPRESSED_H
#include <algorithm>
#include <cstdint>
#include "TraversePST.h"

struct WriteCompressedStats {
    std::uint64_t nWritten = 0;    // particles written
    std::uint64_t nSkipped = 0;    // non-dark particles (not written)
    std::uint64_t nBytes = 0;      // compressed block bytes
    double dEncode = 0;            // thread-seconds in the encoder (sum)
    double dWrite = 0;             // thread-seconds in appends (sum)
    double dMaxThread = 0;         // slowest thread, seconds
    std::uint64_t nScratch = 0;    // writer scratch bytes, summed over threads
    std::uint64_t nScratchMax = 0; // largest per-thread scratch
    WriteCompressedStats& operator+=(const WriteCompressedStats& r) {
        nWritten += r.nWritten;
        nSkipped += r.nSkipped;
        nBytes += r.nBytes;
        dEncode += r.dEncode;
        dWrite += r.dWrite;
        dMaxThread = std::max(dMaxThread, r.dMaxThread);
        nScratch += r.nScratch;
        nScratchMax = std::max(nScratchMax, r.nScratchMax);
        return *this;
    }
};

class ServiceWriteCompressed : public TraverseCount<WriteCompressedStats> {
  public:
    static constexpr int MAX_SUB = 16;
    struct input {
        char achPrefix[1024];          // snapshot name; files are <prefix>.s<pct>.<process>.pkc
        char achInfo[1024];            // free text stored in every header
        double dCenter[3], dPeriod[3]; // periodic box in code units
        double dvFac;                  // code velocity -> dx/dt (as for the other writers)
        double dVelToKms;              // dx/dt -> peculiar km/s
        double dBoxMpch, dRedshift, dTime, dOmega0, dLambda, dHubble;
        double dVelStep; // km/s
        double dPercents[MAX_SUB];
        std::uint64_t nTotal, nLattice, iStep, iSeed;
        std::uint32_t nCell, nSub, nBlock, nPercents;
        std::int32_t iLevel;
    };
    static constexpr std::string_view name = "WriteCompressed";
    explicit ServiceWriteCompressed(PST pst, std::string_view service_name = name) : TraverseCount(pst, service_name, sizeof(input)) {}

  protected:
    virtual int Service(PST pst, void* vin, int nIn, void* vout, int nOut);
};
#endif
