#include "sPHENIXFieldMap.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

void findCenter()
{
    sPHENIXFieldMap fmap("fieldMapFineFullField.csv","fieldMapRoughFullField.csv");

    const int    NZ   = fmap.GetNZ();
    const double zMin = fmap.GetZMin();
    const double dz   = fmap.GetDZ();

    // --- Bz maximum on axis (ir=0) ---
    double BzMax = -1e9;
    int    izBzMax = -1;
    for (int iz = 0; iz < NZ; ++iz) {
        double v = fmap.GetBzGrid(0, iz);
        if (v > BzMax) { BzMax = v; izBzMax = iz; }
    }
    double zBzMax = zMin + izBzMax * dz;
    printf("\nBz maximum on axis: Bz=%.5f T at iz=%d  z=%.1f mm\n",
           BzMax, izBzMax, zBzMax);

    // --- Br sign change at r=100 mm (ir=4) ---
    int irTest = 4;   // r = 100 mm
    printf("\nBr zero crossings at r=100 mm:\n");
    for (int iz = 1; iz < NZ; ++iz) {
        double b0 = fmap.GetBrGrid(irTest, iz-1);
        double b1 = fmap.GetBrGrid(irTest, iz  );
        if (b0 * b1 <= 0.0 && (b0 != 0.0 || b1 != 0.0)) {
            double z0 = zMin + (iz-1) * dz;
            double z1 = zMin +  iz    * dz;
            double zBrZero = z0 - b0 * (z1 - z0) / (b1 - b0);
            printf("  iz=%d..%d  z=%.1f..%.1f mm  Br=%.5f..%.5f T  → zero at z=%.2f mm\n",
                   iz-1, iz, z0, z1, b0, b1, zBrZero);
        }
    }

    // --- dBr/dz sign changes at r=100 mm ---
    printf("\ndBr/dz sign changes at r=100 mm:\n");
    for (int iz = 2; iz < NZ-1; ++iz) {
        double dm = (fmap.GetBrGrid(irTest,iz  ) - fmap.GetBrGrid(irTest,iz-2)) / (2*dz);
        double dp = (fmap.GetBrGrid(irTest,iz+1) - fmap.GetBrGrid(irTest,iz-1)) / (2*dz);
        if (dm * dp < 0.0) {
            double z = zMin + iz * dz;
            printf("  iz=%d  z=%.1f mm  (%.4e -> %.4e T/mm)\n", iz, z, dm, dp);
        }
    }

    // --- Bz and Br vs z in the ±160 mm window around the Bz peak ---
    printf("\n  iz   z[mm]   Bz_axis[T]   Br(r=100mm)[T]\n");
    int i0 = std::max(0, izBzMax - 8);
    int i1 = std::min(NZ-1, izBzMax + 8);
    for (int iz = i0; iz <= i1; ++iz) {
        double z = zMin + iz * dz;
        printf("  %3d  %7.1f  %9.5f    %+10.6f\n",
               iz, z, fmap.GetBzGrid(0,iz), fmap.GetBrGrid(irTest,iz));
    }
}
