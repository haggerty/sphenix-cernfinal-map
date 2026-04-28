#include "sPHENIXFieldMap.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>

// ============================================================
// Default constructor (empty map, for ROOT dictionary)
// ============================================================
sPHENIXFieldMap::sPHENIXFieldMap()
{
    std::memset(fBz,   0, sizeof(fBz));
    std::memset(fBzWt, 0, sizeof(fBzWt));
    std::memset(fBr,   0, sizeof(fBr));
    std::memset(fNpts, 0, sizeof(fNpts));
}

// ============================================================
// Constructor
// ============================================================
sPHENIXFieldMap::sPHENIXFieldMap(const std::string& fineCSV,
                                 const std::string& roughCSV)
{
    // Zero-initialise all grids
    std::memset(fBz,   0, sizeof(fBz));
    std::memset(fBzWt, 0, sizeof(fBzWt));
    std::memset(fBr,   0, sizeof(fBr));
    std::memset(fNpts, 0, sizeof(fNpts));

    // Fine map fills the grid first; rough fills only empty cells
    LoadCSV(fineCSV,  /*overwriteExisting=*/true);
    LoadCSV(roughCSV, /*overwriteExisting=*/false);

    // Normalise accumulated sums → means
    for (int ir = 0; ir < kNR; ++ir)
        for (int iz = 0; iz < kNZ; ++iz)
            if (fBzWt[ir][iz] > 0.0)
                fBz[ir][iz] /= fBzWt[ir][iz];

    // Derive Maxwell-consistent Br from the averaged Bz field
    EnforceMaxwell();

    // Report coverage
    int filled = 0;
    for (int ir = 0; ir < kNR; ++ir)
        for (int iz = 0; iz < kNZ; ++iz)
            if (fNpts[ir][iz] > 0) ++filled;

    std::printf("sPHENIXFieldMap: %d / %d grid nodes filled  "
                "(r=[0,900] mm, z=[-2700,2100] mm)\n",
                filled, kNR * kNZ);
}

// ============================================================
// LoadCSV
// Columns: x_s  y_s  z_s  |B|  Bx_s  By_s  Bz_s   (mm, T)
//
// Coordinate transforms to sPHENIX:
//   x_phx = x_s          y_phx = z_s         z_phx = -y_s
//   Bx_phx = Bx_s        By_phx = Bz_s       Bz_phx = -By_s
//
// Cylindrical (around z_phx = beam axis):
//   r   = sqrt(x_phx^2 + y_phx^2)
//   phi = atan2(y_phx, x_phx)
//   Bz_cyl = Bz_phx
//   Br     = Bx_phx*cos(phi) + By_phx*sin(phi)   (stored only for diagnostics)
// ============================================================
void sPHENIXFieldMap::LoadCSV(const std::string& filename,
                               bool overwriteExisting)
{
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Cannot open: " + filename);

    std::printf("Loading %s …\n", filename.c_str());

    long long nread = 0, nused = 0;
    std::string line;

    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;

        double xs, ys, zs, Bmag, Bxs, Bys, Bzs;
        if (std::sscanf(line.c_str(),
                        "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                        &xs, &ys, &zs, &Bmag, &Bxs, &Bys, &Bzs) != 7)
            continue;
        ++nread;

        // ── Coordinate transform ──────────────────────────────────────────
        const double xp  =  xs;
        const double yp  =  zs;
        const double zp  = -ys;

        const double r   = std::sqrt(xp*xp + yp*yp);
        const double phi = std::atan2(yp, xp);

        // Only Bz enters the Maxwell enforcement; we accumulate it.
        // Bz_phx = -By_s
        const double Bz_phx = -Bys;

        // ── Grid placement ────────────────────────────────────────────────
        // Skip points outside the grid domain
        if (r  < kRMin || r  > kRMax + kdR) continue;
        if (zp < kZMin || zp > kZMax + kdZ) continue;

        // Nearest node (bilinear deposit: split weight among 4 neighbours)
        const double rFrac = (r  - kRMin) / kdR;
        const double zFrac = (zp - kZMin) / kdZ;

        const int ir0 = static_cast<int>(rFrac);
        const int iz0 = static_cast<int>(zFrac);
        const int ir1 = ir0 + 1;
        const int iz1 = iz0 + 1;

        if (ir0 < 0 || ir1 >= kNR || iz0 < 0 || iz1 >= kNZ) continue;

        const double wr1 = rFrac - ir0;   // weight towards ir1
        const double wr0 = 1.0 - wr1;
        const double wz1 = zFrac - iz0;
        const double wz0 = 1.0 - wz1;

        // Four corner weights
        const double w[2][2] = {
            {wr0*wz0, wr0*wz1},
            {wr1*wz0, wr1*wz1}
        };
        const int IRs[2] = {ir0, ir1};
        const int IZs[2] = {iz0, iz1};

        for (int a = 0; a < 2; ++a) {
            for (int b = 0; b < 2; ++b) {
                const int ir = IRs[a];
                const int iz = IZs[b];
                const double wab = w[a][b];
                if (wab == 0.0) continue;

                if (!overwriteExisting && fNpts[ir][iz] > 0) continue;

                fBz  [ir][iz] += wab * Bz_phx;
                fBzWt[ir][iz] += wab;
                fNpts[ir][iz]++;   // coarse count; precise normalisation uses fBzWt
            }
        }

        ++nused;
        (void)phi;   // phi is not needed for the cylindrical-symmetric average
    }

    std::printf("  read %lld rows, used %lld\n", nread, nused);
}

// ============================================================
// EnforceMaxwell
//
// From ∇·B = 0 in cylindrical coordinates (azimuthal symmetry):
//
//   (1/r) ∂(r Br)/∂r + ∂Bz/∂z = 0
//
//   → r·Br(r,z) = -∫₀ʳ r′ (∂Bz/∂z)(r′,z) dr′
//
//   → Br(r,z) = -(1/r) ∫₀ʳ r′ (∂Bz/∂z)(r′,z) dr′
//
// Discretisation (trapezoidal rule, vertex-centred grid):
//
//   Br[ir][iz] = -(dr / r[ir]) * Σ_{k=0}^{ir} w_k * r[k] * (∂Bz/∂z)[k][iz]
//
// where w_0 = w_ir = 0.5, w_k = 1 (k interior), r[0] = 0 → w_0 · r[0] = 0.
//
// ∂Bz/∂z is computed by central differences; forward/backward at the edges.
//
// Cells with no measurements are filled by nearest-available-r extrapolation
// before computing the z-derivative, so the integral is always contiguous.
// ============================================================
void sPHENIXFieldMap::EnforceMaxwell()
{
    // ── Step 1a: linear z-interpolation BETWEEN data points only ────────
    // For each r-shell with any measurements, fill the interior gaps by
    // linear interpolation.  Leading/trailing zeros are left as-is so that
    // step 1b can choose a better r-neighbour for those z positions instead
    // of propagating a stale fine-map edge value into the rough-map-only
    // z region.
    for (int ir = 0; ir < kNR; ++ir) {
        int first = -1, last_iz = -1;
        for (int iz = 0; iz < kNZ; ++iz) {
            if (fNpts[ir][iz] > 0) {
                if (first < 0) first = iz;
                last_iz = iz;
            }
        }
        if (first < 0) continue;   // no data for this r-shell

        // Linear interpolation between consecutive data points
        int prev = first;
        for (int iz = first + 1; iz <= last_iz; ++iz) {
            if (fNpts[ir][iz] > 0) {
                double b0 = fBz[ir][prev], b1 = fBz[ir][iz];
                for (int jz = prev + 1; jz < iz; ++jz) {
                    double t = double(jz - prev) / double(iz - prev);
                    fBz[ir][jz] = (1.0-t)*b0 + t*b1;
                }
                prev = iz;
            }
        }
        // fBz[ir][iz] = 0 for iz < first and iz > last_iz (left for step 1b)
    }

    // ── Step 1b: fill shells in the r direction, per z slice ─────────────
    // Outward scan then inward scan propagate the nearest non-zero value to
    // any zero cell.  Key effect: in the z range where only the rough map
    // (r = 100, 200, … mm) has data, small-r shells (r = 50, 75 mm) pick
    // up the rough-map Bz rather than a stale fine-map constant extrapolation.
    for (int iz = 0; iz < kNZ; ++iz) {
        // Outward: fills gaps between arm positions on the outer side
        double last = 0.0;
        for (int ir = 0; ir < kNR; ++ir) {
            if (fBz[ir][iz] != 0.0) last = fBz[ir][iz];
            else if (last != 0.0)   fBz[ir][iz] = last;
        }
        // Inward: fills ir < innermost arm (e.g. ir = 0, 1)
        last = 0.0;
        for (int ir = kNR - 1; ir >= 0; --ir) {
            if (fBz[ir][iz] != 0.0) last = fBz[ir][iz];
            else if (last != 0.0)   fBz[ir][iz] = last;
        }
    }

    // ── Step 1c: z-fill for cells still zero after steps 1a + 1b ─────────
    // These are z positions beyond the full map coverage (|z| > 2457 mm).
    // Propagate the nearest edge value flat into the extrapolated tails.
    for (int ir = 0; ir < kNR; ++ir) {
        double last = 0.0;
        for (int iz = 0; iz < kNZ; ++iz) {
            if (fBz[ir][iz] != 0.0) last = fBz[ir][iz];
            else if (last != 0.0)   fBz[ir][iz] = last;
        }
        last = 0.0;
        for (int iz = kNZ - 1; iz >= 0; --iz) {
            if (fBz[ir][iz] != 0.0) last = fBz[ir][iz];
            else if (last != 0.0)   fBz[ir][iz] = last;
        }
    }

    // ── Step 1d: Gaussian smoothing in z (σ ≈ 50 mm = 2.5 cells) ────────
    // Applied after all fill steps so that the fine/rough boundary kink and
    // any staircase left by the 100-mm rough-map z spacing are smoothed out
    // before computing ∂Bz/∂z.  The field varies on scales >> 50 mm so
    // physics is not affected.
    {
        const int    half  = 6;        // truncate at ±6 cells
        const double sigma = 2.5;      // cells  (= 50 mm)
        const double inv2s2 = 1.0 / (2.0 * sigma * sigma);

        // Pre-compute kernel weights
        double kern[2*half+1];
        for (int d = -half; d <= half; ++d)
            kern[d + half] = std::exp(-d*d * inv2s2);

        for (int ir = 0; ir < kNR; ++ir) {
            double tmp[kNZ];
            for (int iz = 0; iz < kNZ; ++iz) {
                double sum = 0.0, wsum = 0.0;
                for (int d = -half; d <= half; ++d) {
                    int jz = iz + d;
                    if (jz < 0 || jz >= kNZ) continue;
                    double w = kern[d + half];
                    sum  += w * fBz[ir][jz];
                    wsum += w;
                }
                tmp[iz] = sum / wsum;
            }
            for (int iz = 0; iz < kNZ; ++iz)
                fBz[ir][iz] = tmp[iz];
        }
    }

    // ── Step 2: compute ∂Bz/∂z at every (ir, iz) ─────────────────────────
    double dBzdz[kNR][kNZ];
    const double inv2dz = 1.0 / (2.0 * kdZ);

    for (int ir = 0; ir < kNR; ++ir) {
        // Central differences in the interior
        for (int iz = 1; iz < kNZ - 1; ++iz)
            dBzdz[ir][iz] = (fBz[ir][iz+1] - fBz[ir][iz-1]) * inv2dz;

        // One-sided differences at the boundaries
        dBzdz[ir][0]       = (fBz[ir][1]       - fBz[ir][0]      ) / kdZ;
        dBzdz[ir][kNZ - 1] = (fBz[ir][kNZ - 1] - fBz[ir][kNZ - 2]) / kdZ;
    }

    // ── Step 3: trapezoidal integration from r=0 outward ─────────────────
    // Build a running sum S[iz] = Σ_{k=0}^{ir} w_k * r_k * dBzdz[k][iz]
    // in a single O(kNR·kNZ) sweep instead of O(kNR²·kNZ).
    //
    // Br[0][iz] = 0  (boundary: B_r = 0 on axis)
    // Br[ir][iz] = -(dr / r_ir) * S_ir[iz]
    //
    // Trapezoidal weights: w_0 = 0.5, w_interior = 1, w_ir = 0.5
    // Since r_0 = 0, the k=0 term always vanishes.

    for (int iz = 0; iz < kNZ; ++iz)
        fBr[0][iz] = 0.0;

    // S[iz] accumulates the trapezoidal integral up to the previous ir.
    // Initialise to the half-weight contribution of k=0: r_0=0 → zero.
    double S[kNZ] = {};

    for (int ir = 1; ir < kNR; ++ir) {
        const double r_ir   = kRMin + ir * kdR;
        const double r_prev = kRMin + (ir-1) * kdR;

        for (int iz = 0; iz < kNZ; ++iz) {
            // Add the interior contribution of the previous node (full weight)
            // and the half-weight endpoint of the current node.
            // The step from ir-1 to ir uses the trapezoidal rule:
            //   ΔS = 0.5*(r_{ir-1}*dBzdz[ir-1] + r_ir*dBzdz[ir]) * dr
            // We keep S as the cumulative sum of full-weight interior terms
            // and add the endpoint half-weight when computing Br.
            S[iz] += r_prev * dBzdz[ir-1][iz];   // interior weight for ir-1

            // Br[ir] = -(dr/r_ir) * [S + 0.5*r_ir*dBzdz[ir]]
            fBr[ir][iz] = -(kdR / r_ir) *
                           (S[iz] + 0.5 * r_ir * dBzdz[ir][iz]);
        }
    }
}

// ============================================================
// Interp – bilinear interpolation on a (kNR × kNZ) grid
// Returns 0 outside [kRMin,kRMax] × [kZMin,kZMax].
// ============================================================
double sPHENIXFieldMap::Interp(const double grid[kNR][kNZ],
                                double r, double z) const
{
    if (r < kRMin || r > kRMax || z < kZMin || z > kZMax) return 0.0;

    const double rFrac = (r - kRMin) / kdR;
    const double zFrac = (z - kZMin) / kdZ;

    int ir = static_cast<int>(rFrac);
    int iz = static_cast<int>(zFrac);

    // Clamp so that ir+1 and iz+1 are valid
    if (ir >= kNR - 1) ir = kNR - 2;
    if (iz >= kNZ - 1) iz = kNZ - 2;

    const double fr = rFrac - ir;
    const double fz = zFrac - iz;

    return (1.0-fr)*(1.0-fz)*grid[ir  ][iz  ]
         + (    fr)*(1.0-fz)*grid[ir+1][iz  ]
         + (1.0-fr)*(    fz)*grid[ir  ][iz+1]
         + (    fr)*(    fz)*grid[ir+1][iz+1];
}

// ============================================================
// GetField – public interface, cylindrical coordinates
// ============================================================
void sPHENIXFieldMap::GetField(double r, double /*phi*/, double z,
                                double& Br, double& Bphi, double& Bz) const
{
    Bphi = 0.0;
    Bz   = Interp(fBz, r, z);
    Br   = Interp(fBr, r, z);
}

// ============================================================
// GetFieldXYZ – public interface, sPHENIX Cartesian coordinates
// ============================================================
void sPHENIXFieldMap::GetFieldXYZ(double x, double y, double z,
                                   double& Bx, double& By, double& Bz) const
{
    const double r   = std::sqrt(x*x + y*y);
    const double phi = std::atan2(y, x);

    double Br, Bphi, Bzcyl;
    GetField(r, phi, z, Br, Bphi, Bzcyl);

    if (r > 0.0) {
        const double cos_phi = x / r;
        const double sin_phi = y / r;
        Bx = Br * cos_phi - Bphi * sin_phi;
        By = Br * sin_phi + Bphi * cos_phi;
    } else {
        Bx = 0.0;
        By = 0.0;
    }
    Bz = Bzcyl;
}
