#pragma once

#include <string>
#include <array>

// Magnetic field map for the sPHENIX solenoid.
//
// Input CSVs are in the surveyor coordinate system:
//   surveyor (x,y,z) = sPHENIX (x, -z, y)
//   i.e.  x_phx = x_s,  y_phx = z_s,  z_phx = -y_s
//
// This class works internally in sPHENIX cylindrical coordinates:
//   r   = sqrt(x_phx^2 + y_phx^2)  [mm]
//   phi = atan2(y_phx, x_phx)       [rad]
//   z   = z_phx = -y_s              [mm]  (beam axis)
//
// The "best" field is built by:
//   1. Averaging all azimuthal measurements into a regular (r,z) grid.
//      Fine-map data takes precedence; rough-map fills the extended z tails.
//   2. Enforcing ∇·B = 0 numerically:
//      Br(r,z) = -(1/r) ∫₀ʳ r′ ∂Bz/∂z dr′
//      This replaces the measured Br with a physically consistent one.
//      For a current-free region ∇×B = 0 is also satisfied to the degree
//      that the averaged Bz(r,z) is accurate.
//   3. Bφ is set to zero (azimuthal symmetry of the solenoid).
//
// All lengths in mm, field components in Tesla.

class sPHENIXFieldMap {
public:
    // Load and build the map.  Both CSV files are required.
    sPHENIXFieldMap(const std::string& fineCSV,
                    const std::string& roughCSV);

    // Default constructor (required for ROOT dictionary; produces an empty map).
    sPHENIXFieldMap();

    // Field in sPHENIX cylindrical coordinates.
    // Points outside the grid boundary return zero.
    void GetField(double r, double phi, double z,
                  double& Br, double& Bphi, double& Bz) const;

    // Field in sPHENIX Cartesian coordinates.
    void GetFieldXYZ(double x, double y, double z,
                     double& Bx, double& By, double& Bz) const;

    // Coverage of the combined map after building
    double GetRMin()  const { return kRMin; }
    double GetRMax()  const { return kRMax; }
    double GetZMin()  const { return kZMin; }
    double GetZMax()  const { return kZMax; }
    double GetDR()    const { return kdR;   }
    double GetDZ()    const { return kdZ;   }
    int    GetNR()    const { return kNR;   }
    int    GetNZ()    const { return kNZ;   }

    // Read-only access to the internal grids (for diagnostics / ROOT histograms)
    double GetBzGrid (int ir, int iz) const { return fBz  [ir][iz]; }
    double GetBrGrid (int ir, int iz) const { return fBr  [ir][iz]; }
    int    GetNGrid  (int ir, int iz) const { return fNpts[ir][iz]; }

private:
    // ── Grid layout ──────────────────────────────────────────────────────────
    // r vertices:  0, 25, 50, … , 900 mm        (37 nodes, 25 mm step)
    // z vertices: -2700, -2680, … , 2100 mm     (241 nodes, 20 mm step)
    // Node (0,*) is the axis; Br is forced to zero there.
    static constexpr int    kNR   = 37;
    static constexpr int    kNZ   = 241;
    static constexpr double kRMin =    0.0;   // mm
    static constexpr double kRMax =  900.0;   // mm
    static constexpr double kZMin = -2700.0;  // mm
    static constexpr double kZMax =  2100.0;  // mm
    static constexpr double kdR   =   25.0;   // mm
    static constexpr double kdZ   =   20.0;   // mm

    // Grid arrays [r-index][z-index]
    double fBz  [kNR][kNZ];   // averaged measured Bz
    double fBzWt[kNR][kNZ];   // accumulated weight (for averaging)
    double fBr  [kNR][kNZ];   // Maxwell-consistent Br (computed)
    int    fNpts[kNR][kNZ];   // number of contributing measurements

    // ── Internal helpers ─────────────────────────────────────────────────────
    // Load one CSV file.  If overwriteExisting=false the rough map only fills
    // cells that the fine map left empty (fNpts == 0).
    void LoadCSV(const std::string& filename, bool overwriteExisting);

    // Replace fBr with the Maxwell-consistent values derived from fBz.
    void EnforceMaxwell();

    // Bilinear interpolation of an (kNR × kNZ) grid at (r, z).
    // Returns 0 if (r,z) is outside the grid.
    double Interp(const double grid[kNR][kNZ],
                  double r, double z) const;

    // Nearest grid indices (clamped)
    static int IROf(double r) {
        int ir = static_cast<int>((r - kRMin) / kdR);
        if (ir < 0) ir = 0;
        if (ir >= kNR - 1) ir = kNR - 2;
        return ir;
    }
    static int IZOf(double z) {
        int iz = static_cast<int>((z - kZMin) / kdZ);
        if (iz < 0) iz = 0;
        if (iz >= kNZ - 1) iz = kNZ - 2;
        return iz;
    }
};
