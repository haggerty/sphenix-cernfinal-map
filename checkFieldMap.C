// checkFieldMap.C
// ROOT macro: build the sPHENIX field map and make diagnostic plots.
//
// Usage (from the directory containing the CSV and source files):
//   root -l -b -q 'sPHENIXFieldMap.cxx+' 'checkFieldMap.C'
//
// The first argument compiles sPHENIXFieldMap into a shared library
// before the macro runs.
//
// Plots produced (saved as PDF):
//   1. Bz(r,z)  – averaged measured axial field
//   2. Br(r,z)  – Maxwell-consistent radial field (from ∇·B = 0)
//   3. ∇·B(r,z) – should be ≈ 0 everywhere
//   4. Bz on axis vs z
//   5. Bz and Br at z = 0 vs r

#include "sPHENIXFieldMap.h"

#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <cmath>
#include <cstdio>
#include <algorithm>

void checkFieldMap(const char* dir = ".")
{
    sPHENIXFieldMap fmap(
        Form("%s/fieldMapFineFullField.csv",  dir),
        Form("%s/fieldMapRoughFullField.csv", dir));

    const int    NR   = fmap.GetNR();
    const int    NZ   = fmap.GetNZ();
    const double dr   = fmap.GetDR();
    const double dz   = fmap.GetDZ();
    const double rMin = fmap.GetRMin();
    const double rMax = fmap.GetRMax();
    const double zMin = fmap.GetZMin();
    const double zMax = fmap.GetZMax();

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    // ── 2-D histograms ────────────────────────────────────────────────────
    TH2D* hBz   = new TH2D("hBz",
        "B_{z} [T] (#phi-averaged, sPHENIX coords);z [mm];r [mm]",
        NZ, zMin, zMax, NR, rMin, rMax);
    TH2D* hBr   = new TH2D("hBr",
        "B_{r} [T] (Maxwell-consistent, from #nabla#cdotB = 0);z [mm];r [mm]",
        NZ, zMin, zMax, NR, rMin, rMax);
    TH2D* hDivB = new TH2D("hDivB",
        "#nabla#cdotB [T/mm];z [mm];r [mm]",
        NZ, zMin, zMax, NR, rMin, rMax);

    for (int ir = 1; ir < NR - 1; ++ir) {
        double r = rMin + ir * dr;
        for (int iz = 1; iz < NZ - 1; ++iz) {
            double z = zMin + iz * dz;

            hBz->Fill(z, r, fmap.GetBzGrid(ir, iz));
            hBr->Fill(z, r, fmap.GetBrGrid(ir, iz));

            // ∇·B = (1/r) d(r Br)/dr  +  dBz/dz
            double dBzdz   = (fmap.GetBzGrid(ir, iz+1) -
                              fmap.GetBzGrid(ir, iz-1)) / (2.0*dz);
            double rBr_p   = (r + dr) * fmap.GetBrGrid(ir+1, iz);
            double rBr_m   = (r - dr) * fmap.GetBrGrid(ir-1, iz);
            double d_rBr_dr = (rBr_p - rBr_m) / (2.0*dr);

            hDivB->Fill(z, r, d_rBr_dr/r + dBzdz);
        }
    }

    // ── 1-D slices ────────────────────────────────────────────────────────
    TH1D* hBzAxis = new TH1D("hBzAxis",
        "B_{z} on axis (r = 0);z [mm];B_{z} [T]",
        NZ, zMin, zMax);
    TH1D* hBzMid  = new TH1D("hBzMid",
        "B_{z} at z = 0;r [mm];B_{z} [T]",
        NR, rMin, rMax);
    TH1D* hBrMid  = new TH1D("hBrMid",
        "B_{r} at z = 0;r [mm];B_{r} [T]",
        NR, rMin, rMax);

    for (int iz = 0; iz < NZ; ++iz) {
        double z = zMin + iz * dz;
        double Br, Bphi, Bz;
        fmap.GetField(0.0, 0.0, z, Br, Bphi, Bz);
        hBzAxis->SetBinContent(iz + 1, Bz);
    }
    for (int ir = 0; ir < NR; ++ir) {
        double r = rMin + ir * dr;
        double Br, Bphi, Bz;
        fmap.GetField(r, 0.0, 0.0, Br, Bphi, Bz);
        hBzMid->SetBinContent(ir + 1, Bz);
        hBrMid->SetBinContent(ir + 1, Br);
    }

    // ── Draw ─────────────────────────────────────────────────────────────
    TCanvas* c1 = new TCanvas("c1", "sPHENIX Field Map", 1400, 900);
    c1->Divide(2, 2);

    c1->cd(1);
    hBz->Draw("COLZ");

    c1->cd(2);
    hBr->Draw("COLZ");

    c1->cd(3);
    hDivB->GetZaxis()->SetRangeUser(-5e-4, 5e-4);
    hDivB->Draw("COLZ");

    c1->cd(4);
    hBzAxis->SetLineColor(kBlue + 1);
    hBzAxis->SetLineWidth(2);
    hBzAxis->Draw("HIST");

    TCanvas* c2 = new TCanvas("c2", "Mid-plane slices", 800, 400);
    c2->Divide(2, 1);
    c2->cd(1);
    hBzMid->SetLineColor(kBlue + 1); hBzMid->SetLineWidth(2);
    hBzMid->Draw("HIST");
    c2->cd(2);
    hBrMid->SetLineColor(kRed + 1); hBrMid->SetLineWidth(2);
    hBrMid->Draw("HIST");

    c1->SaveAs("fieldMap_overview.pdf");
    c2->SaveAs("fieldMap_midplane.pdf");
    printf("Saved fieldMap_overview.pdf and fieldMap_midplane.pdf\n");

    // ── Maxwell check ─────────────────────────────────────────────────────
    // Compute mean and RMS of ∇·B by iterating over filled bins directly.
    double sum_d = 0, sum2_d = 0, max_d = 0;
    int    n_d = 0;
    for (int ix = 1; ix <= hDivB->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= hDivB->GetNbinsY(); ++iy) {
            double v = hDivB->GetBinContent(ix, iy);
            if (v == 0.0) continue;
            sum_d  += v;
            sum2_d += v*v;
            max_d   = std::max(max_d, std::abs(v));
            ++n_d;
        }
    double mean_d = (n_d > 0) ? sum_d / n_d : 0.0;
    double rms_d  = (n_d > 0) ? std::sqrt(sum2_d/n_d - mean_d*mean_d) : 0.0;

    printf("\n=== Maxwell check: ∇·B [T/mm] over %d grid cells ===\n", n_d);
    printf("  mean  = %+.3e\n", mean_d);
    printf("  RMS   =  %.3e\n", rms_d);
    printf("  |max| =  %.3e\n", max_d);
    printf("  (Bz ~ 1.4 T, map size ~ 4.8 m  →  scale ~ 3e-4 T/mm)\n");

    // ── Spot-check GetFieldXYZ ────────────────────────────────────────────
    printf("\n=== GetFieldXYZ spot check [mm → T] ===\n");
    struct Pt { double x, y, z; };
    Pt pts[] = {
        {  0,   0,    0},
        {200,   0,    0},
        {  0,   0,  500},
        {300,   0,  500},
        {600,   0, 1000},
    };
    for (auto& p : pts) {
        double Bx, By, Bz;
        fmap.GetFieldXYZ(p.x, p.y, p.z, Bx, By, Bz);
        printf("  (%5.0f,%5.0f,%5.0f)  Bx=%8.5f  By=%8.5f  Bz=%8.5f  |B|=%.5f T\n",
               p.x, p.y, p.z, Bx, By, Bz,
               std::sqrt(Bx*Bx + By*By + Bz*Bz));
    }
}
