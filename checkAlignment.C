// checkAlignment.C
// Checks alignment of the sPHENIX solenoid axis with the z axis.
//
// Strategy: for a perfectly aligned solenoid the transverse field components
// (Bx_phx, By_phx) average to zero over phi at every (r,z).  A rigid tilt
// by angle theta about some transverse axis appears as a systematic non-zero
// mean:   <B_transverse> = B_total * sin(theta) ~ B_total * theta
//
// We bin the raw fine-map measurements in z and accumulate the phi-averaged
// transverse field.  The tilt angle and its azimuthal direction (in the
// sPHENIX x-y plane) are then:
//
//   theta_x = atan2(<Bx>, <Bz>)   (tilt that rotates z toward x)
//   theta_y = atan2(<By>, <Bz>)   (tilt that rotates z toward y)
//
// Usage:
//   root -l -b -q 'checkAlignment.C'

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

void checkAlignment(const char* fineCSV = "fieldMapFineFullField.csv")
{
    // z binning – use same 20 mm grid as the field map
    const double zMin = -2700.0, zMax = 2100.0, dz = 20.0;
    const int NZ = int((zMax - zMin) / dz) + 1;

    // Accumulators per z bin: sum of Bx, By, Bz, and count
    std::vector<double> sumBx(NZ,0), sumBy(NZ,0), sumBz(NZ,0);
    std::vector<int>    npts(NZ,0);

    // Also keep inner-arm (r < 150 mm) separately
    std::vector<double> sumBx_in(NZ,0), sumBy_in(NZ,0), sumBz_in(NZ,0);
    std::vector<int>    npts_in(NZ,0);

    std::ifstream f(fineCSV);
    if (!f.is_open()) {
        printf("Cannot open %s\n", fineCSV); return;
    }
    printf("Reading %s …\n", fineCSV);

    long long nread = 0;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        double xs, ys, zs, Bmag, Bxs, Bys, Bzs;
        if (std::sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                        &xs,&ys,&zs,&Bmag,&Bxs,&Bys,&Bzs) != 7) continue;
        ++nread;

        // Coordinate transform to sPHENIX
        const double xp = xs, yp = zs, zp = -ys;
        const double r  = std::sqrt(xp*xp + yp*yp);

        // Field components in sPHENIX Cartesian
        const double Bx_p =  Bxs;   // Bx_phx = Bx_s
        const double By_p =  Bzs;   // By_phx = Bz_s
        const double Bz_p = -Bys;   // Bz_phx = -By_s

        int iz = int((zp - zMin) / dz);
        if (iz < 0 || iz >= NZ) continue;

        sumBx[iz] += Bx_p;
        sumBy[iz] += By_p;
        sumBz[iz] += Bz_p;
        npts[iz]  += 1;

        if (r < 150.0) {
            sumBx_in[iz] += Bx_p;
            sumBy_in[iz] += By_p;
            sumBz_in[iz] += Bz_p;
            npts_in[iz]  += 1;
        }
    }
    printf("  read %lld rows\n", nread);

    // ── Build histograms ────────────────────────────────────────────────────
    TH1D* hBx  = new TH1D("hBx", ";z [mm];#LT B_{x} #GT [T]", NZ, zMin, zMax);
    TH1D* hBy  = new TH1D("hBy", ";z [mm];#LT B_{y} #GT [T]", NZ, zMin, zMax);
    TH1D* hTx  = new TH1D("hTx", ";z [mm];#theta_{x} [mrad]", NZ, zMin, zMax);
    TH1D* hTy  = new TH1D("hTy", ";z [mm];#theta_{y} [mrad]", NZ, zMin, zMax);

    // Inner-arm versions
    TH1D* hBx_in = new TH1D("hBx_in","",NZ,zMin,zMax);
    TH1D* hBy_in = new TH1D("hBy_in","",NZ,zMin,zMax);
    TH1D* hTx_in = new TH1D("hTx_in","",NZ,zMin,zMax);
    TH1D* hTy_in = new TH1D("hTy_in","",NZ,zMin,zMax);

    // Print table header
    printf("\n  iz   z[mm]   n   <Bx>[mT]   <By>[mT]   theta_x[mrad]  theta_y[mrad]\n");

    for (int iz = 0; iz < NZ; ++iz) {
        if (npts[iz] == 0) continue;
        double z   = zMin + iz * dz;
        double bx  = sumBx[iz] / npts[iz];
        double by  = sumBy[iz] / npts[iz];
        double bz  = sumBz[iz] / npts[iz];
        double tx  = (bz != 0) ? 1e3 * std::atan2(bx, bz) : 0;
        double ty  = (bz != 0) ? 1e3 * std::atan2(by, bz) : 0;

        hBx->SetBinContent(iz+1, bx * 1e3);  // mT
        hBy->SetBinContent(iz+1, by * 1e3);
        hTx->SetBinContent(iz+1, tx);
        hTy->SetBinContent(iz+1, ty);

        if (std::abs(z) < 1500)   // print only the central region
            printf("  %3d  %6.0f  %4d  %+8.3f   %+8.3f   %+8.4f       %+8.4f\n",
                   iz, z, npts[iz], bx*1e3, by*1e3, tx, ty);

        if (npts_in[iz] > 0) {
            double bx2 = sumBx_in[iz]/npts_in[iz];
            double by2 = sumBy_in[iz]/npts_in[iz];
            double bz2 = sumBz_in[iz]/npts_in[iz];
            hBx_in->SetBinContent(iz+1, bx2*1e3);
            hBy_in->SetBinContent(iz+1, by2*1e3);
            hTx_in->SetBinContent(iz+1, (bz2!=0)?1e3*std::atan2(bx2,bz2):0);
            hTy_in->SetBinContent(iz+1, (bz2!=0)?1e3*std::atan2(by2,bz2):0);
        }
    }

    // ── Global phi-averaged tilt (all z, all r) ────────────────────────────
    double gBx=0,gBy=0,gBz=0; long long gN=0;
    for (int iz=0;iz<NZ;++iz){ if(npts[iz]==0) continue;
        gBx+=sumBx[iz]; gBy+=sumBy[iz]; gBz+=sumBz[iz]; gN+=npts[iz]; }
    if (gN>0) {
        gBx/=gN; gBy/=gN; gBz/=gN;
        printf("\nGlobal phi-averaged transverse field (all z, all r, N=%lld):\n", gN);
        printf("  <Bx> = %+.4f mT    <By> = %+.4f mT    <Bz> = %.4f T\n",
               gBx*1e3, gBy*1e3, gBz);
        printf("  theta_x = %+.4f mrad    theta_y = %+.4f mrad\n",
               1e3*std::atan2(gBx,gBz), 1e3*std::atan2(gBy,gBz));
        double theta = 1e3*std::sqrt(gBx*gBx+gBy*gBy)/std::abs(gBz);
        double phi_t = std::atan2(gBy,gBx)*180.0/M_PI;
        printf("  |theta| = %.4f mrad    tilt direction phi = %.1f deg\n", theta, phi_t);
    }

    // ── Draw ───────────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1","Solenoid alignment",1200,800);
    c1->Divide(2,2);

    auto style = [](TH1D* h, int col) {
        h->SetLineColor(col); h->SetLineWidth(2);
        h->SetMarkerColor(col); h->SetMarkerStyle(20); h->SetMarkerSize(0.4);
    };

    c1->cd(1);
    style(hBx, kBlue+1); style(hBx_in, kRed+1);
    hBx->SetTitle("#LT B_{x} #GT vs z (all r = blue, r<150 mm = red);z [mm];#LT B_{x} #GT [mT]");
    hBx->Draw("HIST"); hBx_in->Draw("HIST SAME");
    TLine* l0a = new TLine(zMin,0,zMax,0); l0a->SetLineStyle(2); l0a->Draw();

    c1->cd(2);
    style(hBy, kBlue+1); style(hBy_in, kRed+1);
    hBy->SetTitle("#LT B_{y} #GT vs z (all r = blue, r<150 mm = red);z [mm];#LT B_{y} #GT [mT]");
    hBy->Draw("HIST"); hBy_in->Draw("HIST SAME");
    TLine* l0b = new TLine(zMin,0,zMax,0); l0b->SetLineStyle(2); l0b->Draw();

    c1->cd(3);
    style(hTx, kBlue+1); style(hTx_in, kRed+1);
    hTx->SetTitle("#theta_{x} = atan2(#LT B_{x} #GT, #LT B_{z} #GT);z [mm];#theta_{x} [mrad]");
    hTx->Draw("HIST"); hTx_in->Draw("HIST SAME");
    TLine* l0c = new TLine(zMin,0,zMax,0); l0c->SetLineStyle(2); l0c->Draw();

    c1->cd(4);
    style(hTy, kBlue+1); style(hTy_in, kRed+1);
    hTy->SetTitle("#theta_{y} = atan2(#LT B_{y} #GT, #LT B_{z} #GT);z [mm];#theta_{y} [mrad]");
    hTy->Draw("HIST"); hTy_in->Draw("HIST SAME");
    TLine* l0d = new TLine(zMin,0,zMax,0); l0d->SetLineStyle(2); l0d->Draw();

    c1->SaveAs("fieldMap_alignment.pdf");
    printf("\nSaved fieldMap_alignment.pdf\n");
}
