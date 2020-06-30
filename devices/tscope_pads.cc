// Simon Spannagel (DESY) January 2016

#include "assembly.h"
#include "propagate.h"
#include "materials.h"
#include "constants.h"
#include "log.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"


using namespace std;
using namespace gblsim;
using namespace unilog;

/** ----------------------------------------------------------------------------
 ** Telescope properties */
double m_cms_pixel = 285e-3 / X0_Si + 500e-3 / X0_Si + 700e-3 / X0_PCB;  // without cover
//double m_cms_pixel = 285e-3 / X0_Si + 500e-3 / X0_Si + 700e-3 / X0_PCB + 30 / X0_PVC + 10 / X0_Kapton;
double m_diamond_pad = 20e-3 / X0_Al + 500e-3 / X0_Diamond + 20e-3 / X0_Al;
double BEAM = 0.260;  // Beam: 260 MeV Pi at PSI
/** ---------------------------------------------------------------------------- */

vector<plane> build_telescope(char mode, unsigned n_planes, float spacing, float dut_distance){
  vector<plane> tmp;
  for (unsigned i(0); i < n_planes; i++) {
    float dut_spacing = i >= n_planes / 2 ? dut_distance - spacing : 0;
    tmp.emplace_back(spacing * i + dut_spacing, m_cms_pixel, true, mode == 'x' ? r_cms_pix_x : r_cms_pix_y); }
  return tmp;
}

void draw_plane_spacing(unsigned n_planes, float dut_distance, float dut_spacing, float step=.5){

  vector<float> res_x;
  vector<float> res_y;
  vector<float> spacings;
  for (float s(1); s < 60;){
    for (auto mode: {'x', 'y'}){
      auto tel = build_telescope(mode, n_planes, s, dut_spacing);
      float tel_center = (n_planes / 2. - 1) * s + dut_spacing / 2;
      tel.insert(tel.begin() + n_planes / 2, plane(tel_center + dut_distance / 2, m_diamond_pad, false));
      tel.insert(tel.begin() + n_planes / 2, plane(tel_center - dut_distance / 2, m_diamond_pad, false));
      telescope mytel(tel, BEAM);
      if (mode == 'x') {
        res_x.emplace_back(mytel.getResolution(n_planes / 2)); }
      else {
        res_y.emplace_back(mytel.getResolution(n_planes / 2)); }
    }
    spacings.emplace_back(s);
    s += step;
  }

  auto c = new TCanvas("c","bla", 1600, 800);
  c->Divide(2);
  vector<TGraph> graphs;
  for (int i(0); i < 2; i++){
    auto pad = c->cd(i + 1);
    auto g = new TGraph(spacings.size(), spacings.data(), i ? res_y.data() : res_x.data());
    g->SetNameTitle("g", TString::Format("Resolution in %s", i ? "Y" : "X"));
    g->SetMarkerSize(.5);
    g->SetMarkerStyle(20);
    g->SetLineColor(2);
    g->Draw("al");
    auto x_ax = g->GetXaxis();
    x_ax->SetTitle("Z_{t} [mm]");
    auto y_ax = g->GetYaxis();
    y_ax->SetTitle("Resolution [#mum]");
    y_ax->SetTitleOffset(1.6);
    pad->SetLeftMargin(.12);
  }
  c->SaveAs("res_spacing.pdf");
}

void draw_dut_gap(unsigned n_planes, float dut_distance, float spacing, float step=.5){

  vector<float> res_x;
  vector<float> res_y;
  vector<float> distances;
  for (float d(30); d < 200;){
    for (auto mode: {'x', 'y'}){
      auto tel = build_telescope(mode, n_planes, spacing, d);
      float tel_center = (n_planes / 2. - 1) * spacing + d / 2;
      tel.insert(tel.begin() + n_planes / 2, plane(tel_center + dut_distance / 2, m_diamond_pad, false));
      tel.insert(tel.begin() + n_planes / 2, plane(tel_center - dut_distance / 2, m_diamond_pad, false));
      telescope mytel(tel, BEAM);
      if (mode == 'x') {
        res_x.emplace_back(mytel.getResolution(n_planes / 2)); }
      else {
        res_y.emplace_back(mytel.getResolution(n_planes / 2)); }
    }
    distances.emplace_back(d);
    d += step;
  }

  auto c = new TCanvas("c0","bla", 1600, 800);
  c->Divide(2);
  vector<TGraph> graphs;
  for (int i(0); i < 2; i++){
    auto pad = c->cd(i + 1);
    auto g = new TGraph(distances.size(), distances.data(), i ? res_y.data() : res_x.data());
    g->SetNameTitle("g", TString::Format("Resolution in %s", i ? "Y" : "X"));
    g->SetMarkerSize(.5);
    g->SetMarkerStyle(20);
    g->SetLineColor(2);
    g->Draw("al");
    auto x_ax = g->GetXaxis();
    x_ax->SetTitle("Z_{d} [mm]");
    auto y_ax = g->GetYaxis();
    y_ax->SetTitle("Resolution [#mum]");
    y_ax->SetTitleOffset(1.6);
    pad->SetLeftMargin(.12);
  }
  c->SaveAs("res_gap.pdf");
}

int main(int argc, char* argv[]) {

  /** ----------------------------------------------------------------------------
   ** VERBOSITY */
  Log::ReportingLevel() = Log::FromString("INFO");
  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-v") {
      Log::ReportingLevel() = Log::FromString(std::string(argv[++i])); } }

  
  /** ----------------------------------------------------------------------------
   ** Build the trajectory through the telescope device: */
  LOG(logRESULT) << "Four-plane tracking:";
  LOG(logINFO) << "Tracking Plane thickness [%X0]: " << setprecision(3) << m_cms_pixel * 100;
  LOG(logINFO) << "DUT thickness [%X0]: " << setprecision(3) << m_diamond_pad * 100;

  draw_plane_spacing(4, 15, 5 * 20);
  draw_dut_gap(4, 15, 30);

  return 0;
}
