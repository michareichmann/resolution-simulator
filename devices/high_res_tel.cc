// Simon Spannagel (DESY) January 2016

#include "assembly.h"
#include "propagate.h"
#include "materials.h"
#include "constants.h"
#include "log.h"

using namespace std;
using namespace gblsim;
using namespace unilog;

int main(int argc, char* argv[]) {

  /*
   * Create points on initial trajectory, create trajectory from points,
   * fit and write trajectory to MP-II binary file,
   * get track parameter corrections and covariance matrix at points.
   *
   * Equidistant measurement layers and thin scatterers, propagation
   * with simple jacobian (quadratic in arc length differences).
   */

  
  Log::ReportingLevel() = Log::FromString("INFO");

  for (int i = 1; i < argc; i++) {
    // Setting verbosity:
    if (std::string(argv[i]) == "-v") { 
      Log::ReportingLevel() = Log::FromString(std::string(argv[++i]));
      continue;
    } 
  }
  
  //----------------------------------------------------------------------------
  // Preparation of the particle trajectory:
  
  // Telescope properties:
  double analog_plane = 285e-3 / X0_Si + 500e-3 / X0_Si + 700e-3 / X0_PCB;
  double diamond_plane = 1550e-3 / X0_PCB + 700e-3 / X0_Si + 500e-3 / X0_Diamond; // + 40e-3 / X0_Au + 40e-3 / X0_Au + 10e-3 / X0_Au;
  double digital_plane = 1550e-3 / X0_PCB + 700e-3 / X0_Si + 285e-3 / X0_Si;

  // Beam: 260 MeV Pi at PSI
  double BEAM = 0.260;
  double spacing = 20.32;
  
  //----------------------------------------------------------------------------
  // Build the trajectory through the telescope device:

  plane pl0(0,analog_plane,true,resolution_analog);
  plane pl1(20.32,analog_plane,true,resolution_analog);

  plane diamond1(60.96,diamond_plane,false);
//  plane diamond2(81.28,diamond_plane,false);
  plane silicon(101.6,digital_plane,false);
  //plane silicon(101.6,digital_plane,false);
    
  plane pl2(142.24,analog_plane,true,resolution_analog);
  plane pl3(162.56,analog_plane,true,resolution_analog);

  std::vector<plane> planes;
  planes.push_back(pl0);
  planes.push_back(pl1);
  planes.push_back(diamond1);
//  planes.push_back(diamond2);
  planes.push_back(silicon);
  planes.push_back(pl2);
  planes.push_back(pl3);

  telescope mytel(planes, BEAM);
  LOG(logRESULT) << "Track resolution (X) at Diamond 1: " << mytel.getResolution(2);
  LOG(logRESULT) << "Track resolution (X) at Silicon 1: " << mytel.getResolution(3);


  plane ypl0(0,analog_plane,true,resolution_analog_y);
  plane ypl1(spacing,analog_plane,true,resolution_analog_y);

  plane dia1(spacing * 3, diamond_plane, false);
  plane sil1(spacing * 5, digital_plane, false);

  plane ypl3(spacing * 7, analog_plane,true,resolution_analog_y);
  plane ypl4(spacing * 8, analog_plane,true,resolution_analog_y);

  std::vector<plane> yplanes;
  yplanes.push_back(ypl0);
  yplanes.push_back(ypl1);
  yplanes.push_back(dia1);
  yplanes.push_back(sil1);
  yplanes.push_back(ypl3);
  yplanes.push_back(ypl4);

  telescope ymytel(yplanes, BEAM);
  LOG(logRESULT) << "Track resolution (Y) at Diamond: " << ymytel.getResolution(2);
  LOG(logRESULT) << "Track resolution (Y) at Silicon: " << ymytel.getResolution(3);

  return 0;
}
