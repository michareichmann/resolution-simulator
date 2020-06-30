// Simon Spannagel (DESY) January 2016

#include "assembly.h"
#include "propagate.h"
#include "materials.h"
#include "constants.h"
#include "log.h"
#include "math.h"

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
  float angle_fac = sqrt(1 + pow(tan(30 * M_PI / 180), 2) + pow(tan(35 * M_PI / 180), 2));
  double diamond_plane = (1550e-3 / X0_PCB + 700e-3 / X0_Si + 500e-3 / X0_Diamond) * angle_fac; // + 40e-3 / X0_Au + 40e-3 / X0_Au + 10e-3 / X0_Au;
  double digital_plane = (1550e-3 / X0_PCB + 700e-3 / X0_Si + 285e-3 / X0_Si) * angle_fac;

  // Beam: 200 GeV Pi at CERN
  double BEAM = 200;
//  double spacing = 20.32;
  double spacing = 60;

  //----------------------------------------------------------------------------
  // Build the trajectory through the telescope device:

  std::vector<plane> planes;
  for (int i(0); i < 3; i++)
    planes.push_back(plane(spacing * i, digital_plane, true, resolution_analog));

  for (int i(6); i < 9; i++)
    planes.push_back(plane(spacing * i, diamond_plane, false));

  for (int i(12); i < 15; i++)
    planes.push_back(plane(spacing * i, digital_plane, true, resolution_analog));

  telescope mytel(planes, BEAM);

  LOG(logRESULT) << "Track resolution (X) at Diamond 1: " << mytel.getResolution(3);
  LOG(logRESULT) << "Track resolution (X) at Diamond 2: " << mytel.getResolution(4);
  LOG(logRESULT) << "Track resolution (X) at Diamond 3: " << mytel.getResolution(5);


  std::vector<plane> yplanes;
  for (int i(0); i < 3; i++)
    yplanes.push_back(plane(spacing * i, digital_plane, true, resolution_analog_y));

  for (int i(3); i < 6; i++)
    yplanes.push_back(plane(spacing * i, diamond_plane, false));

  for (int i(6); i < 9; i++)
    yplanes.push_back(plane(spacing * i, digital_plane, true, resolution_analog_y));

  telescope ymytel(yplanes, BEAM);
  LOG(logRESULT) << "Track resolution (Y) at Diamond 1: " << ymytel.getResolution(3);
  LOG(logRESULT) << "Track resolution (Y) at Diamond 2: " << ymytel.getResolution(4);
  LOG(logRESULT) << "Track resolution (Y) at Diamond 3: " << ymytel.getResolution(5);

  return 0;
}
