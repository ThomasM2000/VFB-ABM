/*
 * parameters.h
 *
 * Calibration output helpers.
 *
 * Author: Kimberley Trickey
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>

using namespace std;
namespace util {

/*
 * Model outputs as defined in Manuscript Section 2.1, written one line per
 * call for the calibration / sensitivity driver.
 *
 * Columns are chosen to match what Section 2.2 actually measures:
 *   cells          <- PicoGreen dsDNA          (2.2.6)
 *   collagen       <- Sircol soluble collagen  (2.2.5)
 *   total protein  <- Bradford                 (2.2.4)
 *   viability      <- LIVE/DEAD                (2.2.7)
 * at days 0, 3, 6 and 9 (2.2).
 */
 void outputTotalChem(BMWorld* myWorld, string filename) {
  #ifdef CALIBRATION
    cout << "Outputting model outputs to: " << filename << endl;
    ofstream output_file;
    output_file.open(filename.c_str(), ios::app);
  
    const int totalPatches = myWorld->nx * myWorld->ny * myWorld->nz;
  
    /* ECM protein amounts - Section 2.1, in ug. */
    double col = 0, eln = 0, ha = 0, fha = 0;
    for (int in = 0; in < totalPatches; in++) {
      col += myWorld->worldECM[in].ncollagen[read_t];
      eln += myWorld->worldECM[in].nelastin[read_t];
      ha  += myWorld->worldECM[in].HA[read_t];
      fha += myWorld->worldECM[in].fHA[read_t];
    }
  
    /* Cytokine concentrations - Section 2.1, in pg. */
    const double tnf  = myWorld->world_total_tnf();
    const double tgf  = myWorld->world_total_tgf();
    const double fgf  = myWorld->world_total_fgf();
    const double il6  = myWorld->world_total_il6();
    const double il8  = myWorld->world_total_il8();
    const double il10 = myWorld->world_total_il10();
  
    /* Bradford total protein (2.2.4): proteinaceous species only, in ug.
     * HA is a glycosaminoglycan, not a protein, and is excluded. */
    const double PG_PER_UG = 1.0e6;
    const double totalProtein =
        col + eln + (tnf + tgf + fgf + il6 + il8 + il10) / PG_PER_UG;
  
    output_file << fixed << setprecision(5);
    output_file << myWorld->reportDay()        << "\t";  // t, days
    output_file << myWorld->cells.actualSize() << "\t";  // PicoGreen
    output_file << col  << "\t" << eln << "\t" << ha << "\t" << fha << "\t";
    output_file << totalProtein << "\t";                 // Bradford
    output_file << tnf << "\t" << tgf << "\t" << fgf << "\t"
                << il6 << "\t" << il8 << "\t" << il10 << "\t";
    /* Mechanical properties - Section 2.1. */
    output_file << BMWorld::E << "\t" << BMWorld::Q << "\t"
                << BMWorld::massLoss << "\t";
    output_file << endl;
  #endif  // ifdef CALIBRATION
    return;
  }
  
}  // namespace util

#endif /* PARAMETERS_H_ */
