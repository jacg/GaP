#include "nain4.hh"
#include "g4-mandatory.hh"
#include "GasProperties.hh"

#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4OpticalPhysics.hh>
#include <G4SystemOfUnits.hh>
#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>

#include <iostream>
#include <assert.h>
#include <vector>

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////

G4double XenonRefractiveIndex(G4double energy, G4double density)
{
  // Formula for the refractive index taken from
  // A. Baldini et al., "Liquid Xe scintillation calorimetry
  // and Xe optical properties", arXiv:physics/0401072v1 [physics.ins-det]

  // The Lorentz-Lorenz equation (also known as Clausius-Mossotti equation)
  // relates the refractive index of a fluid with its density:
  // (n^2 - 1) / (n^2 + 2) = - A · d_M,     (1)
  // where n is the refractive index, d_M is the molar density and
  // A is the first refractivity viral coefficient:
  // A(E) = \sum_i^3 P_i / (E^2 - E_i^2),   (2)
  // with:
  G4double P[3] = {71.23, 77.75, 1384.89}; // [eV^3 cm3 / mole]
  G4double E[3] = {8.4, 8.81, 13.2};       // [eV]

  // Note.- Equation (1) has, actually, a sign difference with respect
  // to the one appearing in the reference. Otherwise, it yields values
  // for the refractive index below 1.

  // Let's calculate the virial coefficient.
  // We won't use the implicit system of units of Geant4 because
  // it results in loss of numerical precision.

  energy = energy / eV;
  G4double virial = 0.;

  for (G4int i=0; i<3; i++)
  virial = virial + P[i] / (energy*energy - E[i]*E[i]);

  // Need to use g/cm3
  density = density / g * cm3;

  G4double mol_density = density / 131.29;
  G4double alpha = virial * mol_density;

  // Isolating now the n2 from equation (1) and taking the square root
  G4double n2 = (1. - 2*alpha) / (1. + alpha);
  if (n2 < 1.) {
    //      G4String msg = "Non-physical refractive index for energy "
    // + bhep::to_string(energy) + " eV. Use n=1 instead.";
    //      G4Exception("[XenonProperties]", "RefractiveIndex()",
    // 	  JustWarning, msg);
    n2 = 1.;
  }

  return sqrt(n2);
  
 }
 
 //////////////////////////////////////////////////////////////////////
 
 G4double GXeScintillation(G4double energy, G4double pressure)
{
  // FWHM and peak of emission extracted from paper:
  // Physical review A, Volume 9, Number 2,
  // February 1974. Koehler, Ferderber, Redhead and Ebert.
  // Pressure must be in atm = bar
  // XXX Check if there is some newest results.

  pressure = pressure / atmosphere;

  G4double Wavelength_peak  = (0.05 * pressure + 169.45) * nm;

  G4double Wavelength_sigma = 0.;
  if (pressure < 4.)
  Wavelength_sigma = 14.3 * nm;
  else
  Wavelength_sigma = (-0.117 * pressure + 15.16) * nm / (2.*sqrt(2*log(2)));

  G4double Energy_peak  = (h_Planck * c_light / Wavelength_peak);
  G4double Energy_sigma = (h_Planck * c_light * Wavelength_sigma / pow(Wavelength_peak,2));

  G4double intensity = exp(-pow(Energy_peak/eV-energy/eV,2) /
  (2*pow(Energy_sigma/eV, 2))) /
  (Energy_sigma/eV*sqrt(pi*2.));

  return intensity;
}

/////////////////////////////////////////////////////////////////////////

G4double GXeDensity(G4double pressure)
{
  // Computes Xe (gas) density at T = 293 K
  // Values are taken from the reference file nexus/data/gxe_density_table.txt
  // (which, in turn, is downloaded from https://webbook.nist.gov/chemistry/fluid).
  // We assume a linear interpolation between any pair of values in the database.

  G4double density;

  const G4int n_pressures = 6;
  G4double data[n_pressures][2] = {{  1.0 * bar,   5.419 * kg/m3},
                                   {  5.0 * bar,  27.721 * kg/m3},
                                   { 10.0 * bar,  57.160 * kg/m3},
                                   { 13.5 * bar,  78.949 * kg/m3},
                                   { 20.0 * bar, 122.510 * kg/m3},
                                   { 30.0 * bar, 199.920 * kg/m3}};
  G4bool found = false;

  for (G4int i=0; i<n_pressures-1; ++i) {
    if  (pressure >= data[i][0] && pressure < data[i+1][0]) {
      G4double x1 = data[i][0];
      G4double x2 = data[i+1][0];
      G4double y1 = data[i][1];
      G4double y2 = data[i+1][1];
      density = y1 + (y2-y1)*(pressure-x1)/(x2-x1);
      found = true;
      break;
    }
  }

  if (!found) {
    if (pressure == data[n_pressures-1][0]) {
      density = data[n_pressures-1][1];
    }
    else {
      throw "Unknown xenon density for this pressure!";
    }
  }

  return density;
}

////////////////////////////////////////////////////////////////////////



