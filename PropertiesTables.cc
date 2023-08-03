#include "GasProperties.hh"
#include "PropertiesTables.hh"

#include "nain4.hh"
#include "g4-mandatory.hh"
#include "n4-constants.hh"
#include "n4-utils.hh"

#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4OpticalPhysics.hh>
#include <G4SystemOfUnits.hh>
#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>

using namespace CLHEP;

using vecd = std::vector<G4double>;

const G4double optPhotMinE_   =  0.2  * eV;
const G4double optPhotMaxE_   = 11.5  * eV;
const G4double optPhotMaxWL_  = optPhotMinE_ * nm / c4::hc;
const G4double optPhotMinWL_  = optPhotMaxE_ * nm / c4::hc;
const G4double noAbsLength_   = 1.e8  * m;
const vecd     optPhotRangeE_ = {optPhotMinE_, optPhotMaxE_};

////////////////////////////////////////////////////////////////////////
G4MaterialPropertiesTable* peek_properties(){ return n4::material_properties().done(); }


////////////////////////////////////////////////////////////////////////
G4MaterialPropertiesTable* GXe_properties(G4double pressure,
                         [[maybe_unused]] G4double temperature,
                                          G4int    sc_yield,
                                          G4double e_lifetime){

  // REFRACTIVE INDEX
  const size_t ri_entries = 200;
  auto ri_energy = n4::linspace(optPhotMinE_, optPhotMaxE_, ri_entries);
  auto density   = GXeDensity(pressure);
  auto rIndex    = n4::map<G4double>([=] (auto e) { return XenonRefractiveIndex(e, density); }
                                    , ri_energy);

  // for (size_t i=0; i<ri_entries; i++) {
  //   G4cout << "* GXe rIndex:  " << std::setw(7)
  //          << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
  // }


  // EMISSION SPECTRUM
  // Sampling from ~150 nm to 200 nm <----> from 6.20625 eV to 8.20625 eV
  const size_t sc_entries = 200;
  auto sc_energy = n4::linspace(6.20625 * eV, 8.20625 * eV, sc_entries);
  auto intensity = n4::map<G4double>([=] (auto e) { return GXeScintillation(e, pressure); },
                                     sc_energy);
  //for (int i=0; i<sc_entries; i++) {
  //  G4cout << "* GXe Scint:  " << std::setw(7) << sc_energy[i]/eV
  //         << " eV -> " << intensity[i] << G4endl;
  //}


  return n4::material_properties()
    .add("RINDEX"                    , ri_energy, rIndex)
    .add("ABSLENGTH"                 , optPhotRangeE_, noAbsLength_)
    .add("SCINTILLATIONCOMPONENT1"   , sc_energy, intensity) // Not sure if this makes sense
    .add("SCINTILLATIONCOMPONENT2"   , sc_energy, intensity) // Not sure if this makes sense
    .NEW("ELSPECTRUM"                , sc_energy, intensity)
    .add("SCINTILLATIONYIELD"        , sc_yield)
    .add("RESOLUTIONSCALE"           ,   1.0)
    .add("SCINTILLATIONTIMECONSTANT1",   4.5 * ns)
    .add("SCINTILLATIONTIMECONSTANT2", 100.  * ns)
    .add("SCINTILLATIONYIELD1"       ,    .1)
    .add("SCINTILLATIONYIELD2"       ,    .9)
    .NEW("ATTACHMENT"                , e_lifetime)
    .done();

}

////////////////////////////////////////////////////////////////////////
G4MaterialPropertiesTable* quartz_properties(){

  // REFRACTIVE INDEX
  // The range is chosen to be up to ~10.7 eV because Sellmeier's equation
  // for fused silica is valid only in that range
  const size_t ri_entries = 200;
  auto ri_energy = n4::linspace(optPhotMinE_, optPhotMaxE_, ri_entries);

  // The following values for the refractive index have been calculated
  // using Sellmeier's equation:
  //    n^2 - 1 = B_1 * \lambda^2 / (\lambda^2 - C_1) +
  //            + B_2 * \lambda^2 / (\lambda^2 - C_2) +
  //            + B_3 * \lambda^2 / (\lambda^2 - C_3),
  // with wavelength \lambda in micrometers and
  //    B_1 = 4.73E-1, B_2 = 6.31E-1, B_3 = 9.06E-1
  //    C_1 = 1.30E-2, C_2 = 4.13E-3, C_3 = 9.88E+1.

  auto ref_index_sellmeier = [] (auto e) {
    auto B_1 = 4.73e-1;
    auto B_2 = 6.31e-1;
    auto B_3 = 9.06e-1;
    auto C_1 = 1.30e-2;
    auto C_2 = 4.13e-3;
    auto C_3 = 9.88e+1;

    auto lambda  = c4::hc / e / nm * um; // in micron
    auto lambda2 = std::pow(lambda, 2);
    auto n2 = 1 + lambda2 * ( B_1 / (lambda2 - C_1)
                            + B_2 / (lambda2 - C_2)
                            + B_3 / (lambda2 - C_3));
    return std::sqrt(n2);
  };

  auto rIndex = n4::map<G4double>(ref_index_sellmeier, ri_energy);

  // for (int i=0; i<ri_entries; i++) {
  //   G4cout << "* FusedSilica rIndex:  " << std::setw(5) << ri_energy[i]/eV
  //         << " eV -> " << rIndex[i] << G4endl;
  // }


  // ABSORPTION LENGTH
  auto abs_energy = n4::scale_by(eV, {
    optPhotMinE_ / eV,
              6.46499, 6.54000, 6.59490, 6.64000, 6.72714, 6.73828, 6.75000,
              6.82104, 6.86000, 6.88000, 6.89000, 7.00000, 7.01000, 7.01797,
              7.05000, 7.08000, 7.08482, 7.30000, 7.36000, 7.40000, 7.48000,
              7.52000, 7.58000, 7.67440, 7.76000, 7.89000, 7.93000, 8.00000,
    optPhotMaxE_ / eV
    });

  auto absLength = n4::scale_by(cm, {
    noAbsLength_ / cm,
    noAbsLength_ / cm,  200.0 ,  200.0 ,   90.0 ,   45.0 ,  45.0  ,   30.0 ,
                 24.0,   21.0 ,   20.0 ,   19.0 ,   16.0 ,  14.0  ,   13.0 ,
                  8.5,    8.0 ,    6.0 ,    1.5 ,    1.2 ,   1.0  ,     .65,
                   .4,     .37,     .32,     .28,     .22,    .215,   5e-5 ,
    5e-5
    });

  return n4::material_properties()
    .add("RINDEX"   ,  ri_energy, rIndex)
    .add("ABSLENGTH", abs_energy, absLength)
    .done();

}

////////////////////////////////////////////////////////////////////////


G4MaterialPropertiesTable* TPB_properties() {
  // WLS ABSORPTION LENGTH (Version NoSecWLS)
  // The NoSecWLS is forced by setting the WLS_absLength to a very large value
  // for wavelengths higher than 380 nm where the WLS emission spectrum starts.
  auto infinite = noAbsLength_ / nm; // ~6200 nm
  auto WLS_abs_energy = n4::factor_over(c4::hc/nm, {optPhotMaxWL_, 380, 370, 360, 330, 320, 310, 300, 270, 250, 230, 210, 190, 170, 150, optPhotMinWL_});
  auto WLS_absLength  = n4::scale_by   (       nm, {infinite, infinite,  50,  30,  30,  50,  80, 100, 100, 400, 400, 350, 250, 350, 400, 400          });

  //for (int i=0; i<WLS_abs_energy.size(); i++)
  //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
  //         << " eV  ==  "              << std::setw(8) << (c4::hc / WLS_abs_energy[i]) / nm
  //         << " nm  ->  "              << std::setw(6) << WLS_absLength[i] / nm << " nm" << G4endl;

  // WLS EMISSION SPECTRUM
  // Implemented with formula (7), with parameter values in table (3)
  // Sampling from ~380 nm to 600 nm <--> from 2.06 to 3.26 eV
  const G4int WLS_emi_entries = 120;
  vecd WLS_emi_energy(WLS_emi_entries);
  for (int i=0; i<WLS_emi_entries; i++) { WLS_emi_energy.push_back( (2.06 + 0.01 * i) * eV); }

  auto tpb_emission_spectrum = [] (G4double e) {
     auto A      =   0.782;
     auto alpha  =   0.037;
     auto sigma1 =  15.43 ; // in nm
     auto mu1    = 418.10 ; // in nm
     auto sigma2 =   9.72 ; // in nm
     auto mu2    = 411.2  ; // in nm

     auto wl       = c4::hc / e / nm;
     auto exponent = alpha * (mu1 + alpha * pow(sigma1, 2)/2 - wl);
     auto erfc_arg = (mu1 + alpha * pow(sigma1, 2) - wl) / (sqrt(2) * sigma1);
     auto gaussian = 1 / sqrt(CLHEP::twopi) / sigma2
                   * exp(-pow(wl - mu2, 2) / 2 / pow(sigma2, 2));

     return A * alpha/2 * exp(exponent) * erfc(erfc_arg) + (1-A) * gaussian;
  };

  vecd WLS_emiSpectrum(WLS_emi_energy.size());
  std::transform(begin(WLS_emi_energy), end(WLS_emi_energy), begin(WLS_emiSpectrum), tpb_emission_spectrum);

  // G4cout << "* TPB WLSemi:  " << std::setw(4)
  //        << wl << " nm -> " << WLS_emiSpectrum[i] << G4endl;

  return n4::material_properties()
    .add("RINDEX"               , optPhotRangeE_, 1.67)
    // Assuming no absorption except WLS
    // .add("ABSLENGTH"            , abs_energy, abs_energy) // ??????????????????????
    .add("ABSLENGTH"            , optPhotRangeE_, noAbsLength_) // Wrong argument replaced (see previous line)
    .add("WLSABSLENGTH"         , WLS_abs_energy, WLS_absLength)
    .add("WLSCOMPONENT"         , WLS_emi_energy, WLS_emiSpectrum)
    .add("WLSTIMECONSTANT"      , 1.2 * ns) // WLS Delay
    .add("WLSMEANNUMBERPHOTONS" , 0.65) // WLS Quantum Efficiency
    // According to the paper, the QE of TPB depends on the incident wavelength.
    // As Geant4 doesn't allow this possibility, it is set to the value corresponding
    // to Xe scintillation spectrum peak.
    .done();
}

///////////////////////////////////////////////////////////////////////

G4MaterialPropertiesTable* FakeDielectric_properties(G4double pressure,
                                                     G4double temperature,
                                                     G4double transparency,
                                                     G4double thickness,
                                                     G4int    sc_yield,
                                                     G4double e_lifetime,
                                                     G4double photoe_p){
   // ABSORPTION LENGTH
   G4double abs_length = -thickness/log(transparency);

   // PHOTOELECTRIC REEMISSION
   // https://aip.scitation.org/doi/10.1063/1.1708797
   G4double stainless_wf = 4.3 * eV; // work function

   G4MaterialPropertiesTable* xenon_pt = GXe_properties(pressure, temperature, sc_yield, e_lifetime);
   return n4::material_properties()
     .add("ABSLENGTH"                   , optPhotRangeE_, abs_length)
     .NEW("WORK_FUNCTION"               ,                 stainless_wf)
     .NEW("OP_PHOTOELECTRIC_PROBABILITY",                 photoe_p)
     //     .copy_from(xenon_pt, {"cosas"})
     .done();

   //return n4::material_properties()
    //.add("ABSLENGTH", abs_energy, absLength)
    //.add("WORK_FUNCTION", stainless_wf) //TRUE NEEDED
    //.add("OP_PHOTOELECTRIC_PROBABILITY", photoe_p, true) //TRUE NEEDED


    //.add("RINDEX", prop_vec->GetEnergyVector(), prop_vec->GetPropertyVector() )
    //.add("SCINTILLATIONCOMPONENT1", xenon_pt->GetProperty("SCINTILLATIONCOMPONENT1"))
    //.add("SCINTILLATIONCOMPONENT2", xenon_pt->GetProperty("SCINTILLATIONCOMPONENT2"))
    //.add("SCINTILLATIONYIELD", xenon_pt->GetProperty("SCINTILLATIONYIELD"))
    //.add("RESOLUTIONSCALE", xenon_pt->GetProperty("RESOLUTIONSCALE"))
    //.add("SCINTILLATIONTIMECONSTANT1", xenon_pt->GetProperty("SCINTILLATIONTIMECONSTANT1"))
    //.add("SCINTILLATIONTIMECONSTANT2", xenon_pt->GetProperty("SCINTILLATIONTIMECONSTANT2"))
    //.add("SCINTILLATIONYIELD1", xenon_pt->GetProperty("SCINTILLATIONYIELD1"))
    //.add("SCINTILLATIONYIELD2", xenon_pt->GetProperty("SCINTILLATIONYIELD2"))
    //.add("ATTACHMENT", xenon_pt->GetProperty("ATTACHMENT"),1)  //TRUE NEEDED
    //.done();

}

////////////////////////////////////////////////////////////////////////

G4MaterialPropertiesTable* GAr_properties(G4double sc_yield, G4double e_lifetime){

    //e_lifetime = 1000.*ms;

    // An argon gas proportional scintillation counter with UV avalanche photodiode scintillation
    // readout C.M.B. Monteiro, J.A.M. Lopes, P.C.P.S. Simoes, J.M.F. dos Santos, C.A.N. Conde

    // REFRACTIVE INDEX
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;


    vecd abs_energy_Ar;
    vecd absLength_Ar;
    vecd sc_energy_Ar;
    vecd intensity_Ar;
    const G4int sc_entries = 380;


    vecd ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    vecd rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double wl = h_Planck * c_light / ri_energy[i] * 1000; // in micron
      // From refractiveindex.info
      rIndex.push_back(1 + 0.012055*(0.2075*pow(wl,2)/(91.012*pow(wl,2)-1) +
                                     0.0415*pow(wl,2)/(87.892*pow(wl,2)-1) +
                                     4.3330*pow(wl,2)/(214.02*pow(wl,2)-1)));
      //G4cout << "* GAr rIndex:  " << std::setw(5) << ri_energy[i]/eV
      //       << " eV -> " << rIndex[i] << G4endl;


    // ABSORPTION LENGTH
    abs_energy_Ar = {optPhotMinE_, optPhotMaxE_};
    absLength_Ar  = {noAbsLength_, noAbsLength_};


     // EMISSION SPECTRUM
//    G4double Wavelength_peak  = 128.000 * nm;
    G4double Wavelength_peak  = 128.000 * nm; // Xe, to be changed back
    G4double Wavelength_sigma =   2.929 * nm;
    G4double Energy_peak  = (h_Planck*c_light / Wavelength_peak);
    G4double Energy_sigma = (h_Planck*c_light * Wavelength_sigma / pow(Wavelength_peak,2));
    //G4cout << "*** GAr Energy_peak: " << Energy_peak/eV << " eV   Energy_sigma: "
    //       << Energy_sigma/eV << " eV" << G4endl;

    // Sampling from ~110 nm to 150 nm <----> from ~11.236 eV to 8.240 eV

    sc_energy_Ar;
    intensity_Ar;
    for (int i=0; i<sc_entries; i++){
      sc_energy_Ar.push_back(8.240*eV + 0.008*i*eV);
      intensity_Ar.push_back(exp(-pow(Energy_peak/eV-sc_energy_Ar[i]/eV,2) /
                              (2*pow(Energy_sigma/eV, 2)))/(Energy_sigma/eV*sqrt(pi*2.)));
      //G4cout << "* GAr energy: " << std::setw(6) << sc_energy[i]/eV << " eV  ->  "
      //       << std::setw(6) << intensity[i] << G4endl;
    }

    }

    return n4::material_properties()
        .add("RINDEX", ri_energy, rIndex)
        .add("ABSLENGTH", abs_energy_Ar, absLength_Ar)
        .add("SCINTILLATIONCOMPONENT1", sc_energy_Ar, intensity_Ar)
        .add("SCINTILLATIONCOMPONENT2", sc_energy_Ar, intensity_Ar)
        .NEW("ELSPECTRUM",  sc_energy_Ar, intensity_Ar)
        .add("SCINTILLATIONYIELD", sc_yield)
        .add("SCINTILLATIONTIMECONSTANT1",6.*ns)  // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
        .add("SCINTILLATIONTIMECONSTANT2", 3480.*ns) // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
        .add("SCINTILLATIONYIELD1", .136) // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
        .add("SCINTILLATIONYIELD2", .864) // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
        .add("RESOLUTIONSCALE", 1.0)
        .NEW("ATTACHMENT", e_lifetime)
    .done();
}
