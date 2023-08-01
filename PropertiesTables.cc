#include "nain4.hh"
#include "g4-mandatory.hh"
#include "GasProperties.hh"
#include "PropertiesTables.hh"

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

G4double optPhotMinE_ =  0.2  * eV;
G4double optPhotMaxE_ = 11.5  * eV;
G4double noAbsLength_ = 1.e8  * m;


////////////////////////////////////////////////////////////////////////

G4MaterialPropertiesTable* peek_properties(){

    return n4::material_properties()
    .done();   
    
}
////////////////////////////////////////////////////////////////////////

G4MaterialPropertiesTable* GXe_properties(G4double pressure,
                               G4double temperature,
                               G4int    sc_yield,
                               G4double e_lifetime){
   
    // REFRACTIVE INDEX
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;
    
    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    G4double density = GXeDensity(pressure);
    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double ri = XenonRefractiveIndex(ri_energy[i], density);
      rIndex.push_back(ri);
      // G4cout << "* GXe rIndex:  " << std::setw(7)
      //        << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    
    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    
    // EMISSION SPECTRUM
    // Sampling from ~150 nm to 200 nm <----> from 6.20625 eV to 8.20625 eV
    const G4int sc_entries = 200;
    std::vector<G4double> sc_energy;
    for (int i=0; i<sc_entries; i++){
      sc_energy.push_back(6.20625 * eV + 0.01 * i * eV);
    }
    std::vector<G4double> intensity;
    for (G4int i=0; i<sc_entries; i++) {
      G4double inten = GXeScintillation(sc_energy[i], pressure);
      intensity.push_back(inten);
    }
    //for (int i=0; i<sc_entries; i++) {
    //  G4cout << "* GXe Scint:  " << std::setw(7) << sc_energy[i]/eV
    //         << " eV -> " << intensity[i] << G4endl;
    //}    
    
    
    return n4::material_properties()
        .add("RINDEX", ri_energy, rIndex)
        .add("ABSLENGTH", abs_energy, absLength)
        .add("SCINTILLATIONCOMPONENT1", sc_energy, intensity)
        .add("SCINTILLATIONCOMPONENT2", sc_energy, intensity)
        .add("ELSPECTRUM", sc_energy, intensity, 1)   //Pq1?
        .add("SCINTILLATIONYIELD", sc_yield)
        .add("RESOLUTIONSCALE", 1.0)
        .add("SCINTILLATIONTIMECONSTANT1", 4.5 * ns)
        .add("SCINTILLATIONTIMECONSTANT2", 100. * ns)
        .add("SCINTILLATIONYIELD1", .1)
        .add("SCINTILLATIONYIELD2", .9)
        .add("ATTACHMENT", e_lifetime, 1)
        .done();

}

////////////////////////////////////////////////////////////////////////

G4MaterialPropertiesTable* quartz_properties(){

// REFRACTIVE INDEX
  // The range is chosen to be up to ~10.7 eV because Sellmeier's equation
  // for fused silica is valid only in that range

  const G4int ri_entries = 200;
  G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

  std::vector<G4double> ri_energy;
  for (int i=0; i<ri_entries; i++) {
    ri_energy.push_back(optPhotMinE_ + i * eWidth);
  }

  // The following values for the refractive index have been calculated
  // using Sellmeier's equation:
  //    n^2 - 1 = B_1 * \lambda^2 / (\lambda^2 - C_1) +
  //            + B_2 * \lambda^2 / (\lambda^2 - C_2) +
  //            + B_3 * \lambda^2 / (\lambda^2 - C_3),
  // with wavelength \lambda in micrometers and
  //    B_1 = 4.73E-1, B_2 = 6.31E-1, B_3 = 9.06E-1
  //    C_1 = 1.30E-2, C_2 = 4.13E-3, C_3 = 9.88E+1.

  G4double B_1 = 4.73e-1;
  G4double B_2 = 6.31e-1;
  G4double B_3 = 9.06e-1;
  G4double C_1 = 1.30e-2;
  G4double C_2 = 4.13e-3;
  G4double C_3 = 9.88e+1;

  std::vector<G4double> rIndex;
  for (int i=0; i<ri_entries; i++) {
    G4double lambda = h_Planck*c_light/ri_energy[i]*1000; // in micron
    G4double n2 = 1 + B_1*pow(lambda,2)/(pow(lambda,2)-C_1)
      + B_2*pow(lambda,2)/(pow(lambda,2)-C_2)
      + B_3*pow(lambda,2)/(pow(lambda,2)-C_3);
    rIndex.push_back(sqrt(n2));
    // G4cout << "* FusedSilica rIndex:  " << std::setw(5) << ri_energy[i]/eV
    //       << " eV -> " << rIndex[i] << G4endl;
  }


  // ABSORPTION LENGTH
  std::vector<G4double> abs_energy = {
    optPhotMinE_,  6.46499 * eV,
    6.54000 * eV,  6.59490 * eV,  6.64000 * eV,  6.72714 * eV,
    6.73828 * eV,  6.75000 * eV,  6.82104 * eV,  6.86000 * eV,
    6.88000 * eV,  6.89000 * eV,  7.00000 * eV,  7.01000 * eV,
    7.01797 * eV,  7.05000 * eV,  7.08000 * eV,  7.08482 * eV,
    7.30000 * eV,  7.36000 * eV,  7.40000 * eV,  7.48000 * eV,
    7.52000 * eV,  7.58000 * eV,  7.67440 * eV,  7.76000 * eV,
    7.89000 * eV,  7.93000 * eV,  8.00000 * eV,
    optPhotMaxE_
  };

  std::vector<G4double> absLength = {
    noAbsLength_, noAbsLength_,
    200.0 * cm,   200.0 * cm,  90.0 * cm,  45.0 * cm,
    45.0 * cm,    30.0 * cm,  24.0 * cm,  21.0 * cm,
    20.0 * cm,    19.0 * cm,  16.0 * cm,  14.0 * cm,
    13.0 * cm,     8.5 * cm,   8.0 * cm,   6.0 * cm,
    1.5 * cm,     1.2 * cm,   1.0 * cm,   .65 * cm,
     .4 * cm,     .37 * cm,   .32 * cm,   .28 * cm,
     .22 * cm,    .215 * cm,  .00005*cm,
     .00005* cm
  };
  
   return n4::material_properties()
    .add("RINDEX", ri_energy, rIndex)
    .add("ABSLENGTH", abs_energy, absLength)
    .done();

}

////////////////////////////////////////////////////////////////////////


G4MaterialPropertiesTable* TPB_properties(){

  // REFRACTIVE INDEX
  std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
  std::vector<G4double> TPB_rIndex      = {1.67    , 1.67};
  
  // ABSORPTION LENGTH
  // Assuming no absorption except WLS
  std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
  std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
  
   // WLS ABSORPTION LENGTH (Version NoSecWLS)
   // The NoSecWLS is forced by setting the WLS_absLength to noAbsLength_
   // for wavelengths higher than 380 nm where the WLS emission spectrum starts.
   std::vector<G4double> WLS_abs_energy = {
     optPhotMinE_,
     h_Planck * c_light / (380. * nm),  h_Planck * c_light / (370. * nm),
     h_Planck * c_light / (360. * nm),  h_Planck * c_light / (330. * nm),
     h_Planck * c_light / (320. * nm),  h_Planck * c_light / (310. * nm),
     h_Planck * c_light / (300. * nm),  h_Planck * c_light / (270. * nm),
     h_Planck * c_light / (250. * nm),  h_Planck * c_light / (230. * nm),
     h_Planck * c_light / (210. * nm),  h_Planck * c_light / (190. * nm),
     h_Planck * c_light / (170. * nm),  h_Planck * c_light / (150. * nm),
     optPhotMaxE_
   };

   std::vector<G4double> WLS_absLength = {
     noAbsLength_,                 // ~6200 nm
     noAbsLength_,   50. * nm,     // 380 , 370 nm
     30. * nm,      30. * nm,      // 360 , 330 nm
     50. * nm,      80. * nm,      // 320 , 310 nm
     100. * nm,     100. * nm,     // 300 , 270 nm
     400. * nm,     400. * nm,     // 250 , 230 nm
     350. * nm,     250. * nm,     // 210 , 190 nm
     350. * nm,     400. * nm,     // 170 , 150 nm
     400. * nm                     // ~108 nm
   };

   //for (int i=0; i<WLS_abs_energy.size(); i++)
   //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
   //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
   //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / nm << " nm" << G4endl;
  
   // WLS EMISSION SPECTRUM
   // Implemented with formula (7), with parameter values in table (3)
   // Sampling from ~380 nm to 600 nm <--> from 2.06 to 3.26 eV
   const G4int WLS_emi_entries = 120;
   std::vector<G4double> WLS_emi_energy;
   for (int i=0; i<WLS_emi_entries; i++)
      WLS_emi_energy.push_back(2.06 * eV + 0.01 * i * eV);

      std::vector<G4double> WLS_emiSpectrum;
      G4double A      = 0.782;
      G4double alpha  = 3.7e-2;
      G4double sigma1 = 15.43;
      G4double mu1    = 418.10;
      G4double sigma2 = 9.72;
      G4double mu2    = 411.2;

    for (int i=0; i<WLS_emi_entries; i++) {
      G4double wl = (h_Planck * c_light / WLS_emi_energy[i]) / nm;
      WLS_emiSpectrum.push_back(A * (alpha/2.) * exp((alpha/2.) *
                          (2*mu1 + alpha*pow(sigma1,2) - 2*wl)) *
                          erfc((mu1 + alpha*pow(sigma1,2) - wl) / (sqrt(2)*sigma1)) +
                          (1-A) * (1 / sqrt(2*pow(sigma2,2)*3.1416)) *
                                exp((-pow(wl-mu2,2)) / (2*pow(sigma2,2))));
      // G4cout << "* TPB WLSemi:  " << std::setw(4)
      //        << wl << " nm -> " << WLS_emiSpectrum[i] << G4endl;
    };
   
  return n4::material_properties()
    .add("RINDEX", rIndex_energies, TPB_rIndex)
    .add("ABSLENGTH", abs_energy, abs_energy)
    .add("WLSABSLENGTH", WLS_abs_energy, WLS_absLength)
    .add("WLSCOMPONENT", WLS_emi_energy, WLS_emiSpectrum)
    .add("WLSTIMECONSTANT",  1.2 * ns)  // WLS Delay
    .add("WLSMEANNUMBERPHOTONS",  0.65)      // WLS Quantum Efficiency
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
   G4double abs_length   = -thickness/log(transparency);
   std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
   std::vector<G4double> absLength  = {abs_length, abs_length};
   
   // PHOTOELECTRIC REEMISSION
   // https://aip.scitation.org/doi/10.1063/1.1708797
   G4double stainless_wf = 4.3 * eV; // work function
   
   G4MaterialPropertiesTable* xenon_pt = GXe_properties(pressure, temperature, sc_yield, e_lifetime);
   xenon_pt ->  AddProperty("ABSLENGTH", abs_energy, absLength);
   xenon_pt ->  AddConstProperty("WORK_FUNCTION", stainless_wf, true);
   xenon_pt ->  AddConstProperty ("OP_PHOTOELECTRIC_PROBABILITY", photoe_p, true);
   
   return xenon_pt;
   
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
    
    
    std::vector<G4double> abs_energy_Ar;
    std::vector<G4double> absLength_Ar;
    std::vector<G4double> sc_energy_Ar;
    std::vector<G4double> intensity_Ar;
    const G4int sc_entries = 380;


    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
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
        .add("ELSPECTRUM",  sc_energy_Ar, intensity_Ar, sc_entries)
        .add("SCINTILLATIONYIELD", sc_yield)
        .add("SCINTILLATIONTIMECONSTANT1",6.*ns)  // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
        .add("SCINTILLATIONTIMECONSTANT2", 3480.*ns) // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
        .add("SCINTILLATIONYIELD1", .136) // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
        .add("SCINTILLATIONYIELD2", .864) // From https://dspace.mit.edu/bitstream/handle/1721.1/129347/1903.06706.pdf?sequence=2&isAllowed=y
        .add("RESOLUTIONSCALE", 1.0)
        .add("ATTACHMENT", e_lifetime, 1)
    .done();   

}   

///////////////////////////////////////////////////////////////////////








