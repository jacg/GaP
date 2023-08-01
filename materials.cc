#include "nain4.hh"
#include "g4-mandatory.hh"
#include "materials.hh"
#include "PropertiesTables.hh"
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

G4Material* peek_with_properties() {
  auto peek = n4::material_from_elements_N("peek", 1.30*g/cm3, kStateSolid,{{"H", 12},{"C" , 18},{"O", 3}});
  peek -> SetMaterialPropertiesTable(peek_properties()) ;
  return peek;
}

////////////////////////////////////////////////////////////////////////

G4Material* quartz_with_properties() {
  auto quartz = n4::material("G4_SILICON_DIOXIDE");  
  quartz -> SetMaterialPropertiesTable(quartz_properties());
  return quartz;
}

////////////////////////////////////////////////////////////////////////

G4Material* TPB_with_properties() {
  auto tpb = n4::material_from_elements_N("TPB", 1.*g/cm3, kStateSolid,{{"H", 22},{"C" , 28}});
  tpb -> SetMaterialPropertiesTable(TPB_properties());
  return tpb;
}

////////////////////////////////////////////////////////////////////////

G4Material* CopyMaterial(G4Material* original, const G4String& newname) {
  auto newmat = n4::material(newname);  
    if (newmat == 0) {

      G4double density     = original->GetDensity();
      G4double temperature = original->GetTemperature();
      G4double pressure    = original->GetPressure();
      G4State  state       = original->GetState();
      G4int    n_elem      = original->GetNumberOfElements();

      if (n_elem == 1) {
        G4double z = original->GetZ();
        G4double a = original->GetA();
        newmat = new G4Material(newname, z, a, density, state, temperature, pressure);
      }
      else {
        const G4double* fractions = original->GetFractionVector();
        newmat = new G4Material(newname, density, n_elem, state, temperature, pressure);
        for (G4int i = 0; i < n_elem; ++i)
          newmat->AddElement(new G4Element(original->GetElement(i)->GetName(),
                                          original->GetElement(i)->GetSymbol(),
                                          original->GetElement(i)->GetZ(),
                                          original->GetElement(i)->GetA()),
                                          fractions[i]);
        }
    }

    return newmat;   
}

////////////////////////////////////////////////////////////////////////

G4Material* FakeDielectric_with_properties(G4Material* model_mat, G4String name,
                      G4double pressure,
                                      G4double temperature,
                                      G4double transparency,
                                      G4double thickness,
                                      G4int    sc_yield,
                                      G4double e_lifetime,
                                      G4double photoe_p) {
  
  auto mesh_mat = CopyMaterial(model_mat, name);  
  mesh_mat -> SetMaterialPropertiesTable (FakeDielectric_properties(pressure,
                                        temperature,
                                        transparency,
                                        thickness,
                                        sc_yield,
                                        e_lifetime,
                                        photoe_p));    
 return mesh_mat;
 
}               
                 
////////////////////////////////////////////////////////////////////////
 
 G4Material* GAr_with_properties(G4double pressure, G4double temperature, G4double sc_yield, G4double e_lifetime){

  G4double density = 1.60279 * (pressure / bar) * (300 * kelvin) / temperature * kg/m3; // 1.60729 is the pressure at 1 bar and 300K

  auto GAr = n4::material_from_elements_N("GAr", density, kStateGas,{{"Ar", 1}});
  GAr -> SetMaterialPropertiesTable(GAr_properties(sc_yield, e_lifetime)) ;
  return GAr; 

}

////////////////////////////////////////////////////////////////////////



