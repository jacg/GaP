#ifndef MATERIALS_H
#define MATERIALS_H

#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>
#include "PropertiesTables.hh"


G4Material* peek_with_properties();
G4Material* quartz_with_properties();
G4Material* TPB_with_properties();
G4Material* CopyMaterial(G4Material* original, const G4String& newname);
G4Material* FakeDielectric_with_properties(G4Material* model_mat, G4String name,
									  G4double pressure,
                                      G4double temperature,
                                      G4double transparency,
                                      G4double thickness,
                                      G4int    sc_yield,
                                      G4double e_lifetime,
                                      G4double photoe_p);
G4Material* GAr_with_properties(G4double pressure, G4double temperature, G4double sc_yield, G4double e_lifetime);                       
//G4Material* GKr_with_properties(G4double pressure, G4double temperature, G4double sc_yield, G4double e_lifetime);                       


#endif
