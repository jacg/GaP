#ifndef PROPERTIESTABLES_H
#define PROPERTIESTABLES_H

#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>
#include "GasProperties.hh"


G4MaterialPropertiesTable* peek_properties();
G4MaterialPropertiesTable* quartz_properties();
G4MaterialPropertiesTable* TPB_properties();
G4MaterialPropertiesTable* GXe_properties(G4double pressure,
                               G4double temperature,
                               G4int    sc_yield,
                               G4double e_lifetime);
G4MaterialPropertiesTable* GAr_properties(G4double sc_yield, G4double e_lifetime);                                                         
G4MaterialPropertiesTable* FakeDielectric_properties(G4double pressure,
                                    G4double temperature,
                                       G4double transparency,
                                       G4double thickness,
                                       G4int    sc_yield,
                                       G4double e_lifetime,
                                       G4double photoe_p);
                              


#endif
