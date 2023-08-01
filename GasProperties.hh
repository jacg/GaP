#ifndef GASPROPERTIES_H
#define GASPROPERTIES_H

#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>

G4double XenonRefractiveIndex(G4double energy, G4double density);
G4double GXeScintillation(G4double energy, G4double pressure);
G4double GXeDensity(G4double pressure);


#endif
