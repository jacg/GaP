#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <G4Types.hh>
#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>

void build_mesh_holder(G4double meshBracket_rad_, G4LogicalVolume* vessel, G4Material* peek, G4Material* steel);

G4PVPlacement* geometry();


#endif
