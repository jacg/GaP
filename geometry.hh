#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <G4Types.hh>
#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>

void place_mesh_holder_in(G4LogicalVolume* vessel);

G4PVPlacement* geometry();


#endif
