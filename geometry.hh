#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <G4Types.hh>
#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>

#include "detector.hh"

G4LogicalVolume* get_world();

void place_mesh_holder_in(G4LogicalVolume* vessel);
void place_quartz_window_holder_in(G4LogicalVolume* vessel);
void place_pmt_holder_in(G4LogicalVolume* vessel);

struct field_cage_parameters {
  G4double cathode_z;
  G4double cathBracket_z;
  G4double drift_length;
  G4double el_length;
  G4double drift_z;
  G4double drift_r;
};

field_cage_parameters model_something_old();
field_cage_parameters model_something_new();

void place_rings_in         (G4LogicalVolume* vessel, field_cage_parameters const & fcp);
void place_anode_el_gate_in (G4LogicalVolume* vessel, field_cage_parameters const & fcp);

G4PVPlacement* geometry();


#endif
