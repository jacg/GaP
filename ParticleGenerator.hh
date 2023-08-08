#ifndef PARTICLEGEN_H
#define PARTICLEGEN_H

void banner(G4String msg);
void generate_Co57(G4Event* event, G4ThreeVector position, G4double /*time*/);
void generate_Ba133(G4Event* event, G4ThreeVector position, G4double /*time*/);
void generate_ion_decay(G4Event* event, G4ThreeVector position, G4double /*time*/);
std::vector<std::tuple<G4ParticleDefinition*, G4double>> generate_partilces_and_energies_tuples();
void generate_particles_in_event(
    G4Event* event,
    std::function<G4ThreeVector()> generate_position,
    std::vector<std::tuple<G4ParticleDefinition*, G4double>> const & particles_and_energies
);
void generate_particles_in_event(
    G4Event* event,
    G4ThreeVector position,
    std::vector<std::tuple<G4ParticleDefinition*, G4double>> const & particles_and_energies
);

#endif
