#ifndef KR83_H
#define KR83_H

G4ThreeVector random_generator_inside_drift(std::optional<G4double> fixed_z);
void kr83_generator(G4Event* event, G4double energy_32_ , G4double energy_9_ , G4double probGamma_9_, G4double lifetime_9_);


#endif
