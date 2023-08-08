#include "nain4.hh"
#include "g4-mandatory.hh"
#include "n4_ui.hh"
#include "n4-utils.hh"
#include "n4-volumes.hh"

#include "ParticleGenerator.hh"
#include "kr83.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <FTFP_BERT.hh>
#include <G4DecayPhysics.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4Event.hh>
#include <G4IonTable.hh>
#include <G4OpticalPhysics.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4RandomDirection.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>

#include "G4OpticalPhoton.hh"

G4double vessel_out_rad_    = 288./2  *mm;
G4double vessel_out_length_ = 46.679  *cm;
G4double angle =  0 * rad;
G4double source_pos_z = 0 * cm; //range: [-vessel_out_length_/2,vessel_out_length_]
G4double cathode_z = 4.505*mm; //cathode surface, not center

void add_particle_to_vertex(G4PrimaryVertex* vertex, G4ParticleDefinition* particle, G4double energy) {
    auto primary = new G4PrimaryParticle(particle);
    primary -> SetKineticEnergy(energy);
    primary -> SetMomentumDirection(G4RandomDirection());
    // primary -> SetCharge(...)  if needed eventually
    vertex -> SetPrimary(primary);
};

std::unique_ptr<G4PrimaryVertex> generate_vertex(G4ThreeVector pos) {
    auto vertex = std::make_unique<G4PrimaryVertex>();
    vertex -> SetPosition(pos.x(), pos.y(), pos.z());
    return vertex;
};

void generate_particles_in_event(
    G4Event* event,
    std::function<G4ThreeVector()> generate_position,
    std::vector<std::tuple<G4ParticleDefinition*, G4double>> const & particles_and_energies
) {
    auto vertex = generate_vertex(generate_position());
    for (auto [particle, energy] : particles_and_energies) {
        add_particle_to_vertex(vertex.get(), particle, energy);
    }
    event -> AddPrimaryVertex(vertex.release());
}
void generate_particles_in_event(
    G4Event* event,
    G4ThreeVector position,
    std::vector<std::tuple<G4ParticleDefinition*, G4double>> const & particles_and_energies
) { return generate_particles_in_event(event, [&position] { return position; }, particles_and_energies); }

void banner(G4String msg) {
    G4cout << "********************************** :) " << msg << " (: **********************************" << G4endl;
}

void generate_Co57(G4Event* event, G4ThreeVector position, G4double /*time*/){

    //G4double lifetime = 271.8 * year;
    banner("time");

    auto gamma          = n4::find_particle("gamma");

    auto p_xray_122 = 0.86;
    auto p_no_xray_122  = 1 - p_xray_122;
    static auto distribution_1 = n4::random::biased_choice{{p_xray_122, p_no_xray_122}};

    auto p_xray_136 = 0.11;
    auto p_no_xray_136  = 1 - p_xray_136;
    static auto distribution_2 = n4::random::biased_choice{{p_xray_136, p_no_xray_136}};


    auto random_event_1 = distribution_1();
    auto random_event_2 = distribution_2();

    if (random_event_1 == 0) { banner("DECAY 1"); generate_particles_in_event(event, position, {{gamma, 122 * keV}}); }
    if (random_event_2 == 0) { banner("DECAY 2"); generate_particles_in_event(event, position, {{gamma, 136 * keV}}); }
}

void generate_Ba133(G4Event* event, G4ThreeVector position, G4double /*time*/){

    //G4double lifetime = XXX * year;
    //banner("time")
    auto gamma = n4::find_particle("gamma");

    auto p_xray = 0.99;
    auto p_no_xray = 1 - p_xray;
    static auto distribution_1 = n4::random::biased_choice{{p_xray, p_no_xray}};

    auto p_gamma_356 = 0.62;
    auto p_no_gamma_356 = 1 - p_gamma_356;
    static auto distribution_2 = n4::random::biased_choice{{p_gamma_356, p_no_gamma_356}};

    auto p_gamma_81 = 0.34;
    auto p_no_gamma_81 = 1 - p_gamma_81;
    static auto distribution_3 = n4::random::biased_choice{{p_gamma_81, p_no_gamma_81}};

    auto p_gamma_303 = 0.18;
    auto p_no_gamma_303 = 1 - p_gamma_303;
    static auto distribution_4 = n4::random::biased_choice{{p_gamma_303, p_no_gamma_303}};


    auto random_event_1 = distribution_1();
    auto random_event_2 = distribution_2();
    auto random_event_3 = distribution_3();
    auto random_event_4 = distribution_4();

    if (random_event_1 == 0) { /*banner("DECAY 1");*/ generate_particles_in_event(event, position, {{gamma,  31 * keV}}); }
    if (random_event_2 == 0) { /*banner("DECAY 2");*/ generate_particles_in_event(event, position, {{gamma, 356 * keV}}); }
    if (random_event_3 == 0) { /*banner("DECAY 3");*/ generate_particles_in_event(event, position, {{gamma,  81 * keV}}); }
    if (random_event_4 == 0) { /*banner("DECAY 4");*/ generate_particles_in_event(event, position, {{gamma, 303 * keV}}); }
}

void generate_ion_decay(G4Event* event, G4ThreeVector position, G4double /*time*/){
    std::string IonName{"Kr83m"};

    G4int A, Z;
    G4double E = 0;
    //G4double charge = 0 * eplus;
    auto T = 0*keV;

    if      (IonName == "Co57" ) { Z = 27; A =  57; }
    else if (IonName == "Ba133") { Z = 56; A = 133; }
    else if (IonName == "Am241") { Z = 95; A = 241; }
    else if (IonName == "Fe55" ) { Z = 26; A =  55; }
    else if (IonName == "Kr83m") { Z = 36; A =  83; E = 41.557 * keV; }
    else { throw "Unknown ion name: " + IonName; }

    G4ParticleDefinition* ion = G4IonTable:: GetIonTable() -> GetIon(Z, A, E);

    // TODO: is the charge really needed?
    generate_particles_in_event(event, position, {{ion, T}});
}

std::vector<std::tuple<G4ParticleDefinition*, G4double>> generate_particles_and_energies_tuples(){
	std::vector<std::tuple<G4ParticleDefinition*, G4double>> particles_and_energies;
	G4ParticleDefinition* opticalPhoton = G4OpticalPhoton::Definition();
	particles_and_energies.push_back(std::make_tuple(opticalPhoton, 100.0 * keV));
	
    return particles_and_energies;
}
	



