#include "nain4.hh"
#include "g4-mandatory.hh"
#include "n4_ui.hh"
#include "n4-utils.hh"

#include "geometry.hh"
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
#include <G4RunManagerFactory.hh>
#include <G4Step.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4Types.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4UnionSolid.hh>
#include <G4VisAttributes.hh>
#include <G4VisExecutive.hh>
#include <G4VisManager.hh>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <vector>


void generate_gammas(G4Event* event, G4ThreeVector position, G4double time, G4double energy) {
    auto gamma = nain4::find_particle("gamma");
    auto p = energy*keV * G4RandomDirection();
    auto vertex = new G4PrimaryVertex{position, time};
    vertex -> SetPrimary(new G4PrimaryParticle(gamma,  p.x(),  p.y(),  p.z()));
    event -> AddPrimaryVertex(vertex);
}

void generate_electrons(G4Event* event, G4ThreeVector /*position*/, G4double /*time*/) {
    static G4ParticleGun* particleGun = new G4ParticleGun(1);
    auto electron = nain4::find_particle("e-");
    particleGun -> SetParticleDefinition(electron);
    particleGun -> SetParticleEnergy(511*keV);
    particleGun -> SetParticleMomentumDirection(G4RandomDirection());
    particleGun -> GeneratePrimaryVertex(event);
}

namespace find_me_a_good_name {
    G4double cathode_z           = (90.1125 - 15.745) *mm;
    G4double mesh_thickn_        =   0.075            *mm;
    G4double meshBracket_thickn_ =   6                *mm;
    G4double drift_length_       =  96                *mm - meshBracket_thickn_ ;
    G4double drift_z             = cathode_z - mesh_thickn_/2 - drift_length_/2;
    G4double meshBracket_rad_    = 180./2  *mm;
}

void generate_inside(G4Event* event, G4double /*time*/){
    static G4ParticleGun* particleGun = new G4ParticleGun(2);
    auto electron = nain4::find_particle("gamma");
    particleGun -> SetParticleDefinition(electron);
    particleGun -> SetParticleEnergy(511*keV);
    particleGun -> SetParticleMomentumDirection(G4RandomDirection());

    using namespace find_me_a_good_name;

    G4double r     = G4RandFlat::shoot( 0., meshBracket_rad_);
    G4double angle = G4RandFlat::shoot( 0., 2*M_PI);
    G4double z     = G4RandFlat::shoot(-drift_length_/2 + drift_z, drift_length_/2 + drift_z);

    G4double pos_x = r * cos(angle);
    G4double pos_y = r * sin(angle);
    G4double pos_z = z;
    
    particleGun -> SetParticlePosition({pos_x, pos_y, pos_z});
    particleGun -> GeneratePrimaryVertex(event);

}

void add_particle_to_vertex(G4PrimaryVertex* vertex, G4ParticleDefinition* particle, G4double energy) {
    auto primary = new G4PrimaryParticle(particle);
    primary -> SetKineticEnergy(energy);
    primary -> SetMomentumDirection(G4RandomDirection());
    // primary -> SetCharge(...)  if needed eventually
    vertex -> SetPrimary(primary);
};

G4ThreeVector position_active_volume() {
    using namespace find_me_a_good_name;
    auto r     = G4RandFlat::shoot( 0., meshBracket_rad_);
    auto angle = G4RandFlat::shoot( 0., 2*M_PI); // Flat distribution in the radius doesn't make sense
    auto z     = G4RandFlat::shoot(-drift_length_/2 + drift_z, drift_length_/2 + drift_z);
    auto pos_x = r * cos(angle);
    auto pos_y = r * sin(angle);
    auto pos_z = z;
    return G4ThreeVector{pos_x, pos_y, pos_z};
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
    G4cout << "********************************** " << msg << " **********************************" << G4endl;
}

void generate_Kr83m2_decay(G4Event* event, G4double /*time*/){

    auto electron       = n4::find_particle("e-");
    auto optical_photon = n4::find_particle("opticalphoton");
    auto gamma          = n4::find_particle("gamma");

    auto p_t1_ic_1au = 0.76;
    auto p_t1_ic_2au = 0.09;
    auto p_t1_ic_1au_xray = 1 - p_t1_ic_1au - p_t1_ic_2au;
    static auto distribution_1 = n4::random::biased_choice{{p_t1_ic_1au, p_t1_ic_2au, p_t1_ic_1au_xray}};

    // falta tener en cuenta el tiempo medio de las partículas
    auto p_t2_ic_1au = 0.95;
    auto p_t2_xray = 1 - p_t2_ic_1au;
    static auto distribution_2 = n4::random::biased_choice{{p_t2_ic_1au, p_t2_xray}};

    auto random_event_1 = distribution_1();
    auto random_event_2 = distribution_2();

    if (random_event_1 == 0) {
        banner("DECAY 1 = 0.76");
        generate_particles_in_event(event, position_active_volume, {{electron, 30 * keV},
                                                                    {electron,  2 * keV}});
    } else if (random_event_1 == 1) {
        banner("DECAY 1 = 0.09");
        generate_particles_in_event(event, position_active_volume, {{electron, 18 * keV},
                                                                    {electron, 10 * keV},
                                                                    {electron,  2 * keV},
                                                                    {electron,  2 * keV}});
    } else if (random_event_1 == 2) {
        banner("DECAY 1 = 0.15");
        generate_particles_in_event(event, position_active_volume, {{ electron      , 18 * keV},
                                                                    { optical_photon, 12 * keV},   // This energy looks suspicious for an optical photon
                                                                    { electron      ,  2 * keV}});
    }

    if (random_event_2 == 0) {
        banner("DECAY 2 = 0.95");
        generate_particles_in_event(event, position_active_volume, {{electron, 7.6 * keV},
                                                                    {electron, 1.8 * keV}});
    } else if (random_event_2 == 1) {
        banner("DECAY 2 = 0.05");
        generate_particles_in_event(event, position_active_volume, {{gamma, 9.4 * keV}});
    }
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




auto get_pre_volume_name(G4Step const * const step) {
    return step -> GetPreStepPoint() -> GetTouchableHandle() -> GetVolume() -> GetLogicalVolume() -> GetName();
};

int main(int argc, char *argv[]) {
    std::cout << "Hello World!" << std::endl;
    
    //G4double time;
    //auto get_time[](auto step){ G4double time = step -> GetDeltaTime();};
    
    // G4double vessel_out_rad_    = 288./2  *mm;
    // G4double vessel_out_length_ = 46.679  *cm;
    // G4double angle =  0 * rad;
    // G4double source_pos_z = 0 * cm; //range: [-vessel_out_length_/2,vessel_out_length_]
    // G4double cathode_z = 4.505*mm; //cathode surface, not center

    auto kr83m = [](auto event){generate_ion_decay(event, random_generator_inside_drift(), 0);};
    //auto kr83m= [](auto event){ kr83_generator(event, 32.1473*keV, 9.396*keV,  0.0490, 154.*ns); }; //From the box_source
    //auto gammas = [](auto event){ generate_gammas(event, {0., 0., 167.6775*mm + 50.*mm}, 0); }; //From the box_source
    //auto electrons = [](auto event){generate_electrons(event, {0., 0., 0.}, 0); };
    //auto electrons = [](auto event){generate_electrons(event, {initial_pos()}, 0); };
    //auto electrons = [](auto event){generate_inside(event, 0); };
    //auto Kr83m2 = [](auto event){generate_Kr83m2_decay(event, 0); };
    //auto Co57 = [vessel_out_rad_, angle, source_pos_z](auto event){generate_gammas(event, {vessel_out_rad_*cos(angle), vessel_out_rad_*sin(angle), source_pos_z}, 0, 122.);};  //From the surface
    //auto Co57 = [vessel_out_rad_, angle, source_pos_z](auto event){generate_ion_decay(event, {vessel_out_rad_*cos(angle), vessel_out_rad_*sin(angle), source_pos_z}, 0);};  //From the surface
    //auto Am241 = [cathode_z](auto event){generate_ion_decay(event, {0., 0., cathode_z}, 0);};  //From the surface of the cathode

    G4double  energy_deposit_total;
    G4double  energy_deposit_total_1;
    G4double  energy_deposit_total_2;
    G4int  counts;
    G4int  eventCounter;
    std::vector<int> trackIDVector;

    //auto get_energy = [](G4Step const* step){G4cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " << step -> GetPreStepPoint() -> GetKineticEnergy() << G4endl; }
    [[maybe_unused]]
    auto get_energy_old = [&energy_deposit_total](G4Step const* step) {
        if (step -> GetPreStepPoint() -> GetTouchableHandle() -> GetVolume() -> GetLogicalVolume() -> GetName() == "gas_drift") {
            energy_deposit_total += step -> GetTotalEnergyDeposit();
        }
        //G4cout << "********************************** " << energy_deposit_total << " **********************************" << G4endl;
    } ;
    
    [[maybe_unused]]
    auto print_energy = [& energy_deposit_total](G4Event const*){G4cout << "*******************************/// " << energy_deposit_total << " ///*******************************"  << G4endl;} ;

    const std::string& filename_event = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation_G4.txt";
    const std::string& filename_step = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation.txt";

    //const std::string& filename_event_1 = "EnergyDepositTotal_Co57_10bar_sum_withoutCompton.txt";
    const std::string& filename_event_1 = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation.txt";
    const std::string& filename_event_2 = "EnergyDepositTotal_Kr83m_10bar_sum_justCompton.txt";

    auto write_energy_event = [& energy_deposit_total, & filename_event, & eventCounter](G4Event const*) {
        if (energy_deposit_total != 0.0) {
            std::ofstream file;
            file.open (filename_event, std::ios::app);
            // file << energy_deposit_total << "\n";
            file << energy_deposit_total * 1000. << "\n";
            file.close();
    
            //G4cout << "*************************************  :)  " << energy_deposit_total << "  (:  *************************************      <---------"  << G4endl;
            eventCounter++;
            if (eventCounter % 100000 == 0) {
                G4cout << "*************************************  :)  " << eventCounter << "  (:  *************************************"  << G4endl;
            }
        }
    } ;

    [[maybe_unused]]
    auto write_energy_event_double = [& energy_deposit_total_1, & energy_deposit_total_2, & filename_event_1, & filename_event_2, & eventCounter](G4Event const*) {
        if (energy_deposit_total_1 != 0.0) {
            std::ofstream file;
            file.open (filename_event_1, std::ios::app);
            file << energy_deposit_total_1 * 1000. << "\n";
            //file << energy_deposit_total_1  << "\n";
            file.close();
        }
        if (energy_deposit_total_2 != 0.0) {
            std::ofstream file;
            file.open (filename_event_2, std::ios::app);
            file << energy_deposit_total_2 * 1000. << "\n";
            //file << energy_deposit_total_2 << "\n";
            file.close();
        }

        eventCounter++;
        if (eventCounter % 10000 == 0) {
            G4cout << "*************************************  :)  " << eventCounter << "  (:  *************************************"  << G4endl;
        }
    };

    auto write_info_and_get_energy_step = [&filename_step, &energy_deposit_total, &counts](G4Step const* step) {

        G4Track* track = step->GetTrack();
    
        if (step -> GetPreStepPoint() -> GetTouchableHandle() -> GetVolume() -> GetLogicalVolume() -> GetName() == "gas_drift") {
            // auto energy_deposit_step_pre  = step -> GetPreStepPoint()  -> GetTotalEnergy();
            // auto energy_deposit_step_post = step -> GetPostStepPoint() -> GetTotalEnergy();
            auto energy_kinetic = step->GetPreStepPoint()->GetKineticEnergy();
            //auto energy_deposit_step = energy_deposit_step_pre - energy_deposit_step_post;
        
            G4double energy_deposit_step = step -> GetTotalEnergyDeposit();
         
            //const G4Track* track = step->GetTrack();
        
            G4int PDGEncoding = track->GetDefinition()->GetPDGEncoding();
            G4int particleID = track->GetParentID();
            G4int trackID = track->GetTrackID();
            G4String particleType = track->GetDefinition()->GetParticleName();
            G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
            G4double time = step->GetPostStepPoint()->GetGlobalTime();
            time = time * ns;
        
            G4double stepLenght = step->GetStepLength();
            auto velocity = step->GetPreStepPoint()->GetVelocity();
            auto  speed = step->GetPreStepPoint()->GetBeta() *  CLHEP::c_light;
            const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
            G4String interactionType = process->GetProcessName();
        
            if (energy_deposit_step != 0.0 && interactionType!= "Transportation") {
                if (energy_deposit_step != 0.0) {
                    energy_deposit_total += energy_deposit_step;
                    counts++;
            
                    std::ofstream file(filename_step, std::ios::app);
                    int colWidth = 20;
                    file << std::left << std::setw(colWidth) << PDGEncoding ;
                    file << std::left << std::setw(colWidth) << particleID;
                    file << std::left << std::setw(colWidth) << trackID;
                    file << std::left << std::setw(colWidth) << particleType;
                    file << std::left << std::setw(colWidth) << time;
                    file << std::left << std::setw(colWidth) << stepLenght;
                    file << std::left << std::setw(colWidth) << velocity;
                    file << std::left << std::setw(colWidth) << speed;
                    file << std::left << std::setw(colWidth) << position.x();
                    file << std::left << std::setw(colWidth) << position.y();
                    file << std::left << std::setw(colWidth) << position.z();
                    file << std::left << std::setw(colWidth) << energy_kinetic ;
                    file << std::left << std::setw(colWidth) << energy_deposit_step ;
                    file << std::left << std::setw(colWidth) << interactionType;
                    file << std::endl;
                    file.close();
                }

            }
        }
    };

    [[maybe_unused]]
    auto get_energy_double = [&energy_deposit_total_1](G4Step const* step) {
        if (get_pre_volume_name(step) == "gas_drift") {
            auto energy_deposit_step  = step -> GetTotalEnergyDeposit();
            const G4VProcess* process = step -> GetPostStepPoint() -> GetProcessDefinedStep();
            G4String interactionType = process -> GetProcessName();

            if (interactionType != "Transportation") { energy_deposit_total_1 += energy_deposit_step; }
        }
    };

    [[maybe_unused]]
    auto get_energy = [&energy_deposit_total](G4Step const* step){
        if (get_pre_volume_name(step) == "gas_drift") {
            auto energy_deposit_step  = step -> GetTotalEnergyDeposit();
            const G4VProcess* process = step -> GetPostStepPoint() -> GetProcessDefinedStep();
            G4String interactionType = process -> GetProcessName();

            if (interactionType != "Transportation") { energy_deposit_total += energy_deposit_step; }
        }
    };

    [[maybe_unused]]
    auto get_energy_and_check_track= [&energy_deposit_total, &trackIDVector](G4Step const* step) {
        if (get_pre_volume_name(step) == "gas_drift") {
            auto energy_deposit_step = step -> GetTotalEnergyDeposit();
            const G4VProcess* process = step -> GetPostStepPoint() -> GetProcessDefinedStep();
            G4String interactionType = process -> GetProcessName();

            G4Track* track = step -> GetTrack();
            G4int trackID = track -> GetTrackID();
            auto trackID_check = std::find(trackIDVector.begin(), trackIDVector.end(), trackID);

            if (trackID_check != trackIDVector.end()) {
                //G4cout << "******************************* " << "El número " << trackID << " se encuentra en el vector." << "*******************************"  << G4endl;
            } else {
                //G4cout << "******************************* " << "El número " << trackID << " NO se encuentra en el vector (INFO GUARDADA)" << "*******************************"  << G4endl;
                if(interactionType != "Transportation"){
                    energy_deposit_total += energy_deposit_step;
                }
            }
        }
    };

    [[maybe_unused]]
    auto get_energy_compton_double = [&energy_deposit_total_1, &energy_deposit_total_2, &trackIDVector](G4Step const* step) {
        if (get_pre_volume_name(step) == "gas_drift") {
            auto energy_deposit_step = step  -> GetTotalEnergyDeposit();
            const G4VProcess* process = step  -> GetPostStepPoint() -> GetProcessDefinedStep();
            G4String interactionType = process->GetProcessName();
            bool savedEnergyStep = false;

            if(interactionType != "Transportation") {

                //G4cout << "-----------------------------------------------------------------------------------------------------------------------------"  << G4endl;

                G4Track* track = step->GetTrack();

                if (interactionType == "compt") {
                    G4int trackID = track->GetTrackID();
                    trackIDVector.push_back(trackID);
                    //G4cout << "******************************* COMPTON *******************************"  << G4endl;
                    //G4cout << "******************************* new trackIDvalue" << trackID << " *******************************"  << G4endl;
                }

                G4int particleID = track->GetParentID();
                G4String particleType = track->GetDefinition()->GetParticleName();
                //G4cout << "******************************* particleID(parent) " << particleID << " *******************************                          <-------------------- PARENT"  << G4endl;

                for (int trackIDInVector : trackIDVector) {
                    //G4cout << "******************************* trackIDInVector " << trackIDInVector << " *******************************"  << G4endl;
                    if (trackIDInVector == particleID && particleType == "e-"){
                        energy_deposit_total_2 += energy_deposit_step; //with transportation and just compton
                        //G4cout << "******************************* COMPTON ELECTRON *******************************"  << G4endl;
                        savedEnergyStep = true;
                    }
                }

                if (savedEnergyStep == false) {
                    energy_deposit_total_1 += energy_deposit_step;  //with transportation and without compton
                    //G4cout << "*******************************  NO COMPTON ELECTRON *******************************"  << G4endl;
                }
            }
        }
    };

    [[maybe_unused]]
    auto delete_track = [](G4Track const* track) {
        // G4int parentID = track -> GetParentID();
        // G4int trackID  = track -> GetTrackID();
        G4ThreeVector momentum = track-> GetMomentumDirection();
        //G4cout << "******************************* Momentum: " << momentum.x() << "  " << "Parent ID: " << parentID << "  " <<"Track ID: " << trackID << " *******************************"  << G4endl;
        if (track -> GetVolume() -> GetLogicalVolume() -> GetName() != "vessel_steel") {
            G4Track* non_const_track = const_cast<G4Track*>(track);
            non_const_track->SetTrackStatus(fStopAndKill);
            //G4cout << "******************************* Momentum: " <<track->GetVolume()->GetLogicalVolume()->GetName() << " *******************************"  << G4endl;
            //auto trackStatus = track->GetTrackStatus();
            //if (trackStatus == fStopAndKill) {
            //G4cout << "******************************* DIED  *******************************"  << G4endl;
            //}
        }
    };

    [[maybe_unused]]
    auto create_trackIDVector = [& trackIDVector](G4Track const* track) {
        G4int trackID = track->GetTrackID();
        trackIDVector.push_back(trackID);
        //G4cout << "******************************* " << "NUEVO TRACK AÑADIDO " << trackID << " *******************************"  << G4endl;
    };

    auto delete_file_long = [&filename_step, &eventCounter](G4Run const*) {

        std::ofstream file(filename_step, std::ios::out);
        int colWidth = 20;
        file << std::left << std::setw(colWidth) << "PDG encoding";
        file << std::left << std::setw(colWidth) << "Particle ID";
        file << std::left << std::setw(colWidth) << "Track ID";
        file << std::left << std::setw(colWidth) << "Particle type";
        file << std::left << std::setw(colWidth) << "Time";
        file << std::left << std::setw(colWidth) << "Step lenght";
        file << std::left << std::setw(colWidth) << "Velocity";
        file << std::left << std::setw(colWidth) << "Speed";
        file << std::left << std::setw(colWidth) << "X";
        file << std::left << std::setw(colWidth) << "Y";
        file << std::left << std::setw(colWidth) << "Z";
        file << std::left << std::setw(colWidth) << "Kinetic Energy";
        file << std::left << std::setw(colWidth) << "Deposit energy";
        file << std::left << std::setw(colWidth) << "Interaction type";
        file << std::endl;
        file.close();
    
        eventCounter = 0;
    };

    auto delete_file_short = [&filename_event](G4Run const*) {
        std::ofstream file(filename_event, std::ios::out); // Implicitly closed when `file` goes out of scope
    };

    auto delete_file_short_and_long = [&delete_file_short, &delete_file_long] (auto run) {
        delete_file_short(run);
        delete_file_long(run);
    };

    [[maybe_unused]]
    auto delete_file_short_double = [&filename_event_1, &filename_event_2](G4Run const*) {
        std::ofstream file1(filename_event_1, std::ios::out); // Implicitly closed when `file1` goes out of scope
        std::ofstream file2(filename_event_2, std::ios::out); // Ditto
    };

    auto reset_energy = [&energy_deposit_total, &counts](G4Event const*){
        energy_deposit_total = 0.0;
        counts = 0.0;
    };

    [[maybe_unused]]
    auto reset_energy_and_trackIDVector = [&energy_deposit_total, &counts, &trackIDVector](G4Event const*){
        energy_deposit_total = 0.0;
        counts = 0.0;

        std::vector<int> trackIDVector_empty;
        trackIDVector = trackIDVector_empty;
    };

    [[maybe_unused]]
    auto reset_energy_double = [&energy_deposit_total_1, &energy_deposit_total_2](G4Event const*){
        energy_deposit_total_1 = 0.0;
        energy_deposit_total_2 = 0.0;
    };

    [[maybe_unused]]
    auto reset_eventCounter = [&eventCounter](G4Run const*){ eventCounter = 0; };

    n4::silence hush{G4cout};

    G4int verbosity = 0;
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list ->  ReplacePhysics(new G4EmStandardPhysics_option4());
    //physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    physics_list -> RegisterPhysics(new G4RadioactiveDecayPhysics);
    physics_list -> RegisterPhysics(new G4DecayPhysics());


    auto run_manager = std::unique_ptr<G4RunManager>
        {G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial)};

    // Physics list must be attached to run manager before instantiating other user action classes
    run_manager -> SetUserInitialization(physics_list);
    run_manager -> SetUserInitialization((new n4::actions{kr83m})
                                                -> set(new n4::stepping_action{write_info_and_get_energy_step})
                                                //-> set((new n4::tracking_action) -> post(create_trackIDVector))
                                                //-> set((new n4::tracking_action) -> pre(delete_track))
                                                //-> set(new n4::stepping_action{get_energy})
                                                //-> set(new n4::stepping_action{get_energy_and_check_track})
                                                //-> set((new n4::event_action) -> begin(reset_energy))
                                                //-> set((new n4::event_action) -> end(write_energy_event) -> begin(reset_energy_and_trackIDVector))
                                                -> set((new n4::event_action) -> end(write_energy_event) -> begin(reset_energy))
                                                -> set((new n4::run_action) -> begin(delete_file_short_and_long)));
                                                //-> set((new n4::run_action) -> end(print_energy)));
                                                //-> set((new n4::run_action) -> end(reset_eventCounter)));
    run_manager -> SetUserInitialization(new n4::geometry{geometry});
    run_manager -> Initialize();

    n4::ui(argc, argv);
}
