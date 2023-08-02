#include "nain4.hh"
#include "g4-mandatory.hh"
#include "n4_ui.hh"

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
    auto dir = G4RandomDirection() * energy / CLHEP::c_light;
    vertex -> SetPrimary(new G4PrimaryParticle(particle, dir.x(), dir.y(), dir.z()));
};

G4ThreeVector position() {
    using namespace find_me_a_good_name;
    auto r     = G4RandFlat::shoot( 0., meshBracket_rad_);
    auto angle = G4RandFlat::shoot( 0., 2*M_PI);
    auto z     = G4RandFlat::shoot(-drift_length_/2 + drift_z, drift_length_/2 + drift_z);
    auto pos_x = r * cos(angle);
    auto pos_y = r * sin(angle);
    auto pos_z = z;
    return G4ThreeVector{pos_x, pos_y, pos_z};
};

std::unique_ptr<G4PrimaryVertex> generate_vertex() {
    auto vertex = std::make_unique<G4PrimaryVertex>();
    auto pos = position();
    vertex -> SetPosition(pos.x(), pos.y(), pos.z());
    return vertex;
};

void generate_particles_in_event(G4Event* event, std::vector<std::tuple<G4ParticleDefinition*, G4double>> particles_and_energies) {
    auto vertex = generate_vertex();
    for (auto [particle, energy] : particles_and_energies) {
        add_particle_to_vertex(vertex.get(), particle, energy);
    }
    event -> AddPrimaryVertex(vertex.release());
}

void generate_Kr83m2_decay(G4Event* event, G4double /*time*/){

    //Decay 1
    std::vector<double> probabilities_m2 = {0.76, 0.09, 0.15};
    std::random_device rd_1;
    std::mt19937 gen_1(rd_1());
    std::discrete_distribution<int> dis_1(probabilities_m2.begin(), probabilities_m2.end());

    //Decay 2 //falta tener en cuenta el tiempo medio de las partículas
    std::vector<double> probabilities_m1 = {0.95, 0.05};
    std::random_device rd_2;
    std::mt19937 gen_2(rd_2());
    std::discrete_distribution<int> dis_2(probabilities_m1.begin(), probabilities_m1.end());

    int random_event_1 = dis_1(gen_1);
    int random_event_2 = dis_2(gen_2);

    auto electron       = n4::find_particle("e-");
    auto optical_photon = n4::find_particle("opticalphoton");
    auto gamma          = n4::find_particle("gamma");

    if (random_event_1 == 0) {
        G4cout << "********************************** DECAY 1 = 0.76 **********************************" << G4endl;
        generate_particles_in_event(event, {{electron, 30 * keV},
                                            {electron,  2 * keV}});
    } else if (random_event_1 == 1) {
        G4cout << "********************************** DECAY 1 = 0.09 **********************************" << G4endl;
        generate_particles_in_event(event, {{electron, 18 * keV},
                                            {electron, 10 * keV},
                                            {electron,  2 * keV},
                                            {electron,  2 * keV}});
    } else if (random_event_1 == 2) {
        G4cout << "********************************** DECAY 1 = 0.15 **********************************" << G4endl;
        generate_particles_in_event(event, {{ electron      , 18 * keV},
                                            { optical_photon, 12 * keV},
                                            { electron      ,  2 * keV}});
    }

    if (random_event_2 == 0) {
        G4cout << "********************************** DECAY 2 = 0.95 **********************************" << G4endl;
        generate_particles_in_event(event, {{electron, 7.6 * keV},
                                            {electron, 1.8 * keV}});
    } else if (random_event_2 == 1) {
        G4cout << "********************************** DECAY 2 = 0.05 **********************************" << G4endl;
        auto vertex = generate_vertex();
        add_particle_to_vertex(vertex.get(), gamma      , 9.4 * keV);
        event -> AddPrimaryVertex(vertex.release());
        }
}


void generate_Co57(G4Event* event, G4ThreeVector position, G4double time){

    std::random_device rd;
    std::mt19937 gen(rd());
    
    G4double lifetime = 271.8 * year;
    //G4cout << "********************************** " << time << " **********************************" << G4endl;
    
    //Decay 1
    std::vector<double> probabilities_1 = {0.86, 1. - 0.86};
    std::discrete_distribution<int> dis_1(probabilities_1.begin(), probabilities_1.end());
    int random_event_1 = dis_1(gen);
    
    //Decay 2
    std::vector<double> probabilities_2 = {0.11, 1. - 0.11};
    std::discrete_distribution<int> dis_2(probabilities_2.begin(), probabilities_2.end());
    int random_event_2 = dis_2(gen);
    
    if (random_event_1 == 0) {
        //G4cout << "********************************** DECAY 1 **********************************" << G4endl;
        auto particle_1 = nain4::find_particle("gamma");
        G4ParticleGun* particleGun_1 = new G4ParticleGun(1);
        particleGun_1 -> SetParticleDefinition(particle_1);
        particleGun_1 -> SetParticleEnergy(122*keV);
        particleGun_1 -> SetParticleMomentumDirection(G4RandomDirection());
        particleGun_1 -> SetParticlePosition(position);
        particleGun_1 -> GeneratePrimaryVertex(event);

    } if (random_event_2 == 0) {
        //G4cout << "********************************** DECAY 2 **********************************" << G4endl;
        auto particle_2 = nain4::find_particle("gamma");
        G4ParticleGun* particleGun_2 = new G4ParticleGun(2);
        particleGun_2 -> SetParticleDefinition(particle_2);
        particleGun_2 -> SetParticleEnergy(136*keV);
        particleGun_2 -> SetParticleMomentumDirection(G4RandomDirection());
        particleGun_2 -> SetParticlePosition(position);
        particleGun_2 -> GeneratePrimaryVertex(event);
    }
}


void generate_Ba133(G4Event* event, G4ThreeVector position, G4double time){

    std::random_device rd;
    std::mt19937 gen(rd());
    
    //G4double lifetime = XXX * year;
    //G4cout << "********************************** " << time << " **********************************" << G4endl;
    
    //Decay 1
    std::vector<double> probabilities_1 = {0.99, 0.01};
    std::discrete_distribution<int> dis_1(probabilities_1.begin(), probabilities_1.end());
    int random_event_1 = dis_1(gen);
    
    //Decay 2
    std::vector<double> probabilities_2 = {0.62, 0.38};
    std::discrete_distribution<int> dis_2(probabilities_2.begin(), probabilities_2.end());
    int random_event_2 = dis_2(gen);
    
    //Decay 3
    std::vector<double> probabilities_3 = {0.34,0.66};
    std::discrete_distribution<int> dis_3(probabilities_3.begin(), probabilities_3.end());
    int random_event_3 = dis_3(gen);    
    
    //Decay 4
    std::vector<double> probabilities_4 = {0.18,0.82};
    std::discrete_distribution<int> dis_4(probabilities_4.begin(), probabilities_4.end());
    int random_event_4 = dis_4(gen);  
    
    if (random_event_1 == 0) {
        //G4cout << "********************************** DECAY 1 **********************************" << G4endl;
        auto particle_1 = nain4::find_particle("gamma");
        G4ParticleGun* particleGun_1 = new G4ParticleGun(1);
        particleGun_1->SetParticleDefinition(particle_1);
        auto T_1 = 31.*keV;
        particleGun_1 -> SetParticleEnergy(T_1);
        auto p_1 = G4RandomDirection();
        particleGun_1-> SetParticleMomentumDirection(p_1);
        particleGun_1->SetParticlePosition(position);
        particleGun_1->GeneratePrimaryVertex(event);

    }
    if (random_event_2 == 0) {
        //G4cout << "********************************** DECAY 2 **********************************" << G4endl;
        auto particle_2 = nain4::find_particle("gamma");
        G4ParticleGun* particleGun_2 = new G4ParticleGun(2);
        particleGun_2->SetParticleDefinition(particle_2);
        auto T_2 = 356.*keV;
        particleGun_2 -> SetParticleEnergy(T_2);
        auto p_2 = G4RandomDirection();
        particleGun_2-> SetParticleMomentumDirection(p_2);
        particleGun_2->SetParticlePosition(position);
        particleGun_2->GeneratePrimaryVertex(event);
    }

    if (random_event_3 == 0) {
        //G4cout << "********************************** DECAY 3 **********************************" << G4endl;
        auto particle_3 = nain4::find_particle("gamma");
        G4ParticleGun* particleGun_3 = new G4ParticleGun(3);
        particleGun_3 -> SetParticleDefinition(particle_3);
        particleGun_3 -> SetParticleEnergy(81*keV);
        particleGun_3 -> SetParticleMomentumDirection(G4RandomDirection());
        particleGun_3 -> SetParticlePosition(position);
        particleGun_3 -> GeneratePrimaryVertex(event);
    }

    if (random_event_4 == 0) {
        //G4cout << "********************************** DECAY 3 **********************************" << G4endl;
        auto particle_4 = nain4::find_particle("gamma");
        G4ParticleGun* particleGun_4 = new G4ParticleGun(4);
        particleGun_4 -> SetParticleDefinition(particle_4);
        particleGun_4 -> SetParticleEnergy(303*keV);
        particleGun_4 -> SetParticleMomentumDirection(G4RandomDirection());
        particleGun_4 -> SetParticlePosition(position);
        particleGun_4 -> GeneratePrimaryVertex(event);
    }
}
   
void generate_ion_decay(G4Event* event, G4ThreeVector position, G4double time){

    G4ParticleGun* particleGun = new G4ParticleGun(1);

    std::string IonName{"Kr83m"};

    G4int A;
    G4int Z;
    G4double E = 0.;
    G4double charge = 0. * eplus;
    auto T= 0.*keV;
    auto p = G4RandomDirection();

    if      (IonName == "Co57" ) { Z = 27; A =  57; }
    else if (IonName == "Ba133") { Z = 56; A = 133; }
    else if (IonName == "Am241") { Z = 95; A = 241; }
    else if (IonName == "Fe55" ) { Z = 26; A =  55; }
    else if (IonName == "Kr83m") { Z = 36; A =  83; E = 41.557 * keV; }

    G4ParticleDefinition* ion = G4IonTable:: GetIonTable() -> GetIon(Z, A, E);

    particleGun -> SetParticleDefinition(ion);
    particleGun -> SetParticlePosition(position);
    particleGun -> SetParticleEnergy(T);
    particleGun -> SetParticleMomentumDirection(p);
    particleGun -> SetParticleCharge(charge);

    particleGun->GeneratePrimaryVertex(event);

    delete particleGun;
}



int main(int argc, char *argv[]) {
    std::cout << "Hello World!" << std::endl;
    
    //G4double time;
    //auto get_time[](auto step){ G4double time = step -> GetDeltaTime();};
    
    G4double vessel_out_rad_    = 288./2  *mm;
    G4double vessel_out_length_ = 46.679  *cm;
    G4double angle =  0 * rad;
    G4double source_pos_z = 0 * cm; //range: [-vessel_out_length_/2,vessel_out_length_]
    G4double cathode_z = 4.505*mm; //cathode surface, not center

    auto kr83m = [cathode_z](auto event){generate_ion_decay(event, random_generator_inside_drift(), 0);};
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
    auto get_energy_old = [&energy_deposit_total](G4Step const* step){
    
    if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "gas_drift")
   {
    auto energy_deposit_step = step->GetTotalEnergyDeposit();
    energy_deposit_total += energy_deposit_step;
   }
  
    //G4cout << "********************************** " << energy_deposit_total << " **********************************" << G4endl;

    } ;
    
    auto print_energy = [& energy_deposit_total](G4Event const* event){G4cout << "*******************************/// " << energy_deposit_total << " ///*******************************"  << G4endl;} ;

    const std::string& filename_event = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation_G4.txt";
    const std::string& filename_step = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation.txt";

    //const std::string& filename_event_1 = "EnergyDepositTotal_Co57_10bar_sum_withoutCompton.txt";
    const std::string& filename_event_1 = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation.txt";
    const std::string& filename_event_2 = "EnergyDepositTotal_Kr83m_10bar_sum_justCompton.txt";

    auto write_energy_event = [& energy_deposit_total, & filename_event, & eventCounter](G4Event const* event){

    if  (energy_deposit_total != 0.0){
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

    auto write_energy_event_double = [& energy_deposit_total_1, & energy_deposit_total_2, & filename_event_1, & filename_event_2, & eventCounter](G4Event const* event){

    if  (energy_deposit_total_1 != 0.0){
    std::ofstream file;
    file.open (filename_event_1, std::ios::app);
    file << energy_deposit_total_1 * 1000. << "\n";
    //file << energy_deposit_total_1  << "\n";
    file.close();
    }
    if  (energy_deposit_total_2 != 0.0){
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
    } ;

auto write_info_and_get_energy_step = [&filename_step, &energy_deposit_total, &counts](G4Step const* step) {

    G4Track* track = step->GetTrack();
    
    if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "gas_drift") {
        auto energy_deposit_step_pre = step->GetPreStepPoint()->GetTotalEnergy();
        auto energy_deposit_step_post = step->GetPostStepPoint()->GetTotalEnergy();    
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

    auto get_energy_double = [&energy_deposit_total_1, &energy_deposit_total_2](G4Step const* step){
        if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "gas_drift"){
            auto energy_deposit_step = step->GetTotalEnergyDeposit();
            const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
            G4String interactionType = process->GetProcessName();


            if(interactionType != "Transportation"){
                energy_deposit_total_1 += energy_deposit_step;
            }
        }
    };

    auto get_energy = [&energy_deposit_total](G4Step const* step){
        if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "gas_drift"){
            auto energy_deposit_step = step->GetTotalEnergyDeposit();
            const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
            G4String interactionType = process->GetProcessName();

            if(interactionType != "Transportation"){
                energy_deposit_total += energy_deposit_step;

            }
        }
    };

    auto  get_energy_and_check_track= [&energy_deposit_total, & trackIDVector](G4Step const* step){
        if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "gas_drift"){
            auto energy_deposit_step = step->GetTotalEnergyDeposit();
            const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
            G4String interactionType = process->GetProcessName();

            G4Track* track = step->GetTrack();
            G4int trackID = track->GetTrackID();
            auto trackID_check = std::find(trackIDVector.begin(), trackIDVector.end(), trackID);

            if (trackID_check != trackIDVector.end()) {
                //G4cout << "******************************* " << "El número " << trackID << " se encuentra en el vector." << "*******************************"  << G4endl;
            }else {
                //G4cout << "******************************* " << "El número " << trackID << " NO se encuentra en el vector (INFO GUARDADA)" << "*******************************"  << G4endl;
                if(interactionType != "Transportation"){
                    energy_deposit_total += energy_deposit_step;
                }
            }
        }
    };

    auto get_energy_compton_double = [&energy_deposit_total_1, &energy_deposit_total_2, & trackIDVector](G4Step const* step){
    
        if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "gas_drift"){
            auto energy_deposit_step = step->GetTotalEnergyDeposit();
            const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
            G4String interactionType = process->GetProcessName();
            bool savedEnergyStep = false;

            if(interactionType != "Transportation"){

                //G4cout << "-----------------------------------------------------------------------------------------------------------------------------"  << G4endl;

                G4Track* track = step->GetTrack();

                if (interactionType == "compt"){
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

                if (savedEnergyStep == false){
                    energy_deposit_total_1 += energy_deposit_step;  //with transportation and without compton
                    //G4cout << "*******************************  NO COMPTON ELECTRON *******************************"  << G4endl;
                }
            }
        }
        };

    auto delete_track = [](G4Track const* track) {
        G4int parentID = track->GetParentID();
        G4int trackID = track->GetTrackID();
        G4ThreeVector momentum = track-> GetMomentumDirection();
        //G4cout << "******************************* Momentum: " << momentum.x() << "  " << "Parent ID: " << parentID << "  " <<"Track ID: " << trackID << " *******************************"  << G4endl;
        if (track->GetVolume()->GetLogicalVolume()->GetName() != "vessel_steel") {
            G4Track* non_const_track = const_cast<G4Track*>(track);
            non_const_track->SetTrackStatus(fStopAndKill);
                //G4cout << "******************************* Momentum: " <<track->GetVolume()->GetLogicalVolume()->GetName() << " *******************************"  << G4endl;
                //auto trackStatus = track->GetTrackStatus();
                //if (trackStatus == fStopAndKill) {
                //G4cout << "******************************* DIED  *******************************"  << G4endl;
                //}
    }
};

    auto create_trackIDVector = [& trackIDVector](G4Track const* track) {
        G4int trackID = track->GetTrackID();
        trackIDVector.push_back(trackID);
        //G4cout << "******************************* " << "NUEVO TRACK AÑADIDO " << trackID << " *******************************"  << G4endl;
    };

    auto delete_file_long = [&filename_step, & eventCounter](G4Run const* run){

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

    auto delete_file_short = [& filename_event](G4Run const* run){

    std::ofstream file(filename_event, std::ios::out);
    file.close();
    };

    auto delete_file_short_and_long = [&delete_file_short, &delete_file_long] (auto run) {
        delete_file_short(run);
        delete_file_long(run);
    };

    auto delete_file_short_double = [& filename_event_1, & filename_event_2](G4Run const* run){

    std::ofstream file1(filename_event_1, std::ios::out);
    file1.close();
    std::ofstream file2(filename_event_2, std::ios::out);
    file2.close();
    };

    auto reset_energy = [& energy_deposit_total, & counts](G4Event const* event){
        energy_deposit_total = 0.0;
        counts = 0.0;
    } ;

    auto reset_energy_and_trackIDVector = [& energy_deposit_total, & counts, & trackIDVector](G4Event const* event){
        energy_deposit_total = 0.0;
        counts = 0.0;

        std::vector<int> trackIDVector_empty;
        trackIDVector = trackIDVector_empty;


    } ;


    auto reset_energy_double = [& energy_deposit_total_1, & energy_deposit_total_2](G4Event const* event){
        energy_deposit_total_1 = 0.0;
        energy_deposit_total_2 = 0.0;
    } ;

    auto reset_eventCounter = [& eventCounter](G4Run const* run){
        eventCounter = 0;
    };

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
