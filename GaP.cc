#include "nain4.hh"
#include "g4-mandatory.hh"
#include "n4_ui.hh"
#include "n4-utils.hh"
#include "n4-volumes.hh"

#include "geometry.hh"
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
#include <tuple>

auto get_pre_volume_name(G4Step const * const step) {
    return step -> GetPreStepPoint() -> GetTouchableHandle() -> GetVolume() -> GetLogicalVolume() -> GetName();
};

int main(int argc, char *argv[]) {
    std::cout << "Hello World!" << std::endl;
    
    G4double  energy_deposit_total;
    G4double  energy_deposit_total_1;
    G4double  energy_deposit_total_2;
    G4int  counts;
    G4int  eventCounter;
    std::vector<int> trackIDVector;
    
    const std::string& filename_event = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation_G4.txt";
    const std::string& filename_step = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation.txt";
    
    const std::string& filename_event_1 = "EnergyDepositTotal_Kr83m_10bar_sum_withoutTransportation.txt";
    const std::string& filename_event_2 = "EnergyDepositTotal_Kr83m_10bar_sum_justCompton.txt";
    
////////////////////////////////////////////////////////////////////////
//tracking actions

   [[maybe_unused]]
    auto create_trackIDVector = [& trackIDVector](G4Track const* track) {
        G4int trackID = track->GetTrackID();
        trackIDVector.push_back(trackID);
        //G4cout << "******************************* " << "NUEVO TRACK AÑADIDO " << trackID << " *******************************"  << G4endl;
    };

////////////////////////////////////////////////////////////////////////
//stepping actions

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
                    //banner("COMPTON");
                    //G4cout << "******************************* new trackIDvalue" << trackID << " *******************************"  << G4endl;
                }

                G4int particleID = track->GetParentID();
                G4String particleType = track->GetDefinition()->GetParticleName();
                //G4cout << "******************************* particleID(parent) " << particleID << " *******************************                          <-------------------- PARENT"  << G4endl;

                for (int trackIDInVector : trackIDVector) {
                    //G4cout << "******************************* trackIDInVector " << trackIDInVector << " *******************************"  << G4endl;
                    if (trackIDInVector == particleID && particleType == "e-"){
                        energy_deposit_total_2 += energy_deposit_step; //with transportation and just compton
                        //banner("COMPTON ELECTRON");
                        savedEnergyStep = true;
                    }
                }

                if (savedEnergyStep == false) {
                    energy_deposit_total_1 += energy_deposit_step;  //with transportation and without compton
                    //banner("NO COMPTON ELECTRON");
                }
            }
        }
    };



////////////////////////////////////////////////////////////////////////
//event actions

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

////////////////////////////////////////////////////////////////////////
//run actions

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

    [[maybe_unused]]
    auto reset_eventCounter = [&eventCounter](G4Run const*){ eventCounter = 0; };
    
    
    auto delete_file_map_and_reset_eventCounter = [& eventCounter](G4Run const* run){
		const std::string& filename_map = "Detector_map_test.txt";
		const std::string& filename_map_primaries = "Detector_map_test_primaries.txt";
		
		std::ofstream file1(filename_map, std::ios::out);
		int colWidth = 20;
		file1 << std::left << std::setw(colWidth) << "X";
		file1 << std::left << std::setw(colWidth) << "Y";
		file1 << std::left << std::setw(colWidth) << "Z";
		file1 << std::left << std::setw(colWidth) << "Hits check";
		file1 << std::endl;
		file1.close();
		
		std::ofstream file2(filename_map_primaries, std::ios::out);
		file2 << std::left << std::setw(colWidth) << "X";
		file2 << std::left << std::setw(colWidth) << "Y";
		file2 << std::left << std::setw(colWidth) << "Z";
		file2 << std::left << std::setw(colWidth) << "Hits check primaries";
		file2 << std::endl;
		file2.close();
    
		eventCounter = 0;
		
	};

////////////////////////////////////////////////////////////////////////

    n4::silence hush{G4cout};

    G4int verbosity = 0;
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list ->  ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    physics_list -> RegisterPhysics(new G4RadioactiveDecayPhysics);
    physics_list -> RegisterPhysics(new G4DecayPhysics());


    auto run_manager = std::unique_ptr<G4RunManager>
        {G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial)};

    // Physics list must be attached to run manager before instantiating other user action classes
    run_manager -> SetUserInitialization(physics_list);
    
    //G4ParticleTable needs to be call after G4VUserPhysicsList is instantiated and assigned to G4RunManager
    auto opticalphoton = [](auto event){generate_particles_in_event(event, random_generator_inside_drift({}), generate_partilces_and_energies_tuples());};   	
		//auto box_source = [](auto event){generate_particles_in_event(event, {0., 0., 167.6775*mm + 50.*mm}, generate_partilces_and_energies_tuples());};  //From the box_source
		//auto kr83m = [](auto event){generate_ion_decay(event, random_generator_inside_drift({}), 0);}; 
		//auto kr83m_nexus= [](auto event){ kr83_generator(event, 32.1473*keV, 9.396*keV,  0.0490, 154.*ns); }; 
		//auto Co57 = [vessel_out_rad_, angle, source_pos_z](auto event){generate_ion_decay(event, {vessel_out_rad_*cos(angle), vessel_out_rad_*sin(angle), source_pos_z}, 0);};  //From the surface
		//auto Am241 = [cathode_z](auto event){generate_ion_decay(event, {0., 0., cathode_z}, 0);};  //From the surface of the cathode

    
    run_manager -> SetUserInitialization((new n4::actions{opticalphoton})
                                                //-> set(new n4::stepping_action{write_info_and_get_energy_step})
                                                //-> set((new n4::tracking_action) -> post(create_trackIDVector) -> pre(delete_track)
                                                //-> set((new n4::event_action) -> end(write_energy_event) -> begin(reset_energy))
                                                -> set((new n4::run_action) -> begin(delete_file_map_and_reset_eventCounter)));
                                                
     run_manager -> SetUserInitialization(new n4::geometry{geometry});

    // auto world = get_world();
    // //auto& place_something_in = place_mesh_holder_in;
    // //auto& place_something_in = place_quartz_window_holder_in;
    // //auto& place_something_in = place_pmt_holder_in;
    // auto place_something_in = [](auto world){ place_rings_in(world, model_new()); };
    // //auto place_something_in = [](auto world){ place_anode_el_gate_in(world, model_xxx_0()); };
    // run_manager -> SetUserInitialization(new n4::geometry{[&] { place_something_in(world); return n4::place(world).now(); }});

    run_manager -> Initialize();

    n4::ui(argc, argv);
}
