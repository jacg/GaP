#include "nain4.hh"
#include "g4-mandatory.hh"

//#include <G4SystemOfUnits.hh>

#include "detector.hh"

bool process_hits(G4Step *step){			
	const std::string& filename_map = "Detector_map_test.txt";
	const std::string& filename_map_primaries = "Detector_map_test_primaries.txt";
	G4Track* track = step -> GetTrack();
	G4String particleType = track->GetDefinition()->GetParticleName();			
	G4int hits_check = 0;
	
	G4ThreeVector position  = track -> GetVertexPosition();
	
	G4cout << "*************************************  :)  " << particleType << "  (:  *************************************"  << G4endl;   
	
	if (particleType == "e-"){ //opticalphoton (not now because it doesn't work, abs)
	track -> SetTrackStatus(fStopAndKill);
	
	hits_check = 1;
	G4cout << "*************************************  :)  OUCH  (:  *************************************"  << G4endl; 
				
	}
	
	std::ofstream file1(filename_map, std::ios::app);
	int colWidth = 20;
	file1 << std::left << std::setw(colWidth) << position.x();
	file1 << std::left << std::setw(colWidth) << position.y();
	file1 << std::left << std::setw(colWidth) << position.z();
	file1 << std::left << std::setw(colWidth) << hits_check;
	file1 << std::endl;
	file1.close();
	
	if (track -> GetTrackID() == 1){
		std::ofstream file2(filename_map_primaries, std::ios::app);
		file2 << std::left << std::setw(colWidth) << position.x();
		file2 << std::left << std::setw(colWidth) << position.y();
		file2 << std::left << std::setw(colWidth) << position.z();
		file2 << std::left << std::setw(colWidth) << hits_check;
		file2 << std::endl;
		file2.close();
	}
	
	return true;
}
