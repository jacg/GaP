#include "materials.hh"

#include "nain4.hh"
#include "g4-mandatory.hh"
#include "geometry.hh"

#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4Material.hh>
#include <G4OpticalPhysics.hh>
#include <G4RandomDirection.hh>
#include <G4RunManagerFactory.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>

#include <iostream>

G4PVPlacement* geometry() {

  auto model_new_= 1;

  G4double world_size = 0.5 * m;

  G4double vessel_out_rad_    = 288./2  *mm;
  G4double vessel_out_length_ = 46.679  *cm;
  G4double vessel_rad_        = 276./2  *mm;
  G4double vessel_length_     = 38.599  *cm; // Adjusted length so that the gas volume is centered. Original length (38.639  *cm),

  G4double mesh_rad_          = 104./2  *mm;
  G4double mesh_thickn_       = 0.075   *mm;
  G4double mesh_transparency_ = 0.95;

  G4double meshBracket_rad_      = 180./2  *mm;
  G4double meshBracket_thickn_   = 6.      *mm;
  G4double anodeBracket_rad_     = 160./2  *mm;
  G4double anodeBracket_thickn_  = 6.975   *mm;

  G4double pmt_rad_ = 25.4/2  *mm;

  G4double enclosure_pmt_rad_        = 120./2  *mm;
  G4double enclosure_pmt_thickn_     = 8.5     *mm;
  G4double enclosure_pmt_length_     = 113.5   *mm;
  G4double enclosurevac_pmt_length_  = 110.5   *mm;

  G4double plate_pmt_rad_        = 105./2  *mm;
  G4double plate_pmt_thickn_     = 105./2  *mm;
  G4double plate_pmt_length_     = 10      *mm;
  G4double plateUp_pmt_length_   = 15      *mm;
  G4double plateUp_pmt_thickn_   = 21.5    *mm;

  G4double pmtHolder_rad_        = 115./2  *mm;
  G4double pmtHolder_length_     = 9.      *mm;

  G4double quartz_window_rad_    = 108./2  *mm;
  G4double quartz_window_thickn_ = 3       *mm;
  G4double tpb_coating_thickn_   = 3       *micrometer;

  G4double photoe_prob_       (0.);

  G4double pressure_          (10.* bar);
  G4double temperature_       (293. * kelvin);
  //G4double sc_yield_          (22222./MeV); // Wsc = 45 eV, fr
  G4double sc_yield_          (1./GeV);
  //G4double sc_yield_          (1000./MeV);
  G4double elifetime_         (1e6* ms);
  G4double drift_vel_         (1. * mm/microsecond);
  G4double drift_transv_diff_ (1. * mm/sqrt(cm));
  G4double drift_long_diff_   (.3 * mm/sqrt(cm));
  G4double el_field_          (16.0 * kilovolt/cm);
  G4double el_vel_            (3. * mm/microsecond);
  G4double el_transv_diff_    (1. * mm/sqrt(cm));
  G4double el_long_diff_      (.3 * mm/sqrt(cm));

  auto Cu = n4::material("G4_Cu");
  auto vacuum = n4::material("G4_Galactic");
  auto steel = n4::material("G4_STAINLESS-STEEL");

  auto gas_ = GAr_with_properties( pressure_, temperature_, sc_yield_, elifetime_);
  auto mesh_mat = FakeDielectric_with_properties(gas_, "mesh_mat",
                                                 pressure_, temperature_, mesh_transparency_, mesh_thickn_,
                                                 sc_yield_, elifetime_, photoe_prob_);
  auto peek = peek_with_properties();
  auto quartz = quartz_with_properties();
  auto tpb = TPB_with_properties();

  //Cylinder, acting as the vessel
  auto world = n4::volume<G4Box>("world", vacuum, world_size/2, world_size/2, world_size/2);

  //Cylinder, acting as the vessel
  auto vessel_steel = n4::volume<G4Tubs>("vessel_steel", steel, 0., vessel_out_rad_, vessel_out_length_/2 , 0., 360.*deg);
  n4::place(vessel_steel).in(world).at({0, 0, 0}).check_overlaps().now();

  //Build inside detector
  //BuildTPC(gas_, mesh_mat, steel, peek, vacuum, quartz, tpb, vessel_steel);


  auto vessel = n4::volume<G4Tubs>("GasVessel", gas_, 0., vessel_rad_, vessel_length_/2 , 0., 360.*deg);
  n4::place(vessel).in(vessel_steel).at({0, 0, 0}).check_overlaps().now();


  G4double cathode_z;
  G4LogicalVolume * gas_el;
  G4double drift_length_  ;
  G4double el_length_    ;

  if (model_new_ == 0) {

    G4double d=3.*mm;
    G4double D=5.*mm;
    G4double ring_rad_int_=130./2*mm;
    G4double ring_rad_out_=140./2*mm;
    G4double ring_thickn_=10.*mm;

    // Cathode
    auto cathode= n4::volume<G4Tubs>("cathode", mesh_mat, 0., mesh_rad_, (mesh_thickn_)/2, 0., 360.*deg);
    //G4double cathode_z = 4.505*mm + mesh_thickn_/2 + 2*D + 5*d + 6*ring_thickn_;  //cathode center from vessel center
    cathode_z = 90.1125*mm - 15.745*mm;  //cathode center from vessel center
    n4::place(cathode).in(vessel).at({0, 0, cathode_z}).check_overlaps().now();

    // Cathode Bracket
    auto cathBracket = n4::volume<G4Tubs>("CathodeBracket", steel, mesh_rad_, meshBracket_rad_, (meshBracket_thickn_)/2, 0., 360.*deg);
    //G4double cathBracket_z = 8.005*mm - meshBracket_thickn_/2 + 2*D + 5*d + 6*ring_thickn_;
    G4double cathBracket_z = cathode_z;
    n4::place(cathBracket).in(vessel).at({0, 0, cathBracket_z}).check_overlaps().now();

    //Cu rings
    auto ring= n4::volume<G4Tubs>("ring", Cu, ring_rad_int_, ring_rad_out_, (ring_thickn_)/2, 0., 360.*deg);
    n4::place(ring).in(vessel).at({0, 0, cathode_z - meshBracket_thickn_/2- D - ring_thickn_/2}).copy_no(0).check_overlaps().now();
    n4::place(ring).in(vessel).at({0, 0, cathode_z - meshBracket_thickn_/2- D - ring_thickn_ - d - ring_thickn_/2}).copy_no(1).check_overlaps().now();
    n4::place(ring).in(vessel).at({0, 0, cathode_z - meshBracket_thickn_/2- D - 2*ring_thickn_ - 2*d - ring_thickn_/2}).copy_no(2).check_overlaps().now();
    n4::place(ring).in(vessel).at({0, 0, cathode_z - meshBracket_thickn_/2- D - 3*ring_thickn_ - 3*d - ring_thickn_/2}).copy_no(3).check_overlaps().now();
    n4::place(ring).in(vessel).at({0, 0, cathode_z - meshBracket_thickn_/2- D - 4*ring_thickn_ - 4*d - ring_thickn_/2}).copy_no(4).check_overlaps().now();
    n4::place(ring).in(vessel).at({0, 0, cathode_z - meshBracket_thickn_/2- D - 5*ring_thickn_ - 5*d - ring_thickn_/2}).copy_no(5).check_overlaps().now();

    // Source box
    G4double source_box_width = 100.*mm;
    G4double source_box_lenght_ = 50.*mm;
    auto source_box = n4::volume<G4Box>("source_box", steel, source_box_width/2, source_box_width/2, source_box_lenght_/2);
    G4double source_box_z = cathode_z + 81.64*mm + source_box_lenght_/2;
    //G4double source_box_z = cathode_z ;
    n4::place(source_box).in(vessel).at({0., 0., source_box_z}).check_overlaps().now();

    //Gas
    drift_length_  = 96.*mm - meshBracket_thickn_ ;
    el_length_     = 15.*mm - anodeBracket_thickn_/2;

    // Drift
    auto gas_drift = n4::volume<G4Tubs>("gas_drift", gas_, 0., meshBracket_rad_, (drift_length_)/2, 0., 360.*deg);
    G4double drift_z = cathode_z - mesh_thickn_/2 - drift_length_/2;
    n4::place(gas_drift).in(vessel).at({0, 0, drift_z}).check_overlaps().now();

    // EL gap
    gas_el = n4::volume<G4Tubs>("gas_el", gas_, 0., anodeBracket_rad_, (el_length_)/2, 0., 360.*deg);
    G4double el_z = drift_z - drift_length_/2 - el_length_/2;
    n4::place(gas_el).in(vessel).at({0, 0, el_z}).check_overlaps().now();
    //el_gen_  = new CylinderPointSampler2020(el_phys_);

  }  else {

    //Cathode
    auto cathode = n4::volume<G4Tubs>("cathode", mesh_mat, 0., mesh_rad_, (mesh_thickn_)/2, 0., 360.*deg);
    cathode_z = 4.505*mm + mesh_thickn_/2;  //cathode center from vessel center
    n4::place(cathode).in(vessel).at({0, 0, cathode_z}).check_overlaps().now();

    //Cathode Bracket
    auto cathBracket = n4::volume<G4Tubs>("CathodeBracket", steel, mesh_rad_, meshBracket_rad_, (meshBracket_thickn_)/2, 0., 360.*deg);
    G4double cathBracket_z = 8.005*mm - meshBracket_thickn_/2;
    n4::place(cathBracket).in(vessel).at({0, 0, cathBracket_z}).check_overlaps().now();

    //Gas
    drift_length_  = 19.825*mm - mesh_thickn_ ;
    el_length_     = 10.775*mm + mesh_thickn_;

    // Drift
    auto gas_drift = n4::volume<G4Tubs>("gas_drift", gas_, 0., anodeBracket_rad_, (drift_length_)/2, 0., 360.*deg);
    G4double drift_z = cathode_z - mesh_thickn_/2 - drift_length_/2;
    n4::place(gas_drift).in(vessel).at({0, 0, drift_z}).check_overlaps().now();


    // EL gap
    //auto gas_el = n4::volume<G4Tubs>("gas_el", gas_, 0, mesh_rad_, (el_length_)/2, 0., 360.*deg);
    gas_el = n4::volume<G4Tubs>("gas_el", gas_, 0., anodeBracket_rad_, (el_length_)/2, 0., 360.*deg);
    G4double el_z = drift_z - drift_length_/2 - el_length_/2;
    n4::place(gas_el).in(vessel).at({0, 0, el_z}).check_overlaps().now();
    //el_gen_  = new CylinderPointSampler2020(el_phys_);
  }

  // Gate
  auto gate = n4::volume<G4Tubs>("gate", mesh_mat, 0., mesh_rad_, (mesh_thickn_)/2, 0., 360.*deg);
  G4double gate_z = el_length_/2 - mesh_thickn_/2;
  n4::place(gate).in(gas_el).at({0, 0, gate_z}).check_overlaps().now();

  // Gate Bracket
  auto gateBracket= n4::volume<G4Tubs>("gateBracket", steel, mesh_rad_, meshBracket_rad_, (meshBracket_thickn_)/2, 0., 360.*deg);
  G4double gateBracket_z = 12.745*mm + meshBracket_thickn_/2;
  n4::place(gateBracket).in(vessel).at({0, 0, -gateBracket_z}).check_overlaps().now();

  G4RotationMatrix* Rot45 = new G4RotationMatrix();
  Rot45->rotateZ(45*deg);
  G4RotationMatrix* Rot_45 = new G4RotationMatrix();
  Rot_45->rotateZ(-45*deg);
  G4RotationMatrix* Rot135 = new G4RotationMatrix();
  Rot135->rotateZ(135*deg);
  G4RotationMatrix* Rot_135 = new G4RotationMatrix();
  Rot_135->rotateZ(-135*deg);

  //Anode
  auto anode = n4::volume<G4Tubs>("Anode", mesh_mat, 0.,  mesh_rad_, (mesh_thickn_)/2, 0., 360.*deg);
  G4double anode_z = - el_length_/2 + mesh_thickn_/2;
  n4::place(anode).in(gas_el).at({0., 0., anode_z}).check_overlaps().now();

  //Anode Bracket
  auto anodeBracket = n4::volume<G4Tubs>("AnodeBracket", steel, mesh_rad_, anodeBracket_rad_, (anodeBracket_thickn_)/2, 0., 360.*deg);
  G4double anodeBracket_z = gateBracket_z + meshBracket_thickn_/2 + 3.775*mm + anodeBracket_thickn_/2;
  n4::place(anodeBracket).in(vessel).at({0., 0., -anodeBracket_z}).check_overlaps().now();

  // Peek Anode Holder
  G4double anodeHolder_length_  = 30        *mm;
  G4double anodeHolder_rad_     = 145.001/2 *mm;
  G4double anodeHolder_thickn_  = 17.       *mm;
  G4double anodeHolder_angle_   = 7.964     *deg;

  G4double anodeHolder_hole_length_ = 20    *mm;
  G4double anodeHolder_hole_rad_    = anodeHolder_rad_ + 8*mm;
  G4double anodeHolder_hole_thickn_ = 10    *mm;
  G4double anodeHolder_hole_angle_  = 7.965 *deg;


  G4Tubs *solid_anodeHolder_block = new G4Tubs("AnodeHolder", anodeHolder_rad_, anodeHolder_rad_+anodeHolder_thickn_, anodeHolder_length_/2, -anodeHolder_angle_/2, anodeHolder_angle_);
  G4Tubs *solid_anodeHolder_hole = new G4Tubs("AnodeHolderHole", anodeHolder_hole_rad_, anodeHolder_hole_rad_+anodeHolder_hole_thickn_, anodeHolder_hole_length_/2, -anodeHolder_hole_angle_/2, anodeHolder_hole_angle_);
  G4VSolid        *solidSub_anodeHolder = new G4SubtractionSolid("AnodeHolderHole_Sub", solid_anodeHolder_block, solid_anodeHolder_hole, 0, G4ThreeVector(0,0,anodeHolder_length_/2-anodeHolder_hole_length_/2) );
  G4LogicalVolume *logic_anodeHolder = new G4LogicalVolume(solidSub_anodeHolder, peek, "AnodeHolder");

  G4double anodeHolder_z = 29.495*mm + anodeHolder_length_/2;
  n4::place(logic_anodeHolder).in(vessel).rotate(*Rot45).at({0., 0., -anodeHolder_z}).copy_no(0).check_overlaps().now();
  n4::place(logic_anodeHolder).in(vessel).rotate(*Rot135).at({0., 0., -anodeHolder_z}).copy_no(1).check_overlaps().now();
  n4::place(logic_anodeHolder).in(vessel).rotate(*Rot_45).at({0., 0., -anodeHolder_z}).copy_no(2).check_overlaps().now();
  n4::place(logic_anodeHolder).in(vessel).rotate(*Rot_135).at({0., 0., -anodeHolder_z}).copy_no(3).check_overlaps().now();


  /////////////////////////

  //Quartz Window
  auto quartz_window = n4::volume<G4Tubs>("QuartzWindow", quartz, 0., quartz_window_rad_, quartz_window_thickn_/2, 0., 360*deg);
  G4double quartz_window_z = 35.495*mm + quartz_window_thickn_/2 ;
  n4::place(quartz_window).in(vessel).at({0., 0., -quartz_window_z}).check_overlaps().now();

  //Evaporated TPB
  auto tpb_coating = n4::volume<G4Tubs>("CoatingTPB", tpb, 0., quartz_window_rad_, tpb_coating_thickn_/2, 0., 360*deg);
  G4double tpb_coating_z = quartz_window_z - quartz_window_thickn_/2 - tpb_coating_thickn_/2 ;
  n4::place(tpb_coating).in(vessel).at({0., 0., -tpb_coating_z}).check_overlaps().now();

  // Peek Quartz Window Holder
  G4double windowHolder_angle_   = 10      *deg;
  G4double windowHolder_length1_ = 4.       *mm;
  G4double windowHolder_length2_ = 19.      *mm;
  G4double windowHolder_length3_ = 5.       *mm;
  G4double windowHolder_rad0_    = 100.4/2  *mm;
  G4double windowHolder_rad1_    = 118.4/2  *mm;
  G4double windowHolder_rad2_    = 122.4/2  *mm;
  G4double windowHolder_rad3_    = 140.4/2  *mm;

  G4double windowHolder_block2_z = windowHolder_length1_/2 + windowHolder_length2_/2 ;
  G4double windowHolder_block3_z = windowHolder_length3_/2 + windowHolder_length2_/2 + windowHolder_block2_z ;

  G4Tubs *solid_windowHolder_block1 = new G4Tubs("QuartzWindowHolder1", windowHolder_rad0_, windowHolder_rad2_, windowHolder_length1_/2, -windowHolder_angle_/2, windowHolder_angle_);
  G4Tubs *solid_windowHolder_block2 = new G4Tubs("QuartzWindowHolder2", windowHolder_rad1_, windowHolder_rad2_, windowHolder_length2_/2, -windowHolder_angle_/2, windowHolder_angle_);
  G4Tubs *solid_windowHolder_block3 = new G4Tubs("QuartzWindowHolder3", windowHolder_rad1_, windowHolder_rad3_, windowHolder_length3_/2, -windowHolder_angle_/2, windowHolder_angle_);
  G4VSolid        *solidUni_windowHolder0 = new G4UnionSolid("WindowHolder_Union1", solid_windowHolder_block1, solid_windowHolder_block2, 0, G4ThreeVector(0,0,-windowHolder_block2_z) );
  G4VSolid        *solidUni_windowHolder = new G4UnionSolid("WindowHolder_Union2", solidUni_windowHolder0, solid_windowHolder_block3, 0, G4ThreeVector(0,0,-windowHolder_block3_z) );
  G4LogicalVolume *logic_windowHolder = new G4LogicalVolume(solidUni_windowHolder, peek, "QuartzWindowHolder");

  G4double quartz_windowHolder_z = quartz_window_z - quartz_window_thickn_/2 -  windowHolder_length1_/2;
  n4::place(logic_windowHolder).in(vessel).rotate(*Rot45).at({0., 0., -quartz_windowHolder_z}).copy_no(0).check_overlaps().now();
  n4::place(logic_windowHolder).in(vessel).rotate(*Rot135).at({0., 0., -quartz_windowHolder_z}).copy_no(1).check_overlaps().now();
  n4::place(logic_windowHolder).in(vessel).rotate(*Rot_45).at({0., 0., -quartz_windowHolder_z}).copy_no(2).check_overlaps().now();
  n4::place(logic_windowHolder).in(vessel).rotate(*Rot_135).at({0., 0., -quartz_windowHolder_z}).copy_no(3).check_overlaps().now();


  G4double windowHolderTop_length1_ = 3.  *mm;
  G4double windowHolderTop_length2_ = 4.  *mm;
  G4double windowHolderTop_rad0_    = windowHolder_rad0_ ;
  G4double windowHolderTop_rad1_    = windowHolderTop_rad0_ + 4.*mm ;
  G4double windowHolderTop_rad2_    = windowHolder_rad1_ ;

  G4double windowHolderTop_block2_z = windowHolderTop_length1_/2 + windowHolderTop_length2_/2 ;

  G4Tubs *solid_windowHolderTop_block1 = new G4Tubs("QuartzWindowHolderTop1", windowHolderTop_rad1_, windowHolderTop_rad2_, windowHolderTop_length1_/2, -windowHolder_angle_/2, windowHolder_angle_);
  G4Tubs *solid_windowHolderTop_block2 = new G4Tubs("QuartzWindowHolderTop2", windowHolderTop_rad0_, windowHolderTop_rad2_, windowHolderTop_length2_/2, -windowHolder_angle_/2, windowHolder_angle_);
  G4VSolid        *solidUni_windowHolderTop = new G4UnionSolid("WindowHolderTop_Union", solid_windowHolderTop_block1, solid_windowHolderTop_block2, 0, G4ThreeVector(0,0,-windowHolderTop_block2_z) );
  G4LogicalVolume *logic_windowHolderTop = new G4LogicalVolume(solidUni_windowHolderTop, peek, "QuartzWindowHolderTop");

  G4double quartz_windowHolderTop_z = quartz_window_z - quartz_window_thickn_/2 +  windowHolderTop_length1_/2;
  n4::place(logic_windowHolderTop).in(vessel).rotate(*Rot45).at({0., 0., -quartz_windowHolderTop_z}).copy_no(0).check_overlaps().now();
  n4::place(logic_windowHolderTop).in(vessel).rotate(*Rot135).at({0., 0., -quartz_windowHolderTop_z}).copy_no(2).check_overlaps().now();
  n4::place(logic_windowHolderTop).in(vessel).rotate(*Rot_45).at({0., 0., -quartz_windowHolderTop_z}).copy_no(3).check_overlaps().now();
  n4::place(logic_windowHolderTop).in(vessel).rotate(*Rot_135).at({0., 0., -quartz_windowHolderTop_z}).copy_no(4).check_overlaps().now();

  //Build PMT
  //G4double pmt_length_ = pmt_.Length();
  G4double pmt_length_ = 43.0 * mm;
  G4double pmt_z  = 42.495*mm + pmt_length_/2;

  G4Tubs *solid_pmt = new G4Tubs("SolidPMT", 0., pmt_rad_, pmt_length_/2, 0., 360.*deg); // Hamamatsu pmt length: 43*mm | STEP pmt gap length: 57.5*mm

  // Position pairs (x,Y) for PMTs
  std::vector <float> pmt_PsX={-15.573, 20.68, -36.253, 0., 36.253, -20.68, 15.573};
  std::vector <float> pmt_PsY={-32.871, -29.922, -2.949, 0., 2.949, 29.922, 32.871};

  // PMT clad
  G4VSolid *solid_enclosure_pmt = new G4Tubs("EnclosurePMT", 0, enclosure_pmt_rad_+ enclosure_pmt_thickn_ , (enclosure_pmt_length_)/2, 0., 360.*deg);
  G4double enclosure_pmt_z = vessel_length_/2 - enclosure_pmt_length_/2;
  G4double relative_pmt_z  = enclosure_pmt_z - pmt_z;
  // Vacuum inside the pmt enclosure
  G4double enclosurevac_pmt_z = vessel_length_/2 - enclosurevac_pmt_length_/2;
  G4double relativevac_pmt_z  = enclosurevac_pmt_z - pmt_z;
  G4VSolid *solid_enclosurevac_pmt = new G4Tubs("EnclosureVacPMT", 0, enclosure_pmt_rad_, (enclosurevac_pmt_length_)/2, 0., 360.*deg);
  // PMT Holder
  G4double pmtHolder_z = enclosurevac_pmt_length_/2 - pmtHolder_length_/2;
  G4VSolid *solid_pmtHolder = new G4Tubs("PMTHolder", 0, pmtHolder_rad_, (pmtHolder_length_)/2, 0., 360.*deg);
  // Steel plate enclosing the pmt tube
  G4double plate1_pmt_z = enclosure_pmt_z - enclosure_pmt_length_/2 - plate_pmt_length_/2;
  G4VSolid *solid_plate1_pmt = new G4Tubs("PMTplateBottom1", 0, plate_pmt_rad_+plate_pmt_thickn_, plate_pmt_length_/2, 0, 360*deg);

  G4ThreeVector pos_pmt = G4ThreeVector(0, 0, 0);
  G4ThreeVector pos_enclosure_pmt = G4ThreeVector(0, 0, 0);
  G4ThreeVector pos_enclosurevac_pmt = G4ThreeVector(0, 0, 0);
  G4ThreeVector pos = G4ThreeVector(0, 0, 0);

  for (G4int i = 0; i < G4int(pmt_PsX.size()); i++) {
    pos_pmt = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, -pmt_z);
    pos_enclosure_pmt = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, relative_pmt_z);
    pos_enclosurevac_pmt = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, relativevac_pmt_z);
    pos = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, 0.);

    G4VSolid*solid_enclosure_pmt = new G4SubtractionSolid("EnclosurePMT_Sub", solid_enclosure_pmt, solid_pmt,  0, pos_enclosure_pmt);
    G4VSolid*solid_enclosurevac_pmt = new G4SubtractionSolid("EnclosureVacPMT_Sub", solid_enclosurevac_pmt, solid_pmt,  0, pos_enclosurevac_pmt);
    G4VSolid*solid_pmtHolder = new G4SubtractionSolid("PMTHolder_Sub", solid_pmtHolder, solid_pmt,  0, pos);
    G4VSolid*solid_plate1_pmt = new G4SubtractionSolid("PMTplateBottom1_Sub", solid_plate1_pmt, solid_pmt,  0, pos);
  }

  G4LogicalVolume *logic_enclosure_pmt = new G4LogicalVolume(solid_enclosure_pmt, steel, "EnclosurePMT");
  n4::place(logic_enclosure_pmt).in(vessel).at({0., 0., -enclosure_pmt_z}).check_overlaps().now();

  G4LogicalVolume *logic_enclosurevac_pmt = new G4LogicalVolume(solid_enclosurevac_pmt, vacuum, "EnclosureVacPMT");
  n4::place(logic_enclosurevac_pmt).in(logic_enclosure_pmt).at({0., 0., enclosure_pmt_z - enclosurevac_pmt_z}).check_overlaps().now();

  G4LogicalVolume *logic_pmtHolder = new G4LogicalVolume(solid_pmtHolder, steel, "PMTHolder");
  n4::place(logic_pmtHolder).in(logic_enclosurevac_pmt).at({0., 0, pmtHolder_z}).check_overlaps().now();

  G4LogicalVolume *logic_plate1_pmt = new G4LogicalVolume(solid_plate1_pmt, steel, "PMTplateBottom1");
  n4::place(logic_plate1_pmt).in(vessel).at({0., 0., -plate1_pmt_z}).check_overlaps().now();

  // Steel plate attached where the peek holders are attached
  auto plate0_pmt = n4::volume<G4Tubs>("PMTplateBottom0", steel, plate_pmt_rad_, plate_pmt_rad_+plate_pmt_thickn_, plate_pmt_length_/2, 0., 360*deg);
  G4double plate0_pmt_z = plate1_pmt_z - plate_pmt_length_;
  n4::place(plate0_pmt).in(vessel).at({0., 0., -plate0_pmt_z}).check_overlaps().now();

  // Upper steel plate at the pmt clad
  G4double plateUp_pmt_rad_ = enclosure_pmt_rad_ + enclosure_pmt_thickn_;
  auto plateUp_pmt = n4::volume<G4Tubs>("PMTplateUp", steel, plateUp_pmt_rad_, plateUp_pmt_rad_+plateUp_pmt_thickn_, plateUp_pmt_length_/2, 0., 360*deg);
  G4double plateUp_pmt_z = vessel_length_/2 - plateUp_pmt_length_/2 ;
  n4::place(plateUp_pmt).in(vessel).at({0., 0., -plateUp_pmt_z}).check_overlaps().now();

  if (model_new_ == 1) {

    // A
    //G4double A_length_ = 56. *mm;
    //G4double A_width_ = 27. *mm;
    //G4double A_thickn_ = 21. *mm;

    //G4double A_hole_length_ = 50. *mm;
    //G4double A_hole_width_ = 16. *mm;

    //G4Box *solid_A_block = new G4Box("ABlock", A_width_/2, A_thickn_/2 , A_length_/2);
    //G4Box *solid_A_hole = new G4Box("AHole", A_hole_width_/2, A_thickn_/2 , A_hole_length_/2);
    //G4VSolid        *solid_A= new G4SubtractionSolid("SolidA", solid_A_block, solid_A_hole, 0, G4ThreeVector(-(A_width_ - A_hole_width_)/2,0.,-(A_length_ - A_hole_length_)/2) );
    //G4LogicalVolume *logic_A = new G4LogicalVolume(solid_A, peek, "A");

    //G4double A_xy = plate_pmt_rad_+plate_pmt_thickn_ - A_length_/2;
    //G4double A_z = plate0_pmt_z + plate_pmt_length_/2 - A_length_/2;

    //n4::place(logic_A).in(vessel).rotate(*Rot_45).at({-A_xy, -A_xy, -A_z}).copy_no(0).check_overlaps().now();
    //n4::place(logic_A).in(vessel).rotate(*Rot_135).at({-A_xy, A_xy, -A_z }).copy_no(1).check_overlaps().now();
    //n4::place(logic_A).in(vessel).rotate(*Rot45).at({A_xy, -A_xy, -A_z}).copy_no(2).check_overlaps().now();
    //n4::place(logic_A).in(vessel).rotate(*Rot135).at({A_xy, A_xy, -A_z}).copy_no(3).check_overlaps().now();
  } else {

    // Peek Mesh Holder (holds cathode and gate)
    G4double meshHolder_length_    = 36.75     *mm;
    G4double meshHolder_width_     = 21.035    *mm;
    G4double meshHolder_thickn_    = 24        *mm;

    G4double meshHolder_hole_length_   = 15.75  *mm;
    G4double meshHolder_holeUp_length_ = 6      *mm;
    G4double meshHolder_hole_thickn_   = 11.667 *mm;

    G4Box *solid_meshHolder_block = new G4Box("MeshHolder", meshHolder_width_/2, meshHolder_thickn_/2 , meshHolder_length_/2);
    G4Box *solid_meshHolder_hole = new G4Box("MeshHolderHoleDown", meshHolder_width_, meshHolder_hole_thickn_ , meshHolder_hole_length_/2);
    G4Box *solid_meshHolder_holeUp = new G4Box("MeshHolderHoleUp", meshHolder_width_, meshHolder_hole_thickn_ , meshHolder_holeUp_length_/2);
    G4VSolid        *solidSub_meshHolder0 = new G4SubtractionSolid("MeshHolderHole_SubDown", solid_meshHolder_block, solid_meshHolder_hole, 0, G4ThreeVector(0., meshHolder_thickn_/2, meshHolder_length_/2-meshHolder_hole_length_/2-5*mm) );
    G4VSolid        *solidSub_meshHolder = new G4SubtractionSolid("MeshHolderHole_SubUp", solidSub_meshHolder0, solid_meshHolder_holeUp, 0, G4ThreeVector(0., meshHolder_thickn_/2, -(meshHolder_length_/2-meshHolder_holeUp_length_/2-5*mm)) );
    G4LogicalVolume *logic_meshHolder = new G4LogicalVolume(solidSub_meshHolder, peek, "MeshHolder");

    G4double meshHolder_x = 5*mm + meshBracket_rad_ * cos(45*deg);
    G4double meshHolder_y = 5*mm + meshBracket_rad_ * sin(45*deg);
    G4double meshHolder_z = -13.005*mm + meshHolder_length_/2 ;

    n4::place(logic_meshHolder).in(vessel).rotate(*Rot_45).at({-meshHolder_x, -meshHolder_y, -meshHolder_z}).copy_no(0).check_overlaps().now();
    n4::place(logic_meshHolder).in(vessel).rotate(*Rot_135).at({-meshHolder_x, meshHolder_y, -meshHolder_z}).copy_no(1).check_overlaps().now();
    n4::place(logic_meshHolder).in(vessel).rotate(*Rot45).at({meshHolder_x, -meshHolder_y, -meshHolder_z}).copy_no(2).check_overlaps().now();
    n4::place(logic_meshHolder).in(vessel).rotate(*Rot135).at({meshHolder_x, meshHolder_y, -meshHolder_z}).copy_no(3).check_overlaps().now();

    //Steel Bar joining Mesh holder and PMT clad
    G4double meshHolderBar_rad_     = 9./2   *mm;
    G4double meshHolderBar_length_  = 35.75  *mm;

    auto meshHolderBar = n4::volume<G4Tubs>("MeshHolderBar", steel, 0., meshHolderBar_rad_, (meshHolderBar_length_)/2, 0., 360.*deg);
    G4double meshHolderBar_xy = 73.769*mm;
    G4double meshHolderBar_z  = meshHolder_z + meshHolder_length_/2 + meshHolderBar_length_/2 ;
    n4::place(meshHolderBar).in(vessel).at({meshHolderBar_xy, meshHolderBar_xy, -meshHolderBar_z}).copy_no(1).check_overlaps().now();
    n4::place(meshHolderBar).in(vessel).at({-meshHolderBar_xy, meshHolderBar_xy, -meshHolderBar_z}).copy_no(2).check_overlaps().now();
    n4::place(meshHolderBar).in(vessel).at({meshHolderBar_xy, -meshHolderBar_xy, -meshHolderBar_z}).copy_no(3).check_overlaps().now();
    n4::place(meshHolderBar).in(vessel).at({-meshHolderBar_xy, -meshHolderBar_xy, -meshHolderBar_z}).copy_no(4).check_overlaps().now();
  }
  return n4::place(world).now();
}
