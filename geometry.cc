#include "materials.hh"

#include "nain4.hh"
#include "g4-mandatory.hh"
#include "geometry.hh"
#include "n4-volumes.hh"
#include "n4-utils.hh"

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

  G4double photoe_prob_       = 0.;
  G4double pressure_          = 10.* bar;
  G4double temperature_       = 293. * kelvin;
  //G4double sc_yield_        =  22222./MeV; // Wsc = 45 eV, fr
  G4double sc_yield_          = 1./GeV;
  //G4double sc_yield_        =  1000./MeV;
  G4double elifetime_         = 1e6* ms;
  // G4double drift_vel_         = 1. * mm/microsecond;
  // G4double drift_transv_diff_ = 1. * mm/sqrt(cm);
  // G4double drift_long_diff_   = .3 * mm/sqrt(cm);
  // G4double el_field_          = 16.0 * kilovolt/cm;
  // G4double el_vel_            = 3. * mm/microsecond;
  // G4double el_transv_diff_    = 1. * mm/sqrt(cm);
  // G4double el_long_diff_      = .3 * mm/sqrt(cm);

  auto Cu     = n4::material("G4_Cu");
  auto vacuum = n4::material("G4_Galactic");
  auto steel  = n4::material("G4_STAINLESS-STEEL");

  auto gas_    = GAr_with_properties( pressure_, temperature_, sc_yield_, elifetime_);
  auto mesh_mat = FakeDielectric_with_properties(gas_, "mesh_mat",
                                                 pressure_, temperature_, mesh_transparency_, mesh_thickn_,
                                                 sc_yield_, elifetime_, photoe_prob_);
  auto peek   = peek_with_properties();
  auto quartz = quartz_with_properties();
  auto tpb    = TPB_with_properties();

  auto world = n4::box("world").cube(world_size).volume(vacuum);

  //Cylinder, acting as the vessel
  auto vessel_steel = n4::tubs("vessel_steel").r(vessel_out_rad_).z(vessel_out_length_).volume(steel);
  n4::place(vessel_steel).in(world).check_overlaps().now();

  //Build inside detector
  //BuildTPC(gas_, mesh_mat, steel, peek, vacuum, quartz, tpb, vessel_steel);

  auto vessel = n4::tubs("GasVessel").r(vessel_rad_).z(vessel_length_).volume(gas_);
  n4::place(vessel).in(vessel_steel).at({0, 0, 0}).check_overlaps().now();

  G4double cathode_z;
  G4double cathBracket_z;
  G4LogicalVolume * gas_el;
  G4double drift_length_  ;
  G4double el_length_    ;
  G4double drift_z;
  G4double drift_r;

  if (model_new_ == 0) {
    cathode_z = 90.1125*mm - 15.745*mm;  //cathode center from vessel center
    //auto cathode_z = 4.505*mm + mesh_thickn_/2 + 2*D + 5*d + 6*ring_thickn_;  //cathode center from vessel center

    // Cathode Bracket
    cathBracket_z = cathode_z;
    //auto cathBracket_z = 8.005*mm - meshBracket_thickn_/2 + 2*D + 5*d + 6*ring_thickn_;

    //Gas
    drift_length_ = 96*mm - meshBracket_thickn_ ;
    el_length_    = 15*mm - anodeBracket_thickn_/2;

    // Drift
    drift_z = cathode_z - mesh_thickn_/2 - drift_length_/2;
    drift_r = meshBracket_rad_;


    G4double d             =   3     * mm;
    G4double D             =   5     * mm;
    G4double ring_rad_int_ = 130 /2. * mm;
    G4double ring_rad_out_ = 140 /2. * mm;
    G4double ring_thickn_  =  10     * mm;

    // ------------------------------------------------------------------------------------------------------------------------
    //Cu rings
    auto ring = n4::tubs("ring").r_inner(ring_rad_int_).r(ring_rad_out_).z(ring_thickn_).place(Cu).in(vessel).check_overlaps();
    auto first_ring_z = cathode_z - meshBracket_thickn_/2 - D - ring_thickn_/2;
    for (auto n : {0,1,2,3,4,5}) {
      ring.clone().at(0, 0, first_ring_z - n * (ring_thickn_ + d)).copy_no(n).now();
    }

    // Source box
    auto source_box_width  = 100*mm;
    auto source_box_length =  50*mm;
    auto source_box_z      = cathode_z + 81.64*mm + source_box_length/2;
    //auto source_box_z = cathode_z ;
    n4::box("source_box").xy(source_box_width).z(source_box_length).place(steel).in(vessel).at(0,0,source_box_z).check_overlaps().now();

  }  else {
    cathode_z = 4.505*mm + mesh_thickn_/2;  //cathode center from vessel center

    //Cathode Bracket
    cathBracket_z = 8.005*mm - meshBracket_thickn_/2;

    //Gas
    drift_length_  = 19.825*mm - mesh_thickn_;
    el_length_     = 10.775*mm + mesh_thickn_;

    // Drift
    drift_z = cathode_z - mesh_thickn_/2 - drift_length_/2;
    drift_r = anodeBracket_rad_;
  }

  // Common volumes
  n4::tubs("cathode"       ).                   r(mesh_rad_       ).z(mesh_thickn_       ).place(mesh_mat).in(vessel).at(0,0,    cathode_z).check_overlaps().now();
  n4::tubs("CathodeBracket").r_inner(mesh_rad_).r(meshBracket_rad_).z(meshBracket_thickn_).place(steel   ).in(vessel).at(0,0,cathBracket_z).check_overlaps().now();
  n4::tubs("gas_drift"     ).                   r(drift_r         ).z(drift_length_      ).place(gas_    ).in(vessel).at(0,0,      drift_z).check_overlaps().now();

  // EL gap
  G4double el_z = drift_z - drift_length_/2 - el_length_/2;
  gas_el = n4::tubs("gas_el").r(anodeBracket_rad_).z(el_length_).volume(gas_);
  n4::place(gas_el).in(vessel).at({0, 0, el_z}).check_overlaps().now();

  // Gate
  G4double gate_z = el_length_/2 - mesh_thickn_/2;
  n4::tubs("gate").r(mesh_rad_).z(mesh_thickn_).place(mesh_mat).in(gas_el).at(0, 0, gate_z).check_overlaps().now();

  // Gate Bracket
  G4double gateBracket_z = 12.745*mm + meshBracket_thickn_/2;
  n4::tubs("gateBracket").r_inner(mesh_rad_).r(meshBracket_rad_).z(meshBracket_thickn_).place(steel).in(vessel).at(0, 0, -gateBracket_z).check_overlaps().now();

  //Anode
  G4double anode_z = - el_length_/2 + mesh_thickn_/2;
  n4::tubs("Anode").r(mesh_rad_).z(mesh_thickn_).place(mesh_mat).in(gas_el).at(0., 0., anode_z).check_overlaps().now();

  //Anode Bracket
  G4double anodeBracket_z = gateBracket_z + meshBracket_thickn_/2 + 3.775*mm + anodeBracket_thickn_/2;
  n4::tubs("AnodeBracket").r_inner(mesh_rad_).r(anodeBracket_rad_).z(anodeBracket_thickn_).place(steel).in(vessel).at(0, 0, -anodeBracket_z).check_overlaps().now();

  // Peek Anode Holder
  G4double anodeHolder_length_  = 30        *mm;
  G4double anodeHolder_rad_     = 145.001/2 *mm;
  G4double anodeHolder_thickn_  = 17.       *mm;
  G4double anodeHolder_angle_   = 7.964     *deg;

  G4double anodeHolder_hole_length_ = 20    *mm;
  G4double anodeHolder_hole_rad_    = anodeHolder_rad_ + 8*mm;
  G4double anodeHolder_hole_thickn_ = 10    *mm;
  G4double anodeHolder_hole_angle_  = 7.965 *deg;


  auto anode_holder = n4::tubs("AnodeHolder"    ).r_inner(anodeHolder_rad_     ).r_delta(anodeHolder_thickn_     ).z(anodeHolder_length_     ).phi_start(-anodeHolder_angle_     /2).phi_delta(anodeHolder_angle_     )
    .subtract(        n4::tubs("AnodeHolderHole").r_inner(anodeHolder_hole_rad_).r_delta(anodeHolder_hole_thickn_).z(anodeHolder_hole_length_).phi_start(-anodeHolder_hole_angle_/2).phi_delta(anodeHolder_hole_angle_))
    .at(0, 0, anodeHolder_length_/2 - anodeHolder_hole_length_/2)
    .volume(peek);


  G4RotationMatrix* Rot45   = new G4RotationMatrix(); Rot45   -> rotateZ(  45*deg);
  G4RotationMatrix* Rot_45  = new G4RotationMatrix(); Rot_45  -> rotateZ( -45*deg);
  G4RotationMatrix* Rot135  = new G4RotationMatrix(); Rot135  -> rotateZ( 135*deg);
  G4RotationMatrix* Rot_135 = new G4RotationMatrix(); Rot_135 -> rotateZ(-135*deg);

  G4double anodeHolder_z = 29.495*mm + anodeHolder_length_/2;
  for (auto [i, angle] : enumerate({45, 135, -45, -135}) ) {
    n4::place(anode_holder).in(vessel).rotate_z(angle*deg).at(0, 0, -anodeHolder_z).copy_no(i).check_overlaps().now();
  }

  /////////////////////////
  //Quartz Window
  G4double quartz_window_z = 35.495*mm + quartz_window_thickn_/2 ;
  n4::tubs("QuartzWindow").r(quartz_window_rad_).z(quartz_window_thickn_).place(quartz).in(vessel).at(0., 0., -quartz_window_z).check_overlaps().now();

  //Evaporated TPB
  G4double tpb_coating_z = quartz_window_z - quartz_window_thickn_/2 - tpb_coating_thickn_/2 ;
  n4::tubs("CoatingTPB").r(quartz_window_rad_).z(tpb_coating_thickn_).place(tpb).in(vessel).at(0., 0., -tpb_coating_z).check_overlaps().now();

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

  auto window_phi = n4::tubs("").phi_start(-windowHolder_angle_/2).phi_delta(windowHolder_angle_);


  auto logic_windowHolder =
   /*   */window_phi        .name("1").r_inner(windowHolder_rad0_).r(windowHolder_rad2_).z(windowHolder_length1_)
    .add( window_phi.clone().name("2").r_inner(windowHolder_rad1_).r(windowHolder_rad2_).z(windowHolder_length2_) ).at(0, 0, -windowHolder_block2_z)
    .add( window_phi.clone().name("3").r_inner(windowHolder_rad1_).r(windowHolder_rad3_).z(windowHolder_length3_) ).at(0, 0, -windowHolder_block3_z)
    .name("QuartzWindowHolder")
    .volume(peek);


  G4double quartz_windowHolder_z = quartz_window_z - quartz_window_thickn_/2 -  windowHolder_length1_/2;
  for (auto [i, angle] : enumerate({45, 135, -45, -135})) {
      n4::place(logic_windowHolder).in(vessel).rotate_z(angle*deg).at(0,0, -quartz_windowHolder_z).copy_no(i).check_overlaps().now();
  }


  G4double windowHolderTop_length1_ = 3.  *mm;
  G4double windowHolderTop_length2_ = 4.  *mm;
  G4double windowHolderTop_rad0_    = windowHolder_rad0_ ;
  G4double windowHolderTop_rad1_    = windowHolderTop_rad0_ + 4.*mm ;
  G4double windowHolderTop_rad2_    = windowHolder_rad1_ ;

  G4double windowHolderTop_block2_z = windowHolderTop_length1_/2 + windowHolderTop_length2_/2 ;

  auto holder_phi = n4::tubs("").phi_start(-windowHolder_angle_/2).phi_delta(windowHolder_angle_);
  auto logic_windowHolderTop =
  /*   */holder_phi.        name("QuartzWindowHolderTop1").r_inner(windowHolderTop_rad1_).r(windowHolderTop_rad2_).z(windowHolderTop_length1_)
    .add(holder_phi.clone().name("QuartzWindowHolderTop2").r_inner(windowHolderTop_rad0_).r(windowHolderTop_rad2_).z(windowHolderTop_length2_)).at(0,0, -windowHolderTop_block2_z)
    .name("WindowHolderTop_Union")
    .volume(peek);


  G4double quartz_windowHolderTop_z = quartz_window_z - quartz_window_thickn_/2 +  windowHolderTop_length1_/2;
  for (auto [i, angle] : enumerate({45, 135, -45, -135})) {
    n4::place(logic_windowHolderTop).in(vessel).rotate_z(angle*deg).at(0, 0, -quartz_windowHolderTop_z).copy_no(i).check_overlaps().now();
  }

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
    pos_pmt              = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm,            -pmt_z);
    pos_enclosure_pmt    = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm,    relative_pmt_z);
    pos_enclosurevac_pmt = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, relativevac_pmt_z);
    pos                  = G4ThreeVector(pmt_PsX[i]*mm, pmt_PsY[i]*mm, 0);

    solid_enclosure_pmt    = new G4SubtractionSolid("EnclosurePMT_Sub"    , solid_enclosure_pmt    , solid_pmt, 0, pos_enclosure_pmt);
    solid_enclosurevac_pmt = new G4SubtractionSolid("EnclosureVacPMT_Sub" , solid_enclosurevac_pmt , solid_pmt, 0, pos_enclosurevac_pmt);
    solid_pmtHolder        = new G4SubtractionSolid("PMTHolder_Sub"       , solid_pmtHolder        , solid_pmt, 0, pos);
    solid_plate1_pmt       = new G4SubtractionSolid("PMTplateBottom1_Sub" , solid_plate1_pmt       , solid_pmt, 0, pos);
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
