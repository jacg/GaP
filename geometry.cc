#include "materials.hh"

#include "nain4.hh"
#include "g4-mandatory.hh"
#include "geometry.hh"
#include "n4-volumes.hh"
#include "n4-utils.hh"

#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4OpticalPhysics.hh>
#include <G4PVPlacement.hh>
#include <G4RandomDirection.hh>
#include <G4RunManagerFactory.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>

#include <G4VSolid.hh>
#include <iostream>


const auto world_size = 0.5 * m;

const auto vessel_out_rad    = 288./2  *mm;
const auto vessel_out_length = 46.679  *cm;
const auto vessel_rad        = 276./2  *mm;
const auto vessel_length     = 38.599  *cm; // Adjusted length so that the gas volume is centered. Original length (38.639  *cm),

const auto mesh_rad          = 104./2  *mm;
const auto mesh_thickn       = 0.075   *mm;
const auto mesh_transparency = 0.95;

const auto meshBracket_rad      = 180./2  *mm;
const auto meshBracket_thickn   = 6.      *mm;
const auto anodeBracket_rad     = 160./2  *mm;
const auto anodeBracket_thickn  = 6.975   *mm;

const auto pmt_rad = 25.4/2  *mm;

const auto enclosure_pmt_rad        = 120./2  *mm;
const auto enclosure_pmt_thickn     = 8.5     *mm;
const auto enclosure_pmt_length     = 113.5   *mm;
const auto enclosurevac_pmt_length  = 110.5   *mm;

const auto plate_pmt_rad        = 105./2  *mm;
const auto plate_pmt_thickn     = 105./2  *mm;
const auto plate_pmt_length     = 10      *mm;
const auto plateUp_pmt_length   = 15      *mm;
const auto plateUp_pmt_thickn   = 21.5    *mm;

const auto pmtHolder_rad        = 115./2  *mm;
const auto pmtHolder_length     = 9       *mm;

const auto quartz_window_rad    = 108./2  *mm;
const auto quartz_window_thickn = 3       *mm;
const auto tpb_coating_thickn   = 3       *micrometer;

const auto photoe_prob       =   0;
const auto pressure          =  10 * bar;
const auto temperature       = 293 * kelvin;
//const auto sc_yield        =  22222./MeV; // Wsc = 45 eV, fr
const auto sc_yield          =   1/GeV;
//const auto sc_yield        =  1000./MeV;
const auto elifetime         = 1e6 * ms;
// const auto drift_vel         = 1. * mm/microsecond;
// const auto drift_transv_diff = 1. * mm/sqrt(cm);
// const auto drift_long_diff   = .3 * mm/sqrt(cm);
// const auto el_field          = 16.0 * kilovolt/cm;
// const auto el_vel            = 3. * mm/microsecond;
// const auto el_transv_diff    = 1. * mm/sqrt(cm);
// const auto el_long_diff      = .3 * mm/sqrt(cm);

G4Material* peek;
G4Material* steel;
G4Material* Cu;
G4Material* vacuum;
G4Material* mesh_mat;
G4Material* quartz;
G4Material* tpb;
G4Material* gas;
G4LogicalVolume* world;

// Must call this at the start of any sub-geometry that you want to visualise
// without having to construct the rest of the geometry.
void ensure_initialized() {
  static bool initialized = false;
  if (initialized) { return; }
  initialized = true;

  Cu     = n4::material("G4_Cu");
  vacuum = n4::material("G4_Galactic");
  steel  = n4::material("G4_STAINLESS-STEEL");

  gas      = GAr_with_properties( pressure, temperature, sc_yield, elifetime);
  mesh_mat = FakeDielectric_with_properties(gas, "mesh_mat",
                                            pressure, temperature, mesh_transparency, mesh_thickn,
                                            sc_yield, elifetime, photoe_prob);
  peek   = peek_with_properties();
  quartz = quartz_with_properties();
  tpb    = TPB_with_properties();
  world = n4::box("world").cube(world_size).volume(vacuum);
}

G4LogicalVolume* get_world() {
  ensure_initialized();
  return world;
}


void place_mesh_holder_in(G4LogicalVolume* vessel) {
  // Ensure this sub-geometry can be built (probably for visualization) without
  // needing to construct the rest of the geometry
  ensure_initialized();

  // Peek Mesh Holder (holds cathode and gate)
  auto meshHolder_length    = 36.75     *mm;
  auto meshHolder_width     = 21.035    *mm;
  auto meshHolder_thickn    = 24        *mm;

  auto meshHolder_hole_length   = 15.75  *mm;
  auto meshHolder_holeUp_length = 6      *mm;
  auto meshHolder_hole_thickn   = 11.667 *mm;

  auto dz_down =   meshHolder_length/2 - meshHolder_hole_length  /2 - 5*mm;
  auto dz_up   = -(meshHolder_length/2 - meshHolder_holeUp_length/2 - 5*mm);

  auto logic_meshHolder =
    /*      */n4::box("MeshHolder"        ).xyz(meshHolder_width  , meshHolder_thickn       , meshHolder_length       )
    .subtract(n4::box("MEshHolderHoleDown").xyz(meshHolder_width*2, meshHolder_hole_thickn*2, meshHolder_hole_length  )).at(0, meshHolder_thickn/2, dz_down)
    .subtract(n4::box("MEshHolderHoleUp"  ).xyz(meshHolder_width*2, meshHolder_hole_thickn*2, meshHolder_holeUp_length)).at(0, meshHolder_thickn/2, dz_up  )
    .name("MeshHolder")
    .volume(peek);

  auto meshHolder_x = 5*mm + meshBracket_rad * cos(45*deg);
  auto meshHolder_y = 5*mm + meshBracket_rad * sin(45*deg);
  auto meshHolder_z = -13.005*mm + meshHolder_length/2 ;

  auto mesh_holder = n4::place(logic_meshHolder).in(vessel)
    .rotate_z(135*deg)                             // Orient the mesh holder
    .at(meshHolder_x, meshHolder_y, -meshHolder_z) // Displace from the center
    .check_overlaps();

  for (auto i : {0,1,2,3}) {
    mesh_holder                                    // Without .clone() rotations are cumulative
      .rotate_z(90 * deg)                          // Add 90 deg
      .copy_no(i)
      .now();
  }

  //Steel Bar joining Mesh holder and PMT clad
  auto meshHolderBar_rad     =  9./2  *mm;
  auto meshHolderBar_length  = 35.75  *mm;
  auto meshHolderBar_xy      = 73.769 *mm;
  auto meshHolderBar_z  = meshHolder_z + meshHolder_length/2 + meshHolderBar_length/2 ;


  auto meshHolderBar = n4::tubs("MeshHolderBar").r(meshHolderBar_rad).z(meshHolderBar_length).volume(steel);
  auto i = 1;
  for   (auto y : {meshHolderBar_xy, -meshHolderBar_xy}) {
    for (auto x : {meshHolderBar_xy, -meshHolderBar_xy}) {
      n4::place(meshHolderBar).in(vessel).at(x, y, -meshHolderBar_z).copy_no(i++).check_overlaps().now();
    }
  }
}

void place_quartz_window_holder_in(G4LogicalVolume* vessel) {
  /////////////////////////
  //Quartz Window and Evaporated TPB
  auto quartz_window_z = 35.495*mm + quartz_window_thickn/2;
  auto tpb_coating_z = quartz_window_z - quartz_window_thickn/2 - tpb_coating_thickn/2;
  n4::tubs("QuartzWindow").r(quartz_window_rad).z(quartz_window_thickn).place(quartz).in(vessel).at(0, 0, -quartz_window_z).check_overlaps().now();
  n4::tubs("CoatingTPB"  ).r(quartz_window_rad).z(  tpb_coating_thickn).place(tpb   ).in(vessel).at(0, 0,   -tpb_coating_z).check_overlaps().now();

  // Peek Quartz Window Holder

  auto windowHolder_angle   =  10      *deg;
  auto windowHolder_length1 =   4      *mm;
  auto windowHolder_length2 =  19      *mm;
  auto windowHolder_length3 =   5      *mm;
  auto windowHolder_rad0    = 100.4/2  *mm;
  auto windowHolder_rad1    = 118.4/2  *mm;
  auto windowHolder_rad2    = 122.4/2  *mm;
  auto windowHolder_rad3    = 140.4/2  *mm;

  auto windowHolder_block2_z = windowHolder_length1/2 + windowHolder_length2/2 ;
  auto windowHolder_block3_z = windowHolder_length3/2 + windowHolder_length2/2 + windowHolder_block2_z ;

  auto window_phi = n4::tubs("").phi_start(-windowHolder_angle/2).phi_delta(windowHolder_angle);

  auto logic_windowHolder =
    /*  */window_phi        .name("1").r_inner(windowHolder_rad0).r(windowHolder_rad2).z(windowHolder_length1)
    .add( window_phi.clone().name("2").r_inner(windowHolder_rad1).r(windowHolder_rad2).z(windowHolder_length2) ).at(0, 0, -windowHolder_block2_z)
    .add( window_phi.clone().name("3").r_inner(windowHolder_rad1).r(windowHolder_rad3).z(windowHolder_length3) ).at(0, 0, -windowHolder_block3_z)
    .name("QuartzWindowHolder")
    .volume(peek);


  auto quartz_windowHolder_z = quartz_window_z - quartz_window_thickn/2 -  windowHolder_length1/2;
  for (auto [i, angle] : enumerate({45, 135, -45, -135})) {
    n4::place(logic_windowHolder).in(vessel).rotate_z(angle*deg).at(0,0, -quartz_windowHolder_z).copy_no(i).check_overlaps().now();
  }

  auto windowHolderTop_length1 = 3  *mm;
  auto windowHolderTop_length2 = 4  *mm;
  auto windowHolderTop_rad0    = windowHolder_rad0;
  auto windowHolderTop_rad1    = windowHolderTop_rad0 + 4*mm;
  auto windowHolderTop_rad2    = windowHolder_rad1;

  auto windowHolderTop_block2_z = windowHolderTop_length1/2 + windowHolderTop_length2/2;

  auto holder_phi = n4::tubs("").phi_start(-windowHolder_angle/2).phi_delta(windowHolder_angle);
  auto logic_windowHolderTop =
    /* */holder_phi.        name("QuartzWindowHolderTop1").r_inner(windowHolderTop_rad1).r(windowHolderTop_rad2).z(windowHolderTop_length1)
    .add(holder_phi.clone().name("QuartzWindowHolderTop2").r_inner(windowHolderTop_rad0).r(windowHolderTop_rad2).z(windowHolderTop_length2)).at(0,0, -windowHolderTop_block2_z)
    .name("WindowHolderTop_Union")
    .volume(peek);


  auto quartz_windowHolderTop_z = quartz_window_z - quartz_window_thickn/2 +  windowHolderTop_length1/2;
  for (auto [i, angle] : enumerate({45, 135, -45, -135})) {
    n4::place(logic_windowHolderTop).in(vessel).rotate_z(angle*deg).at(0, 0, -quartz_windowHolderTop_z).copy_no(i).check_overlaps().now();
  }
};

void place_pmt_holder_in(G4LogicalVolume* vessel) {
  //Build PMT
  //auto pmt_length = pmt_.Length();
  auto pmt_length = 43.0  *mm;
  auto pmt_z      = 42.495*mm + pmt_length/2;

  G4Tubs* solid_pmt = n4::tubs("SolidPMT").r(pmt_rad).z(pmt_length).solid(); // Hamamatsu pmt length: 43*mm | STEP pmt gap length: 57.5*mm

  // Position pairs (x,Y) for PMTs
  std::vector <float> pmt_PsX={-15.573,  20.68 , -36.253, 0, 36.253, -20.68 , 15.573};
  std::vector <float> pmt_PsY={-32.871, -29.922,  -2.949, 0,  2.949,  29.922, 32.871};

  // PMT clad
  G4VSolid *solid_enclosure_pmt = n4::tubs("EnclosurePMT").r(enclosure_pmt_rad + enclosure_pmt_thickn).z(enclosure_pmt_length).solid();
  auto enclosure_pmt_z = vessel_length/2 - enclosure_pmt_length/2;
  auto relative_pmt_z  = enclosure_pmt_z - pmt_z;
  // Vacuum inside the pmt enclosure
  auto enclosurevac_pmt_z = vessel_length/2 - enclosurevac_pmt_length/2;
  auto relativevac_pmt_z  = enclosurevac_pmt_z - pmt_z;
  G4VSolid* solid_enclosurevac_pmt = n4::tubs("EnclosureVacPMT").r(enclosure_pmt_rad).z(enclosurevac_pmt_length).solid();

  // PMT Holder
  auto pmtHolder_z = enclosurevac_pmt_length/2 - pmtHolder_length/2;
  G4VSolid* solid_pmtHolder = n4::tubs("PMTHolder").r(pmtHolder_rad).z(pmtHolder_length).solid();
  // Steel plate enclosing the pmt tube
  auto plate1_pmt_z = enclosure_pmt_z - enclosure_pmt_length/2 - plate_pmt_length/2;
  G4VSolid* solid_plate1_pmt = n4::tubs("PMTplateBottom1").r(plate_pmt_rad+plate_pmt_thickn).z(plate_pmt_thickn).solid();

  G4ThreeVector pos_pmt              = {0, 0, 0};
  G4ThreeVector pos_enclosure_pmt    = {0, 0, 0};
  G4ThreeVector pos_enclosurevac_pmt = {0, 0, 0};
  G4ThreeVector pos                  = {0, 0, 0};

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
  n4::place(logic_enclosure_pmt).in(vessel).at(0, 0, -enclosure_pmt_z).check_overlaps().now();

  G4LogicalVolume *logic_enclosurevac_pmt = new G4LogicalVolume(solid_enclosurevac_pmt, vacuum, "EnclosureVacPMT");
  n4::place(logic_enclosurevac_pmt).in(logic_enclosure_pmt).at(0, 0, enclosure_pmt_z - enclosurevac_pmt_z).check_overlaps().now();

  G4LogicalVolume *logic_pmtHolder = new G4LogicalVolume(solid_pmtHolder, steel, "PMTHolder");
  n4::place(logic_pmtHolder).in(logic_enclosurevac_pmt).at(0, 0, pmtHolder_z).check_overlaps().now();

  G4LogicalVolume *logic_plate1_pmt = new G4LogicalVolume(solid_plate1_pmt, steel, "PMTplateBottom1");
  n4::place(logic_plate1_pmt).in(vessel).at(0, 0, -plate1_pmt_z).check_overlaps().now();

  // Steel plate attached where the peek holders are attached
  auto plate0_pmt = n4::tubs("PMTplateBottom0").r_inner(plate_pmt_rad).r_delta(plateUp_pmt_thickn).z(plateUp_pmt_length).volume(steel);
  auto plate0_pmt_z = plate1_pmt_z - plate_pmt_length;
  n4::place(plate0_pmt).in(vessel).at(0, 0, -plate0_pmt_z).check_overlaps().now();

  // Upper steel plate at the pmt clad
  auto plateUp_pmt_rad = enclosure_pmt_rad + enclosure_pmt_thickn;
  auto plateUp_pmt = n4::tubs("PMTplateUp").r_inner(plateUp_pmt_rad).r_delta(plateUp_pmt_thickn).z(plateUp_pmt_length).volume(steel);
  auto plateUp_pmt_z = vessel_length/2 - plateUp_pmt_length/2 ;
  n4::place(plateUp_pmt).in(vessel).at(0, 0, -plateUp_pmt_z).check_overlaps().now();
};

G4PVPlacement* geometry() {
  ensure_initialized();

  auto model_new= 1;

  //Cylinder, acting as the vessel
  auto vessel_steel = n4::tubs("vessel_steel").r(vessel_out_rad).z(vessel_out_length).volume(steel);
  n4::place(vessel_steel).in(world).check_overlaps().now();

  //Build inside detector
  //BuildTPC(gas, mesh_mat, steel, peek, vacuum, quartz, tpb, vessel_steel);

  auto vessel = n4::tubs("GasVessel").r(vessel_rad).z(vessel_length).volume(gas);
  n4::place(vessel).in(vessel_steel).at(0, 0, 0).check_overlaps().now();

  G4double cathode_z;
  G4double cathBracket_z;
  G4LogicalVolume * gas_el;
  G4double drift_length;
  G4double el_length;
  G4double drift_z;
  G4double drift_r;

  if (model_new == 0) {
    cathode_z = 90.1125*mm - 15.745*mm;  //cathode center from vessel center
    //auto cathode_z = 4.505*mm + mesh_thickn/2 + 2*D + 5*d + 6*ring_thickn;  //cathode center from vessel center

    // Cathode Bracket
    cathBracket_z = cathode_z;
    //auto cathBracket_z = 8.005*mm - meshBracket_thickn/2 + 2*D + 5*d + 6*ring_thickn;

    //Gas
    drift_length = 96*mm - meshBracket_thickn ;
    el_length    = 15*mm - anodeBracket_thickn/2;

    // Drift
    drift_z = cathode_z - mesh_thickn/2 - drift_length/2;
    drift_r = meshBracket_rad;


    auto d             =  3     * mm;
    auto D             =  5     * mm;
    auto ring_rad_int = 130 /2. * mm;
    auto ring_rad_out = 140 /2. * mm;
    auto ring_thickn  =  10     * mm;

    // ------------------------------------------------------------------------------------------------------------------------
    //Cu rings
    auto ring = n4::tubs("ring").r_inner(ring_rad_int).r(ring_rad_out).z(ring_thickn).place(Cu).in(vessel).check_overlaps();
    auto first_ring_z = cathode_z - meshBracket_thickn/2 - D - ring_thickn/2;
    for (auto n : {0,1,2,3,4,5}) {
      ring.clone().at(0, 0, first_ring_z - n * (ring_thickn + d)).copy_no(n).now();
    }

    // Source box
    auto source_box_width  = 100*mm;
    auto source_box_length =  50*mm;
    auto source_box_z      = cathode_z + 81.64*mm + source_box_length/2;
    //auto source_box_z = cathode_z ;
    n4::box("source_box").xy(source_box_width).z(source_box_length).place(steel).in(vessel).at(0,0,source_box_z).check_overlaps().now();

  } else {
    cathode_z = 4.505*mm + mesh_thickn/2;  //cathode center from vessel center

    //Cathode Bracket
    cathBracket_z = 8.005*mm - meshBracket_thickn/2;

    //Gas
    drift_length  = 19.825*mm - mesh_thickn;
    el_length     = 10.775*mm + mesh_thickn;

    // Drift
    drift_z = cathode_z - mesh_thickn/2 - drift_length/2;
    drift_r = anodeBracket_rad;
  }

  // Common volumes
  n4::tubs("cathode"       ).                  r(mesh_rad       ).z(mesh_thickn       ).place(mesh_mat).in(vessel).at(0,0,    cathode_z).check_overlaps().now();
  n4::tubs("CathodeBracket").r_inner(mesh_rad).r(meshBracket_rad).z(meshBracket_thickn).place(steel   ).in(vessel).at(0,0,cathBracket_z).check_overlaps().now();
  n4::tubs("gas_drift"     ).                  r(drift_r        ).z(drift_length      ).place(gas     ).in(vessel).at(0,0,      drift_z).check_overlaps().now();

  // EL gap
  auto el_z = drift_z - drift_length/2 - el_length/2;
  gas_el = n4::tubs("gas_el").r(anodeBracket_rad).z(el_length).volume(gas);
  n4::place(gas_el).in(vessel).at(0, 0, el_z).check_overlaps().now();

  // Gate
  auto gate_z = el_length/2 - mesh_thickn/2;
  n4::tubs("gate").r(mesh_rad).z(mesh_thickn).place(mesh_mat).in(gas_el).at(0, 0, gate_z).check_overlaps().now();

  // Gate Bracket
  auto gateBracket_z = 12.745*mm + meshBracket_thickn/2;
  n4::tubs("gateBracket").r_inner(mesh_rad).r(meshBracket_rad).z(meshBracket_thickn).place(steel).in(vessel).at(0, 0, -gateBracket_z).check_overlaps().now();

  //Anode
  auto anode_z = - el_length/2 + mesh_thickn/2;
  n4::tubs("Anode").r(mesh_rad).z(mesh_thickn).place(mesh_mat).in(gas_el).at(0, 0, anode_z).check_overlaps().now();

  //Anode Bracket
  auto anodeBracket_z = gateBracket_z + meshBracket_thickn/2 + 3.775*mm + anodeBracket_thickn/2;
  n4::tubs("AnodeBracket").r_inner(mesh_rad).r(anodeBracket_rad).z(anodeBracket_thickn).place(steel).in(vessel).at(0, 0, -anodeBracket_z).check_overlaps().now();

  // Peek Anode Holder
  auto anodeHolder_length  =  30       *mm;
  auto anodeHolder_rad     = 145.001/2 *mm;
  auto anodeHolder_thickn  =  17       *mm;
  auto anodeHolder_angle   =   7.964   *deg;

  auto anodeHolder_hole_rad    = anodeHolder_rad + 8*mm;
  auto anodeHolder_hole_length = 20      *mm;
  auto anodeHolder_hole_thickn = 10      *mm;
  auto anodeHolder_hole_angle  =  7.965 *deg;

  auto anode_holder = n4::tubs("AnodeHolder"    ).r_inner(anodeHolder_rad     ).r_delta(anodeHolder_thickn     ).z(anodeHolder_length     ).phi_start(-anodeHolder_angle     /2).phi_delta(anodeHolder_angle     )
    .subtract(        n4::tubs("AnodeHolderHole").r_inner(anodeHolder_hole_rad).r_delta(anodeHolder_hole_thickn).z(anodeHolder_hole_length).phi_start(-anodeHolder_hole_angle/2).phi_delta(anodeHolder_hole_angle))
    .at(0, 0, anodeHolder_length/2 - anodeHolder_hole_length/2)
    .volume(peek);


  G4RotationMatrix* Rot45   = new G4RotationMatrix(); Rot45   -> rotateZ(  45*deg);
  G4RotationMatrix* Rot_45  = new G4RotationMatrix(); Rot_45  -> rotateZ( -45*deg);
  G4RotationMatrix* Rot135  = new G4RotationMatrix(); Rot135  -> rotateZ( 135*deg);
  G4RotationMatrix* Rot_135 = new G4RotationMatrix(); Rot_135 -> rotateZ(-135*deg);

  auto anodeHolder_z = 29.495*mm + anodeHolder_length/2;
  for (auto [i, angle] : enumerate({45, 135, -45, -135}) ) {
    n4::place(anode_holder).in(vessel).rotate_z(angle*deg).at(0, 0, -anodeHolder_z).copy_no(i).check_overlaps().now();
  }

  place_quartz_window_holder_in(vessel);

  place_pmt_holder_in(vessel);
  if (model_new == 1) {

    // A
    //auto A_length = 56 *mm;
    //auto A_width  = 27 *mm;
    //auto A_thickn = 21 *mm;

    //auto A_hole_length = 50 *mm;
    //auto A_hole_width  = 16 *mm;

    //G4Box *solid_A_block = new G4Box("ABlock", A_width/2, A_thickn/2 , A_length/2);
    //G4Box *solid_A_hole = new G4Box("AHole", A_hole_width/2, A_thickn/2 , A_hole_length/2);
    //G4VSolid        *solid_A= new G4SubtractionSolid("SolidA", solid_A_block, solid_A_hole, 0, G4ThreeVector(-(A_width - A_hole_width)/2,0.,-(A_length - A_hole_length)/2) );
    //G4LogicalVolume *logic_A = new G4LogicalVolume(solid_A, peek, "A");

    //auto A_xy = plate_pmt_rad+plate_pmt_thickn - A_length/2;
    //auto A_z = plate0_pmt_z + plate_pmt_length/2 - A_length/2;

    //n4::place(logic_A).in(vessel).rotate(*Rot_45).at(-A_xy, -A_xy, -A_z).copy_no(0).check_overlaps().now();
    //n4::place(logic_A).in(vessel).rotate(*Rot_135).at(-A_xy, A_xy, -A_z).copy_no(1).check_overlaps().now();
    //n4::place(logic_A).in(vessel).rotate(*Rot45).at(A_xy, -A_xy, -A_z).copy_no(2).check_overlaps().now();
    //n4::place(logic_A).in(vessel).rotate(*Rot135).at(A_xy, A_xy, -A_z).copy_no(3).check_overlaps().now();
  } else {
    place_mesh_holder_in(vessel);
  }
  return n4::place(world).now();
}
