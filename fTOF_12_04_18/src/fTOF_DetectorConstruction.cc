//My
#include "fTOF_DetectorConstruction.hh"
#include "fTOF_SensitiveDetector.hh"
#include "MagneticField.hh"

//G4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Color.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "globals.hh"
//magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4UserLimits.hh"
//GDML
//#include <G4GDMLParser.hh>

//root 
#include "TMath.h"

#include "crtConst.hh"


fTOF_DetectorConstruction::fTOF_DetectorConstruction()
{
  //magField = new MagneticField();
  worldVisAtt = new G4VisAttributes();
  quartzVisAtt = new G4VisAttributes();
  sensitiveVisAtt = new G4VisAttributes();
  pmtboxVisAtt = new G4VisAttributes();
  absVisAtt = new G4VisAttributes();
  // Define Materials to be used
  DefineMaterials();
}

fTOF_DetectorConstruction::~fTOF_DetectorConstruction()
{
  //delete magField;
  delete worldVisAtt;
  delete quartzVisAtt;
  delete sensitiveVisAtt;
  delete pmtboxVisAtt;
  delete absVisAtt;
}

void fTOF_DetectorConstruction::DefineMaterials()
{
  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;

  // Define elements
  //G4Element* H = 
  //new G4Element("Hydrogen", symbol = "H", z = 1., a = 1.01*g/mole);
  G4Element* C = 
    new G4Element("Carbon", symbol = "C", z = 6., a = 12.01*g/mole);
  G4Element* N = 
    new G4Element("Nitrogen", symbol = "N", z = 7., a = 14.01*g/mole);
  G4Element* O =
    new G4Element("Oxygen", symbol = "O", z = 8., a = 16.00*g/mole);
  G4Element* Si = 
    new G4Element("Silicon", symbol = "Si", z = 14., a = 28.09*g/mole);
  G4Element* Al = 
    new G4Element("Aluminum", symbol = "Al", z = 13., a = 26.98*g/mole);

  // Quartz Material (SiO2)
  G4Material* SiO2 = 
    new G4Material("quartz", density = 2.200*g/cm3, ncomponents = 2);
  SiO2->AddElement(Si, natoms = 1);
  SiO2->AddElement(O , natoms = 2);
  
  // Air
  //G4Material* Air = 
  //new G4Material("Air", density = 1.290*mg/cm3, ncomponents = 2);
  G4Material* Air = 
    new G4Material("Air", density = 0.000290*mg/cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  // Aluminum
  G4Material* Aluminum =
    new G4Material("Aluminum", density = 2.7*g/cm3, ncomponents = 1);
  Aluminum->AddElement(Al, fractionmass = 1.0);


  // Borosilicate
  // G4Material* Borosilicate = 
  //   new G4Material("Borosilicate glass", density= 2.23*g/cm3, ncomponents=6);
  // Borosilicate->AddElement(B, fractionmass=0.040064);
  // Borosilicate->AddElement(O, fractionmass=0.539562); 
  // Borosilicate->AddElement(Na, fractionmass=0.028191);
  // Borosilicate->AddElement(Al, fractionmass=0.011644);
  // Borosilicate->AddElement(Si, fractionmass=0.377220);
  // Borosilicate->AddElement(K, fractionmass=0.003321);




  // Assign Materials
  world.material = Air;
  sec.material = SiO2;
  sensitive.material = sec.material;
  pmtWin1.material = sec.material;
  pmtWin2.material = sec.material;  
  pmtAbs1.material = Aluminum;
  abs1.material = Aluminum;
  abs2.material = Aluminum;
  //sensitive.material = Aluminum;
  pmtbox.material = Aluminum;
  //abs1.material = Aluminum;
  //abs2.material = Aluminum;


  // fTOF bar 10.04.18
  barBox.material = SiO2;
  trdXY.material = SiO2;
  trdYZ.material = SiO2;
  trdZX.material = SiO2;

  hamWin.material = SiO2;
  planWin.material = SiO2;

  hamChan.material = SiO2;
  planChan.material = SiO2;

  hamBox.material = Aluminum;
  planBox.material = Aluminum;

  //
  // Generate and Add Material Properties Table
  //						
  const G4int num = 36;
  G4double WaveLength[num];
  G4double Absorption[num];      // Default value for absorption
  G4double AirAbsorption[num];
  G4double AirRefractiveIndex[num];
  G4double PhotonEnergy[num];

  // Absorption of quartz per 1m
  G4double QuartzAbsorption[num] =
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541 ,0.998561611,0.998435332,0.998294892,
     0.998138345,0.997963425,0.997767484,0.997547418,0.99729958 ,
     0.99701966 ,0.99670255 ,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,
     0.990610945};

  for (int i=0; i<num; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    Absorption[i] = 100*m;      // Fake number for no absorption
    AirAbsorption[i] = 4.*cm;   // If photon hits air, kill it
    AirRefractiveIndex[i] = 1.;
    PhotonEnergy[num - (i+1)] = twopi*hbarc/WaveLength[i];
    /* Absorption is given per length and G4 needs mean free path
       length, calculate it here
       mean free path length - taken as probablility equal 1/e
       that the photon will be absorbed */
    QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
    //EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*
    //epotekBarJoint.thickness;
  }

  G4double QuartzRefractiveIndex[num] =
    {1.456535,1.456812,1.4571  ,1.457399,1.457712,1.458038,
     1.458378,1.458735,1.459108,1.4595  ,1.459911,1.460344,
     1.460799,1.46128 ,1.461789,1.462326,1.462897,1.463502,
     1.464146,1.464833,1.465566,1.46635 ,1.46719 ,1.468094,
     1.469066,1.470116,1.471252,1.472485,1.473826,1.475289,
     1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};


  // Assign absorption and refraction to materials

  // Quartz
  G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX", PhotonEnergy, QuartzRefractiveIndex, num);
  QuartzMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);
  
  // Assign this material to the bars
  sec.material->SetMaterialPropertiesTable(QuartzMPT);


  // fTOF MPTs 11.04.18
  barBox.material->SetMaterialPropertiesTable(QuartzMPT);
  trdXY.material->SetMaterialPropertiesTable(QuartzMPT);
  trdZX.material->SetMaterialPropertiesTable(QuartzMPT);
  trdYZ.material->SetMaterialPropertiesTable(QuartzMPT);

  hamWin.material->SetMaterialPropertiesTable(QuartzMPT);
  planWin.material->SetMaterialPropertiesTable(QuartzMPT);

  hamChan.material->SetMaterialPropertiesTable(QuartzMPT);
  planChan.material->SetMaterialPropertiesTable(QuartzMPT);




  // Air
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption, num);
  
  // Assign these properties to the world volume
  world.material->SetMaterialPropertiesTable(AirMPT);

}

G4VPhysicalVolume* fTOF_DetectorConstruction::Construct()
{

  /*
  //magnetic field
  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized){
    G4FieldManager* fieldMgr
      = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
    fieldMgr->GetChordFinder()->SetDeltaChord(1.0*mm);
    fieldIsInitialized = true;    
  }
  */

  // 
  // Define World Volume
  //
  world.solid = new G4Box("World",
			  world.sizeX/2,
			  world.sizeY/2,
			  world.sizeZ/2);

  world.logical = new G4LogicalVolume(world.solid,
				      world.material,
				      "World");
  
  world.physical = new G4PVPlacement(0,
				     G4ThreeVector(),
				     world.logical,
				     "World",
				     0,
				     false,
				     0);

  //
  // Define detector
  //
  sec.solid = new G4Box("Sector",
			sec.sizeX/2.0, 
			sec.sizeY/2.0,
		        sec.sizeZ/2.0);

  sec.logical = new G4LogicalVolume(sec.solid, sec.material, "Sector");

  //pmt SiO2 window1
  pmtWin1.solid = new G4Box("pmtWin1",
			    pmtWin1.sizeX/2.0, 
			    pmtWin1.sizeY/2.0,
			    pmtWin1.sizeZ/2.0);

  pmtWin1.logical = new G4LogicalVolume(pmtWin1.solid, pmtWin1.material, "pmtWin1");

  //pmt SiO2 window2
  pmtWin2.solid = new G4Box("pmtWin2",
			    pmtWin2.sizeX/2.0, 
			    pmtWin2.sizeY/2.0,
			    pmtWin2.sizeZ/2.0);

  pmtWin2.logical = new G4LogicalVolume(pmtWin2.solid, pmtWin2.material, "pmtWin2");



  //
  // fTOF bar parts   ///////////////////// 10.04.18 //////////////////////////////////////////////////////////////////////////
  //

  barBox.solid = new G4Box("barBox",
      barBox.sizeX/2.0, 
      barBox.sizeY/2.0,
      barBox.sizeZ/2.0);

  barBox.logical = new G4LogicalVolume(barBox.solid, barBox.material, "barBox");

  trdXY.solid = new G4Trd("XYtrd",
    trdXY.dX1,
    trdXY.dX2,
    trdXY.dY1,
    trdXY.dY2,
    trdXY.dZ);

  trdXY.logical = new G4LogicalVolume(trdXY.solid, trdXY.material, "XYtrd");

  trdZX.solid = new G4Trd("ZXtrd",
    trdZX.dX1,
    trdZX.dX2,
    trdZX.dY1,
    trdZX.dY2,
    trdZX.dZ);

  trdZX.logical = new G4LogicalVolume(trdZX.solid, trdZX.material, "ZXtrd");

  trdYZ.solid = new G4Trd("YZtrd",
    trdYZ.dX1,
    trdYZ.dX2,
    trdYZ.dY1,
    trdYZ.dY2,
    trdYZ.dZ);

  trdYZ.logical = new G4LogicalVolume(trdYZ.solid, trdYZ.material, "YZtrd");



  //
  // MCP PMT Windows  /////////////////// 11.04.18 ////////////////////////////////////////////////
  //

  hamWin.solid = new G4Box("hamWindow",
      hamWin.sizeX/2.0, 
      hamWin.sizeY/2.0,
      hamWin.sizeZ/2.0);
  hamWin.logical = new G4LogicalVolume(hamWin.solid, 
          hamWin.material,"hamWindow");

  planWin.solid = new G4Box("planWindow",
      planWin.sizeX/2.0, 
      planWin.sizeY/2.0,
      planWin.sizeZ/2.0);
  planWin.logical = new G4LogicalVolume(planWin.solid, 
          planWin.material,"planWindow");

  ///////////////////////////////////////////////////////////////////////////////////////////////////


  //
  // MCP PMT Channels
  //


  hamChan.solid = new G4Box("hamChannel",
      hamChan.sizeX/2.0, 
      hamChan.sizeY/2.0,
      hamChan.sizeZ/2.0);
  hamChan.logical = new G4LogicalVolume(hamChan.solid, 
          hamChan.material,"hamChannel");

  planChan.solid = new G4Box("hamChannel",
      planChan.sizeX/2.0, 
      planChan.sizeY/2.0,
      planChan.sizeZ/2.0);
  planChan.logical = new G4LogicalVolume(planChan.solid, 
          planChan.material,"planChannel");


  //
  // MCP PMT Boxes
  //


  hamBox.solid = new G4Box("hamBox",
      hamBox.sizeX/2.0, 
      hamBox.sizeY/2.0,
      hamBox.sizeZ/2.0);
  hamBox.logical = new G4LogicalVolume(hamBox.solid, 
          hamBox.material,"hamBox");

  planBox.solid = new G4Box("planBox",
      planBox.sizeX/2.0, 
      planBox.sizeY/2.0,
      planBox.sizeZ/2.0);
  planBox.logical = new G4LogicalVolume(planBox.solid, 
          planBox.material,"planBox");



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // fTOF multi union 10.04.18

  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;


  G4VSolid* fullBar = new G4Box("Full_Bar", barBox.sizeX/2.0, barBox.sizeY/2.0, barBox.sizeZ/2.0);
  G4LogicalVolume* fullBarLog = new G4LogicalVolume(fullBar,barBox.material,"BarLogical");

//  Ra = G4RotationMatrix();
//  Ta = G4ThreeVector(0.,0.,0.);
//  Tr = G4Transform3D(Ra,Ta);
//  fullBar->AddNode(*barBox.solid, Tr);

/*
  Ra = G4RotationMatrix();
  Ta.setX(0.);
  Ta.setY(0.);
  Ta.setZ(barBox.sizeZ/2. + trdXY.dZ);
  Tr = G4Transform3D(Ra,Ta);
  fullBar->AddNode(*trdXY.solid, Tr);

  Ra = G4RotationMatrix();
  Ra.rotateX(180.*deg);
  Ta.setX(0.);
  Ta.setY(0.);
  Ta.setZ(-barBox.sizeZ/2. - trdXY.dZ);
  Tr = G4Transform3D(Ra,Ta);
  fullBar->AddNode(*trdXY.solid, Tr);


  Ra = G4RotationMatrix();
  Ra.rotateY(90.*deg);
  Ta.setX(barBox.sizeX/2. + trdXY.dZ);
  Ta.setY(0.);
  Ta.setZ(0.);
  Tr = G4Transform3D(Ra,Ta);
  fullBar->AddNode(*trdYZ.solid, Tr);

  Ra = G4RotationMatrix();
  Ra.rotateY(270.*deg);
  Ta.setX(-barBox.sizeX/2. - trdXY.dZ);
  Ta.setY(0.);
  Ta.setZ(0.);
  Tr = G4Transform3D(Ra,Ta);
  fullBar->AddNode(*trdYZ.solid, Tr);


  Ra = G4RotationMatrix();
  Ra.rotateX(270.*deg);
  Ta.setX(0.);
  Ta.setY(barBox.sizeY/2. + trdXY.dZ);
  Ta.setZ(0.);
  Tr = G4Transform3D(Ra,Ta);
  fullBar->AddNode(*trdZX.solid, Tr);

  Ra = G4RotationMatrix();
  Ra.rotateX(90.*deg);
  Ta.setX(0.);
  Ta.setY(-barBox.sizeY/2. - trdXY.dZ);
  Ta.setZ(0.);
  Tr = G4Transform3D(Ra,Ta);
  fullBar->AddNode(*trdZX.solid, Tr);
*/
 // fullBar->Voxelize();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





  //
  // Sensitive volume (Photocathode)
  //
  sensitive.solid = new G4Box("Sensitive",
			      sensitive.sizeX/2.0,
			      sensitive.sizeY/2.0,
			      sensitive.sizeZ/2.0);
  sensitive.logical = new G4LogicalVolume(sensitive.solid, 
					  sensitive.material,"Sensitive");

  //
  // Pmt Box
  //  
  pmtbox.solid = new G4Box("PmtBox",
			   pmtbox.sizeX/2.0,
			   pmtbox.sizeY/2.0,
			   pmtbox.sizeZ/2.0);
  pmtbox.logical = new G4LogicalVolume(pmtbox.solid, 
				       pmtbox.material,"PmtBox");
  
  //
  //Pmt Abs1
  //
  pmtAbs1.solid = new G4Box("pmtAbs1",
			    pmtAbs1.sizeX/2.0,
			    pmtAbs1.sizeY/2.0,
			    pmtAbs1.sizeZ/2.0);
  pmtAbs1.logical = new G4LogicalVolume(pmtAbs1.solid, 
					pmtAbs1.material,"pmtAbs1");

  //
  // abs 1;
  //
  abs1.solid = new G4Box("abs1",
			 abs1.sizeX/2.0,
			 abs1.sizeY/2.0,
			 abs1.sizeZ/2.0);
  abs1.logical = new G4LogicalVolume(abs1.solid, 
				     abs1.material,"abs1");

  //
  // abs 2;
  //
  abs2.solid = new G4Box("abs2",
			 abs2.sizeX/2.0,
			 abs2.sizeY/2.0,
			 abs2.sizeZ/2.0);
  abs2.logical = new G4LogicalVolume(abs2.solid, 
				     abs2.material,"abs2");
  
  
  //-------------------------------------------------------

  G4AssemblyVolume* secAssembly = new G4AssemblyVolume();

  //--------------------------------------------------------



  // Lazy include fTOF bar 10.04.18 ///////////////////////////////////////////

  Ta = G4ThreeVector(0.,0.,0.);
  Ra = G4RotationMatrix(); 
  Tr = G4Transform3D(Ra,Ta);
  secAssembly->AddPlacedVolume(fullBarLog,Tr);

  ///////////////////////////////////////////////////////////////////////////////////////

  // mcp pmt windows 11.04.18 /////////////////////////////////////////////////

//  Ta = G4ThreeVector(barBox.sizeX/2.+2*trdYZ.dZ+hamWin.sizeX/2.,
//    0.,
//    0.);
//  Ra = G4RotationMatrix();
//  Tr = G4Transform3D(Ra,Ta);
//  secAssembly->AddPlacedVolume(hamWin.logical,Tr);

  // Ta = G4ThreeVector(-barBox.sizeX/2.-2*trdYZ.dZ-planWin.sizeX/2.,  // nice and precize
  //   (fTOFConst::hamDownDist - fTOFConst::hamUpDist)/2.,
  //   (fTOFConst::hamLeftDist - fTOFConst::hamRightDist)/2.);
  // Ra = G4RotationMatrix(); 
  // Tr = G4Transform3D(Ra,Ta);
  // secAssembly->AddPlacedVolume(planWin.logical,Tr);

//  Ta = G4ThreeVector(-barBox.sizeX/2.-2*trdYZ.dZ-planWin.sizeX/2.,
//    0.,
//    0.);
//  Ra = G4RotationMatrix();
//  Tr = G4Transform3D(Ra,Ta);
//  secAssembly->AddPlacedVolume(planWin.logical,Tr);



  // MCP boxes

//  Ta = G4ThreeVector(barBox.sizeX/2. + 2*trdYZ.dZ + hamWin.sizeX + hamBox.sizeX/2. + hamChan.sizeX,
//    0.,
//    0.);
//  Ra = G4RotationMatrix();
//  Tr = G4Transform3D(Ra,Ta);
//  secAssembly->AddPlacedVolume(hamBox.logical,Tr);

//  Ta = G4ThreeVector(-barBox.sizeX/2. - 2*trdYZ.dZ - planWin.sizeX - planBox.sizeX/2. - planChan.sizeX,
//    0.,
//    0.);
//  Ra = G4RotationMatrix();
//  Tr = G4Transform3D(Ra,Ta);
//  secAssembly->AddPlacedVolume(planBox.logical,Tr);

  //

  G4int i = 0;
  Ra = G4RotationMatrix(); 
                      ///////////////////// hamamatsu 11.04.18 //////////////////////////////////////////////////////////////////////
    //add sensitive

    Ta.setX(barBox.sizeX/2.0 + hamChan.sizeX/2.);
    Ta.setY(0);
    Ta.setZ(0);
    Tr = G4Transform3D(Ra, Ta);
    secAssembly->AddPlacedVolume(hamChan.logical, Tr);

    //add sensitive
    Ta.setX(-(barBox.sizeX/2.0 + hamChan.sizeX/2.));
    Ta.setY(0);
    Ta.setZ(0);
    Tr = G4Transform3D(Ra, Ta);
    secAssembly->AddPlacedVolume(hamChan.logical, Tr);

//                       ///////////////////// planacon 11.04.18 //////////////////////////////////////////////////////////////////////
//    //add sensitive
//    Ta.setX(-(barBox.sizeX/2.0 + trdYZ.dZ*2. + hamWin.sizeX + hamChan.sizeX/2.));
//    Ta.setY(planChan.sizeY/2.);
//    Ta.setZ( 1.5*planChan.sizeZ - (planChan.sizeZ)*i);
//    Tr = G4Transform3D(Ra, Ta);
//    secAssembly->AddPlacedVolume(planChan.logical, Tr);


                                                                           ////////////////// 10.04.18

  
  //
  //make Imprint
  //
  
  Ra.rotateX(90.0*deg);
  //Ra.rotateY(20.0*deg);
  //Ra.rotateY(180.0*deg);

  //One
  Ta.setX(0.0);
  Ta.setY(0.0);
  Ta.setZ(0.0);
  Tr = G4Transform3D(Ra, Ta);
  secAssembly->MakeImprint(world.logical, Tr, 0, true);

  //Ta.setX(0.0);
  //Ta.setY(-5.0*cm);
  //Ta.setZ(0.0*cm);
  //Tr = G4Transform3D(Ra, Ta);
  //secAssembly->MakeImprint(world.logical, Tr, 0, true);



  //-----------------------------------------------------

  //
  // Set Visualization Attributes
  //
  G4Color blue        = G4Color(0., 0., 1.);
  G4Color green       = G4Color(0., 1., 0.);
  G4Color red         = G4Color(1., 0., 0.);
  G4Color white       = G4Color(1., 1., 1.);
  G4Color cyan        = G4Color(0., 1., 1.);
  G4Color DircColor   = G4Color(0.0, 0.0, 1.0, 0.2);
  G4Color SensColor   = G4Color(0.0, 1.0, 1.0, 0.1);

  worldVisAtt->SetColor(white);
  worldVisAtt->SetVisibility(true);
  quartzVisAtt->SetColor(DircColor);
  quartzVisAtt->SetVisibility(true);
  sensitiveVisAtt->SetColor(SensColor);
  sensitiveVisAtt->SetVisibility(true);
  pmtboxVisAtt->SetColor(red);
  pmtboxVisAtt->SetVisibility(true);
  absVisAtt->SetColor(red);
  absVisAtt->SetVisibility(true);


  // fTOF vis attributes 11.04.18 //////////////////

  fullBarLog->SetVisAttributes(quartzVisAtt);
  hamWin.logical->SetVisAttributes(quartzVisAtt);
  planWin.logical->SetVisAttributes(quartzVisAtt);

  hamChan.logical->SetVisAttributes(sensitiveVisAtt);
  planChan.logical->SetVisAttributes(sensitiveVisAtt);

  hamBox.logical->SetVisAttributes(pmtboxVisAtt);
  planBox.logical->SetVisAttributes(pmtboxVisAtt);

  //////////////////////////////////////////////////




  world.logical->SetVisAttributes(worldVisAtt);

  sec.logical->SetVisAttributes(quartzVisAtt);
  pmtWin1.logical->SetVisAttributes(quartzVisAtt);
  pmtWin2.logical->SetVisAttributes(quartzVisAtt);

  sensitive.logical->SetVisAttributes(sensitiveVisAtt);

  pmtbox.logical->SetVisAttributes(pmtboxVisAtt);

  pmtAbs1.logical->SetVisAttributes(absVisAtt);
  abs1.logical->SetVisAttributes(absVisAtt);
  abs2.logical->SetVisAttributes(absVisAtt);


  //
  // Define Optical Borders
  //

  // Surface for killing photons at borders
  const G4int num1 = 2;
  G4double Ephoton[num1] = {1.5*eV, 5.8*eV};

  G4OpticalSurface* OpVolumeKillSurface =
    new G4OpticalSurface("VolumeKillSurface");
  OpVolumeKillSurface->SetType(dielectric_metal);
  OpVolumeKillSurface->SetFinish(polished);
  OpVolumeKillSurface->SetModel(glisur);

  G4double ReflectivityKill[num1] = {0., 0.};
  G4double EfficiencyKill[num1] = {1., 1.};
  G4MaterialPropertiesTable* VolumeKill = new G4MaterialPropertiesTable();
  VolumeKill->AddProperty("REFLECTIVITY", Ephoton, ReflectivityKill, num1);
  VolumeKill->AddProperty("EFFICIENCY",   Ephoton, EfficiencyKill,   num1);
  OpVolumeKillSurface->SetMaterialPropertiesTable(VolumeKill);
  new G4LogicalSkinSurface("SensitiveSurface", 
			   sensitive.logical, OpVolumeKillSurface);
  
  
  // 
  // Sensitive detector definition
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  fTOF_SensitiveDetector* aSD = new fTOF_SensitiveDetector("fTOF");
  SDman->AddNewDetector(aSD);
  sensitive.logical->SetSensitiveDetector(aSD);
  


  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  printDetectorParameters();


  return world.physical;
}

void fTOF_DetectorConstruction::printDetectorParameters(){
  

}
