#ifndef fTOF_DetectorConstruction_H
#define fTOF_DetectorConstruction_H 1

//My
#include "fTOF_VolumeStructures.hh"

//G4
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"

class MagneticField;

class fTOF_DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  fTOF_DetectorConstruction();
  ~fTOF_DetectorConstruction();
   
public:
  
  G4VPhysicalVolume* Construct();
  void ConstructDetector();

private:
  void DefineMaterials();

public:

  void printDetectorParameters();

private:

  //MagneticField *magField; 

  // Various visibility attributes
  G4VisAttributes* worldVisAtt;
  G4VisAttributes* quartzVisAtt;
  G4VisAttributes* sensitiveVisAtt;
  G4VisAttributes* pmtboxVisAtt;
  G4VisAttributes* absVisAtt;

  WorldStruct world;
  SecStruct sec;
  PmtWindow1Struct pmtWin1;
  PmtWindow2Struct pmtWin2;
  SensitiveStruct sensitive;
  PmtBoxStruct pmtbox;
  PmtAbs1Struct pmtAbs1;
  Abs1Struct abs1;
  Abs2Struct abs2;


  // fTOF bar 10.04.18
  BarBoxStruct barBox;
  BarTrdXYStruct trdXY;
  BarTrdYZStruct trdYZ;
  BarTrdZXStruct trdZX;


 // mcp pmt's 11.04.18
  HamPmtWindowStruct hamWin;
  PlanPmtWindowStruct planWin;

  HamPmtChannelStruct hamChan;
  PlanPmtChannelStruct planChan;

  HamPmtBoxStruct hamBox;
  PlanPmtBoxStruct planBox;

  //LB need to be don stability tests
  //G4UserLimits* stepLimit;  // pointer to user step limits

};

#endif

