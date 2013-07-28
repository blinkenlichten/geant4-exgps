// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
///*
// Part of simulation for use with GEANT4 code.
// Author  Bogdan Maslovskiy <blinkofnight(at,doggy,removeme)gmail(dot)com>
// Taras Schevchenko National University of Kyiv 2012
//****************

#include "DetectorConstruction.hh"
#include "quick_geom.hh"

DetectorConstruction::DetectorConstruction()
{

}

DetectorConstruction::~DetectorConstruction() 
{
  
}

/**
   Create a new object of DetectorSD class and record it's
   pointer value to vector of pointers.
     
   \param name of the sensitive region
     

*/
DetectorSD2 * DetectorConstruction::new_detector_sensitive(const G4String name)
{
  //add a new detector to the container
  vector_DetectorSD.push_back(new DetectorSD2(name));
  //now pull out the pointer:
  std::vector<DetectorSD2*>::iterator iter = 
    vector_DetectorSD.end()-1;
  (*iter)->fill_hist("e-", 0);
  (*iter)->fill_hist("e+", 0);
  (*iter)->fill_hist("gamma", 0);
  (*iter)->fill_hist_deposited("e-", 0);
  (*iter)->fill_hist_deposited("e+", 0);
  (*iter)->fill_hist_deposited("gamma", 0);
  return (*iter);
}

/** 
    Read parameters values from a map<G4String, G4double>;
    Possible G4Sting keys are:
    "tathick", "tadiam", "tgthick", "tgdiam", "tgmass", "althick", "aldiam".
    example1:
    If map contains a pair of values like: "tgmass" and 0.200*g,
    then target detector mass will be set to 0.2g.
    
    //one may actually use this code :
    str_double_map.insert( std::pair<G4String, G4double>
    ("tgmass", 0.200*g));
    //to insert that values into map
    
    example2:    
    If map contains a pair of values: "tathick" and 0.100*cm,
    then Ta-plate thickness will be set to 0.1*cm.
    
    \param std::map reference(a container with G4String and G4double pairs).

*/
void DetectorConstruction::read_parameters(const std::map< G4String, G4double > &str_double_map)
{
  std::map< G4String, G4double >::const_iterator str_double_iterator;

  //add strings to the option's list:
  for( str_double_iterator = str_double_map.begin();
       str_double_iterator != str_double_map.end(); 
       str_double_iterator++)  {
      G4double value = str_double_iterator->second;
      G4String key = str_double_iterator->first;
      // switch(key) {
      // case "some_key": {      }
      // default:
      // }
  }

}

/** Set histogram properties.
    \param minimum value of range. 
    Anything lesser than this setpoint will be ignored.
      
    \param maximum value of range. 
    Anything greater than this setpoint will be ignored.
      
    \param Quantity of histogram bins.
*/
void DetectorConstruction::set_histo(const double min,
				     const double max,
				     const unsigned int n_bins, 
				     const unsigned int E_units)
{
  //lets build it up:
  
  d_hist_min = min;
  d_hist_max = max;
  if(d_hist_max < d_hist_min)//if somehow on Earth..
    {
      d_hist_max = min;  d_hist_min = max;
    }
  d_hist_bins = n_bins;
  d_energy_units = E_units;
}
//-----------------------------------------------------------------------------

/** This macros adds new sensitive detector.
 *  @param detName: (G4String) detector's name.
 *  @param detMaterialPtr: (G4Material*) pointer to detector's material.
 *  @param centerPlacement: (G4ThreeVector) radius vector of the detector's center.
 *  @param cylinderDiameter: (float)diameter of the detector's cylinder.
 *  @param cylinderHeight: (float) height of the detector's cylinder.
**/
DetectorSD2 * DetectorConstruction::CreateNewDetector(G4LogicalVolume *worldLogical,
                                                      G4String detName,
                                                      G4Material * detMaterialPtr,
                                                      G4ThreeVector centerPlacement,
                                                      G4double cylinderDiameter,
                                                      G4double cylinderHeight)
{
  /** detector logical volume.*/
  G4LogicalVolume *detectorLogicalPointer;
  g4solid_object<G4Tubs> *detectorCylinder;

  detectorCylinder =   make_cylinder(worldLogical,
                                     detName + G4String("_cylinder"),
                                     detMaterialPtr,  centerPlacement, cylinderDiameter, cylinderHeight);

  DetectorSD2 *sd2Pointer = new_detector_sensitive(detName);
  detectorLogicalPointer = detectorCylinder->get_logical();
  detectorLogicalPointer->SetSensitiveDetector(sd2Pointer);
  return sd2Pointer;
  }
//-------------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::Construct()
{
 
  //get pointer to sensitive detector manager:
  G4SDManager *det_manager = G4SDManager::GetSDMpointer();

  // --- materials ---
  // первший спосіб, за хімічними елементами:
  G4Element *N    = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element *O    = new G4Element("Oxygen"  , "O", 8, 16.00*g/mole);
  G4Material *Air = new G4Material("Air", 1.29*mg/cm3, 2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);
  
  G4Element *H = new G4Element("Hydrogen","H",1,1*g/mole);
  G4Material *void_dumb_material =  new G4Material("dumb_void", 1e-10*g/cm3, 1);
  void_dumb_material->AddElement(H,1);
  
  // другий спосіб: використовуючи базу готових матеріалів.
  G4NistManager *nistMan = G4NistManager::Instance();
  G4Material *Fe_material= nistMan->FindOrBuildMaterial("G4_Fe");
  G4Material *Poly_material= nistMan->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material *concreteMaterial = nistMan->FindOrBuildMaterial("G4_CONCRETE");
  G4Material *detector_material = nistMan->FindOrBuildMaterial("G4_In");
  
  //---------------------------------------------------------------
  // робимо геометрію
  //---------------------------------------------------------------

  //-------- world box-------------------------
  // світ -- він і є "світ", тут будуть розташовані всі об'єкти:
  G4Box *world_box = new G4Box("WORLD_BOX",
			       //3 параметри паралелепіпеда X,Y,Z
			       20.15*m,
			       20.15*m,
			       20.15*m);

  // заповним повітрям лабораторию:
  G4LogicalVolume *world_logical_volume =
    new G4LogicalVolume(world_box, Air,	"WORLD_LOG");
  
  G4VPhysicalVolume *world_physical_volume = 
    new G4PVPlacement( 0, 
		       G4ThreeVector(),
		       world_logical_volume,
		       "WORLD_PHYS",
		       0, false, 0);
  //------------ Objects ---------------------------------------------------
   g4solid_object<G4Box> *plita =
       make_box(world_logical_volume,
               "plita",
               concreteMaterial,
               /*розташування, радіус-вектор до центру:*/G4ThreeVector(0,0,-4.65*m),
               /*3D-розміри ПОЛОВИННІ коробки: (X/2,Y/2,Z/2) */1.5 *m, 1.5 *m, 0.1 *m);
 //z=-4.65*m

   g4solid_object<G4Box> *febox =
       make_box(world_logical_volume,
               "febox",
               Fe_material,
               G4ThreeVector(0,0, -10.15 *m), 1 *m, 1 *m, 1 *m);
  

  G4ThreeVector polyboxCenter = G4ThreeVector(0,0, -10.15 *m);
  /** This function makes a box with a hole inside (shaped as box too).*/
  make_box_with_hole( world_logical_volume,
		      "polybox",
		      Poly_material,
		      polyboxCenter /*box center*/,
		      /*ЦЕ САМОПАЛЬНА Ф-ЦІЯ, ТУТ ПОВНІ РОЗМІРИ КОРОБОК(НЕ ПОЛОВИННІ)*/
		      G4ThreeVector(900*cm, 900*cm, 900*cm) /*box dimensions(width, height, depth)*/,
		      G4ThreeVector(170*cm, 170*cm, 170*cm)/*hole dimensions(width, height, depth)*/,
		      NULL, NULL);

  /* Detector inside the box. */
  DetectorSD2  *sensDetectorPtr =
      CreateNewDetector(world_logical_volume,
                        "DET.INSIDE",
                        /*void material, this detector is virtual:*/void_dumb_material,
                        /*placement:*/G4ThreeVector(0,0,-10.15*m),
                        /*cylinder diameter*/1.3*m,
                        /*cylinder height*/  10*cm);

  /*Register new detector in DetectorManager:*/
  det_manager->AddNewDetector(sensDetectorPtr);
  /*Register this 'virtual' detector in array, this will let SteppingAction class
    to store information about particles to known detectors.*/
  vector_DetectorSD.push_back(sensDetectorPtr);
  /* Make detector virtual (counts only kinetic energy, no deposited energy)*/
  sensDetectorPtr->DisableDepositedEnergyCount();

  return world_physical_volume;
}

