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
    /*Set energy range from 0 to 310MeV*/
    setHistoRanges(0, 310 * 1e03, 10000);
}

DetectorConstruction::~DetectorConstruction() 
{

}

/**
   Create a new object of DetectorSD class and record it's
   pointer value to vector of pointers.

   \param name of the sensitive region


*/
DetectorSD2 * DetectorConstruction::newRegisteredDetectorSD2(const G4String name)
{
    //add a new detector to the container
    pbDetectorSDArray.push_back(new DetectorSD2(name));
    //now pull out the pointer:
    std::vector<DetectorSD2*>::iterator iter =
            pbDetectorSDArray.end()-1;
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
void DetectorConstruction::readParameters(const std::map< G4String, G4double > &str_double_map)
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
void DetectorConstruction::setHistoRanges(const double min,
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
DetectorSD2 * DetectorConstruction::newCylindricDetector(G4LogicalVolume *worldLogical,
                                                         G4String detName,
                                                         G4Material * detMaterialPtr,
                                                         G4ThreeVector centerPlacement,
                                                         G4double cylinderDiameter,
                                                         G4double cylinderHeight,
                                                         G4RotationMatrix* pRotationMatrix)
{
    /** detector logical volume.*/
    G4LogicalVolume *detectorLogicalPointer;
    g4solid_object<G4Tubs> *detectorCylinder;

    detectorCylinder =   make_cylinder(worldLogical,
                                       detName + G4String("_cylinder"),
                                       detMaterialPtr,  centerPlacement,
                                       cylinderDiameter, cylinderHeight, 0,
                                       pRotationMatrix);

    DetectorSD2 *sd2Pointer = newRegisteredDetectorSD2(detName);
    detectorLogicalPointer = detectorCylinder->get_logical();
    detectorLogicalPointer->SetSensitiveDetector(sd2Pointer);
    return sd2Pointer;
}
//-------------------------------------------------------------------------------
/** This macros adds new sensitive detector.
 *  @param detName: (G4String) detector's name.
 *  @param detMaterialPtr: (G4Material*) pointer to detector's material.
 *  @param centerPlacement: (G4ThreeVector) radius vector of the detector's center.
 *	@param boxHalfDimensions: vector of half sizes of the box (x/2, y/2, z/)
 *	@param pRotationMatrix: (optional) pointer to rotation matrix,
 * the object shall be rotated if not NULL.
 *	@return pointer to newly created DetectorSD2 object.
**/
DetectorSD2 * DetectorConstruction::newBoxDetector(G4LogicalVolume *worldLogical,
                                                   G4String detName,
                                                   G4Material * detMaterialPtr,
                                                   G4ThreeVector centerPlacement,
                                                   G4ThreeVector boxHalfDimensions,
                                                   G4RotationMatrix* pRotationMatrix)
{
    /** detector logical volume.*/
    G4LogicalVolume *detectorLogicalPointer;
    g4solid_object<G4Box> *detectorBox;

    detectorBox =   make_box(worldLogical,
                             detName + G4String("_box"),
                             detMaterialPtr,
                             centerPlacement,
                             boxHalfDimensions.x(),
                             boxHalfDimensions.y(),
                             boxHalfDimensions.z(),
                             pRotationMatrix);

    DetectorSD2 *sd2Pointer = newRegisteredDetectorSD2(detName);
    detectorLogicalPointer = detectorBox->get_logical();
    detectorLogicalPointer->SetSensitiveDetector(sd2Pointer);
    return sd2Pointer;
}
//-------------------------------------------------------------------------------
void DetectorConstruction::setHistoBinsQuantity(const G4int bins)
{
    d_hist_bins = bins;
}
void DetectorConstruction::setHistoEnergyMinimum(const G4double min)
{
    d_hist_min = min;
}
void DetectorConstruction::setHistoEnergyMaximum(const G4double max)
{
    d_hist_max = max;
}
unsigned DetectorConstruction::getHistoEnergyUnits()
{
    return d_energy_units;
}
//-------------------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // --- materials ---
    // создаем материалы
    // первый способ:
    G4Element* Tl = new G4Element("Tantal", "Tl", 73, 180.9479*g/mole);

    G4Material* TargetMaterial = new G4Material("TargetMaterial", 16.65*g/cm3, 1);
    TargetMaterial->AddElement(Tl, 100*perCent);


    // второй способ: использует встроенную в Geant4 базу материалов
    // более простой, но иногда приходится прибегать к первому способу,
    // т.к. не все материалы содержатся в базе
    G4NistManager* nistMan = G4NistManager::Instance();
    G4Material* Air = nistMan->FindOrBuildMaterial("G4_AIR");
    G4Material* saMaterial = nistMan->FindOrBuildMaterial("G4_Al");
    //G4Material* detMaterial = nistMan->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    G4Element *H = new G4Element("Hydrogen","H",1,1*g/mole);

    G4Material *void_dumb_material =  new G4Material("dumb_void", 1e-10*g/cm3, 1);
    void_dumb_material->AddElement(H,1);

    // --- volumes ---
    // создаем геометрию
    G4double saSize = 3*cm; // размер образца
    G4double saThick = 0.5*mm; // толщина образца
    G4double detDiam = 3*cm; // диаметр детектора
    G4double detThick = 1*mm; // толщина детектора
    G4double gap1 = 5*cm; // расстояние от источника до образца
    G4double gap2 = 10*cm; // расстояние от образца до детектора

    // мировой объем в виде параллелепипеда
    G4Box* world_box = new G4Box("world", (saSize + detDiam)/2 + 1*cm, (saSize + detDiam)/2 + 1*cm, gap1 + gap2 + detThick/2 + 6*cm);
    // заполняем его воздухом
    G4LogicalVolume* world_log = new G4LogicalVolume(world_box, Air, "world");
    // и помещаем в начало координат
    G4VPhysicalVolume* world_physical_volume = new G4PVPlacement(0, G4ThreeVector(), world_log, "world", 0, false, 0);

    g4solid_object<G4Box>* sample_box_p =
            make_box(world_log, "sample", TargetMaterial,
                    /*placement*/ G4ThreeVector(0, 0, gap1 + saThick / 2.0),
                   /*half sizes*/G4ThreeVector(saSize/2.0, saSize/2.0, saThick/2.0));

    g4solid_object<G4Tubs>* al_tube_p = make_cylinder(world_log, "poseredini",
                                                      saMaterial,
                                                      G4ThreeVector(0, 0, gap1 + saThick + 5*cm),
                                                      /*outer radius*/ 2 * cm,
                                                      /*height*/ 5 * cm);

    //---------------------------------------------
    // --- visualisation ---
    // abs_log->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
    world_log->SetVisAttributes(G4VisAttributes::Invisible);

    //get pointer to sensitive detector manager:
    G4SDManager *det_manager = G4SDManager::GetSDMpointer();

    /* Detector inside the box.
     * This detector shall be registered in vector: pbDetectorSDArray.*/
    DetectorSD2  *sensDetectorPtr =
            newCylindricDetector(world_log,
                                 "MYDETECTOR",
                                 /*void material, this detector is virtual:*/void_dumb_material,
                                 /*placement:*/G4ThreeVector(0, 0, gap1 + saThick + gap2 + detThick / 2.0),
                                 /*cylinder diameter*/detDiam,
                                 /*cylinder height*/  detThick);

    /*Register new detector in DetectorManager:*/
    det_manager->AddNewDetector(sensDetectorPtr);
    /* Make detector virtual (counts only kinetic energy, no deposited energy)*/
    sensDetectorPtr->DisableDepositedEnergyCount();

    return world_physical_volume;
}

