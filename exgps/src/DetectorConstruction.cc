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
                                                         G4double cylinderHeight)
{
    /** detector logical volume.*/
    G4LogicalVolume *detectorLogicalPointer;
    g4solid_object<G4Tubs> *detectorCylinder;

    detectorCylinder =   make_cylinder(worldLogical,
                                       detName + G4String("_cylinder"),
                                       detMaterialPtr,  centerPlacement, cylinderDiameter, cylinderHeight);

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
                                                   G4RotationMatrix* pRotationMatrix = 0
                                                   )
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

    //Матеріали
    G4NistManager* nistMan = G4NistManager::Instance();
    G4Material* vacuum = nistMan->FindOrBuildMaterial("G4_Galactic");
    // Все в боксі. Бокс має повітря. Шестикутна призма - графіт.
    // G4NistManager* nistMan = G4NistManager::Instance();
    G4Material* Air = nistMan->FindOrBuildMaterial("G4_AIR");
    G4Material* saMaterial = nistMan->FindOrBuildMaterial("G4_GRAPHITE");
    // Детектор - нержавіюча сталь AISI 304 (взято з Вікі)
    G4Element* Fe = new G4Element("Ferum"  , "Fe", 26, 55.847*g/mole);
    G4Element* Cr = new G4Element("Chrom"  , "Cr", 24, 51.9961*g/mole);
    G4Element* Ni = new G4Element("Nikel"  , "Ni", 28, 58.6934*g/mole);
    G4Element* C = new G4Element("Carbon"  , "C", 6, 12.011*g/mole);
    G4Element* Mn = new G4Element("Mangan"  , "Mn", 25, 54.93805*g/mole);
    G4Element* P = new G4Element("Phosphor"  , "P", 15, 30.9737*g/mole);
    G4Element* S = new G4Element("Sirka"  , "S", 16, 32.066*g/mole);
    G4Element* Cu = new G4Element("Cuprum"  , "Cu", 29, 63.546*g/mole);
    G4Material* UnrustIron = new G4Material("UnrustIron", 8*g/cm3, 8);
    UnrustIron->AddElement(Fe, 70*perCent);
    UnrustIron->AddElement(Cr, 18*perCent);
    UnrustIron->AddElement(Ni, 9*perCent);
    UnrustIron->AddElement(C, 0.08*perCent);
    UnrustIron->AddElement(Mn, 2*perCent);
    UnrustIron->AddElement(P, 0.045*perCent);
    UnrustIron->AddElement(S, 0.03*perCent);
    UnrustIron->AddElement(Cu, 0.845*perCent);

    //---------------------------------------------------------------
    // робимо геометрію
    //---------------------------------------------------------------

    //-------- world box-------------------------
    // світ -- він і є "світ", тут будуть розташовані всі об'єкти:
    G4Box *world_box = new G4Box("WORLD_BOX",
                                 //3 параметри паралелепіпеда X,Y,Z
                                 2 * m,
                                 2 * m,
                                 2 * m);

    // заповним повітрям лабораторию:
    G4LogicalVolume *world_logical_volume =
            new G4LogicalVolume(world_box, Air,	"WORLD_LOG");

    G4VPhysicalVolume *world_physical_volume =
            new G4PVPlacement( 0,
                               G4ThreeVector(),
                               world_logical_volume,
                               "WORLD_PHYS",
                               0, false, 0);
    //---------------------------------------------
    G4double zplane[]={-0.25*m, 0.25*m};
    G4double rin[]={0*m,0*m};
    G4double rout[]={135*mm,135*mm};
    G4Polyhedra* sample_box = new G4Polyhedra("sample", 0 , 2*pi,6,2, zplane,rin,rout);
    // образец
    G4LogicalVolume* sample_log = new G4LogicalVolume(sample_box, saMaterial, "sample");
    // помещаем его в мировой объем
    G4VPhysicalVolume* sample_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), sample_log, "sample", world_logical_volume, false, 0);

    //UnrustIron- material of the sample

    // детектор в виде цилиндра
    //  G4Box* det_tube = new G4Box("detector", 0.25*m, 0.1*m, 0.25*m);
    //G4Tubs* det_tube = new G4Tubs("detector", 0, detDiam/2, detThick/2, 0, 360*deg);
    // G4LogicalVolume* det_log = new G4LogicalVolume(det_tube, vacuum, "detector");
    // помещаем его в мировой объем со смещением на gap1+gap2
    // G4VPhysicalVolume* det_phys = new G4PVPlacement(0, G4ThreeVector(0, 0.5*m, 0), det_log, "detector", world_logical_volume, false, 0);


    // калориметр
    G4Box* calor_box = new G4Box("calorimeter", 1*m,1*m,25*cm);
    G4LogicalVolume* calor_log = new G4LogicalVolume(calor_box, Air, "calorimeter");
    G4VPhysicalVolume* calor_phys =
            new G4PVPlacement(0, G4ThreeVector(0, gap + calorThick/2, 0),
                              calor_log,
                              "calorimeter",
                              world_logical_volume, false, 0);

    // слой
    G4Box* layer_box = new G4Box("layer", 25*cm, 5*cm, layerThick/2);
    G4LogicalVolume* layer_log = new G4LogicalVolume(layer_box, vacuum, "layer");

    // помещаем внутрь калориметра nLayers слоев (копий объекта layer_log)
    // задается направление вдоль которого будут размещены слои (kZAxis) и толщина слоя
    G4PVReplica* layer_phys = new G4PVReplica("layer", layer_log, calor_log, kZAxis, nLayers, layerThick);
    G4VPhysicalVolume *det_phys[50][10];
    // детектор
    G4LogicalVolume* det_log;
    if (detThick > 0.) {
        G4Box* det_box = new G4Box("detector", size/2, size/2, detThick/2);
        det_log = new G4LogicalVolume(det_box, UnrustIron, "detector");
        // входит в соfor став слоя (abs_log)
        for(int j=0; j<10; j++)
        { for(int i=0;i<50;i++)
                det_phys[i][j] =
                        new G4PVPlacement(0,
                                          G4ThreeVector( i*size/2-12*size, j*size/2-3*size , 0),
                                          det_log, "detector", layer_log, false, 0);
        }
    }

    //---------------------------------------------
    // --- visualisation ---
    // отключаем отображение мирового объема

    // abs_log->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
    layer_log->SetVisAttributes(G4VisAttributes::Invisible);

    world_log->SetVisAttributes(G4VisAttributes::Invisible);
    calor_log->SetVisAttributes(G4VisAttributes::Invisible);
    sample_log->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));
    //---------------------------------------------

    //get pointer to sensitive detector manager:
    G4SDManager *det_manager = G4SDManager::GetSDMpointer();

    /* Detector inside the box. */
    DetectorSD2  *sensDetectorPtr =
            newCylindricDetector(world_logical_volume,
                                 "DET.INSIDE",
                                 /*void material, this detector is virtual:*/void_dumb_material,
                                 /*placement:*/G4ThreeVector(0,0,-10.15*m),
                                 /*cylinder diameter*/1.3*m,
                                 /*cylinder height*/  10*cm);

    /*Register new detector in DetectorManager:*/
    det_manager->AddNewDetector(sensDetectorPtr);
    /*Register this 'virtual' detector in array, this will let SteppingAction class
    to store information about particles to known detectors.*/
    pbDetectorSDArray.push_back(sensDetectorPtr);
    /* Make detector virtual (counts only kinetic energy, no deposited energy)*/
    sensDetectorPtr->DisableDepositedEnergyCount();

    return world_physical_volume;
}

