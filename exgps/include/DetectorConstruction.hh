/* ========================================================== */
// Code was written by Anton Korneev and Mechinsky Vitaly, 
// Institute of Nuclear Problems, Belarus, Minsk, September 2007
//
// And rewritten by Bogdan Maslovskiy,
// Taras Schevchenko National University of Kyiv, 2011.
/* ========================================================== */

#ifndef DetectorConstruction_H
#define DetectorConstruction_H

#include "DetectorSD2.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "globals.hh"
#include "G4VisAttributes.hh" 
#include "G4SDManager.hh"

#include <vector>
#include "G4VUserDetectorConstruction.hh"
#include "DetectorSD.hh"
#include "DetectorConstructionMessenger.hh"

#include <string>
#include "geom_objects.h"
#include <math.h>

#ifndef M_PI
        #define M_PI 3.141592653589793
#endif

class DetectorConstruction : public G4VUserDetectorConstruction
{                 
  public:
    DetectorConstruction();
    ~DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  
  /** Vector of pointers to DetectorSD objects.
      It will be filled on this->Construct() method call.
      Address of this object may be passed to other places.
   */
  std::vector<DetectorSD2*> pbDetectorSDArray;


  /** 
      Read parameters values from a map<G4String, G4double>;
      Possible G4Sting keys are:
      (empty)
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
  void readParameters(const std::map<G4String, G4double> &str_double_map);
  

protected:
  /**
     Create a new object of DetectorSD class and record it's
     pointer value to vector of pointers.
     
     \param name of the sensitive region
     
  */
  DetectorSD2 * newRegisteredDetectorSD2(const G4String name);

  /** This macros adds new sensitive detector.
   *  @param detName: (G4String) detector's name.
   *  @param detMaterialPtr: (G4Material*) pointer to detector's material.
   *  @param centerPlacement: (G4ThreeVector) radius vector of the detector's center.
   *  @param cylinderDiameter: (float)diameter of the detector's cylinder.
   *  @param cylinderHeight: (float) height of the detector's cylinder.
  **/
  DetectorSD2 * newCylindricDetector(G4LogicalVolume *worldLogical,
                                  G4String detName,
                                  G4Material * detMaterialPtr,
                                  G4ThreeVector centerPlacement,
                                  G4double cylinderDiameter,
                                  G4double cylinderHeight,
                                  G4RotationMatrix* pRotationMatrix = 0);
  /** This macros adds new sensitive detector.
 *  @param detName: (G4String) detector's name.
 *  @param detMaterialPtr: (G4Material*) pointer to detector's material.
 *  @param centerPlacement: (G4ThreeVector) radius vector of the detector's center.
 *	@param boxHalfDimensions: vector of half sizes of the box (x/2, y/2, z/)
 *	@param pRotationMatrix: (optional) pointer to rotation matrix,
 * the object shall be rotated if not NULL.
 *	@return pointer to newly created DetectorSD2 object.
**/
  DetectorSD2 * newBoxDetector(G4LogicalVolume *worldLogical,
                               G4String detName,
                               G4Material * detMaterialPtr,
                               G4ThreeVector centerPlacement,
                               G4ThreeVector boxHalfDimensions,
                               G4RotationMatrix* pRotationMatrix = 0);

  
private:

  DetectorConstructionMessenger *messenger;

};

#endif

