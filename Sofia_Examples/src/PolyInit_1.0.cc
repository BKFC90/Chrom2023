// -*- C++ -*-

// This is a relatively simple code for generating a set of polymers
// in a LAMMPS environment

// To build it you will need the libconfig++ package for parsing the
// input file (copies in each of the Poly folders).

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

#include <libconfig.h++>

using namespace std;
using namespace libconfig;

struct Point3D {
  double x;
  double y;
  double z;
};

struct Bead {
  Point3D Location;
  int Type;
  int ID;
};

struct Polymer {
  Bead* BeadArray;
  void Init(int N)
  { BeadArray = new Bead[N]; };
  void Close()
  { delete[] BeadArray; };
};

int main(int argc, char** argv)
{

  Config Configuration;

  ofstream configstream;
  configstream.open("poly.init");
  
  int PolymerLength;
  
  int NPolymerBaseX;
  int NPolymerBaseY;
  int NPolymerBaseZ;

  double BeadSpacing;
  
  double Xmax;
  double Ymax;
  double Zmax;

  double Mass1;
  double Mass2;

  try
    {
      Configuration.readFile("ConfigPoly.dat");
    }
  catch(const FileIOException &fioex)
    {
      std::cerr << "I/O error while reading file." << std::endl;
      return(EXIT_FAILURE);
    }
  catch(const ParseException &pex)
    {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
		<< " - " << pex.getError() << std::endl;
      return(EXIT_FAILURE);
    }

  // get the parameters

  try
    {

      PolymerLength = Configuration.lookup("PolymerLength");

      NPolymerBaseX = Configuration.lookup("NPolymerBaseX");
      NPolymerBaseY = Configuration.lookup("NPolymerBaseY");
      NPolymerBaseZ = Configuration.lookup("NPolymerBaseZ");

      BeadSpacing = Configuration.lookup("BeadSpacing");
      Xmax = Configuration.lookup("Xmax");
      Ymax = Configuration.lookup("Ymax");
      Zmax = Configuration.lookup("Zmax");
      Mass1 = Configuration.lookup("Mass1");
      Mass2 = Configuration.lookup("Mass2");
    }
    catch(const SettingNotFoundException &MissingParameter)
      {
	cerr << "Parameter " << MissingParameter.getPath()
	     << " not found." << endl;
	return(EXIT_FAILURE);
      }
    catch(const SettingTypeException &BadParameter)
      {
	cerr << "Parameter " << BadParameter.getPath()
	     << " has wrong type." << endl;
	return(EXIT_FAILURE);
      }

  // header

  configstream << "LAMMPS chain data file" << endl;

  configstream << endl;

  int NPolymer = NPolymerBaseX * NPolymerBaseY * NPolymerBaseZ;


  int BeadCount = NPolymer * PolymerLength;
  int BondCount = NPolymer * (PolymerLength - 1);
  int BeadShift = PolymerLength / 2;

  // internal data

  double SimulationVolume = 8 * Xmax * Ymax * Zmax;

  configstream << "# polymer length:  " << PolymerLength << endl;
  configstream << "# " << NPolymer << " polymers" << endl;
  configstream << "# simulation volume: " << SimulationVolume << endl;

  configstream << endl;
  
  configstream << BeadCount << " atoms" << endl;
  configstream << BondCount << " bonds" << endl;
  configstream << 0 << " angles" << endl;
  configstream << 0 << " dihedrals" << endl;
  configstream << 0 << " impropers" << endl;

  configstream << endl;

  configstream << "2 atom types" << endl;
  configstream << "1 bond types" << endl;
  configstream << "0 angle types" << endl;
  configstream << "0 dihedral types" << endl;
  configstream << "0 improper types" << endl;

  configstream << endl;
  
  configstream << -Xmax << " " << +Xmax << " xlo xhi" << endl;
  configstream << -Ymax << " " << +Ymax << " ylo yhi" << endl;
  configstream << -Zmax << " " << +Zmax << " zlo zhi" << endl;

  configstream << endl;

  configstream << "Masses" << endl;

  configstream << endl;

  configstream << "1 " << Mass1 << endl;
  configstream << "2 " << Mass2 << endl;

  configstream << endl;

  // now set up array of polymers

  Polymer* PolymerArray = new Polymer[NPolymer];
  for (int PolymerIndex = 0; PolymerIndex < NPolymer; PolymerIndex++) {
    PolymerArray[PolymerIndex].Init(PolymerLength);
  }    

  double XCell = 2 * Xmax / NPolymerBaseX;
  double YCell = 2 * Ymax / NPolymerBaseY;
  double ZCell = 2 * Zmax / NPolymerBaseZ;

  int BeadIndex = 0;
  double PolymerIndexShiftX = NPolymerBaseX / (double (2));
  double PolymerIndexShiftY = NPolymerBaseY / (double (2));
  double PolymerIndexShiftZ = NPolymerBaseZ / (double (2));

  for (int iz = 0; iz < NPolymerBaseZ; iz++) 
    for (int iy = 0; iy < NPolymerBaseY; iy++)
      for (int ix = 0; ix < NPolymerBaseX; ix++) {

	int PolymerIndex = iz * NPolymerBaseY * NPolymerBaseX + iy * NPolymerBaseX + ix;
	
	// polymer center
	
	double PolymerX = (ix - PolymerIndexShiftX + 0.5) * XCell
	  + (((rand() % 100) / double(100.0)) - 0.5) * 5 * BeadSpacing;
	double PolymerY = (iy - PolymerIndexShiftY + 0.5) * YCell;
	double PolymerZ = (iz - PolymerIndexShiftZ + 0.5) * ZCell;

	// run the length of the polymer
	
	for (int Index = 0; Index < PolymerLength; Index++) {
	  BeadIndex++;
	  PolymerArray[PolymerIndex].BeadArray[Index].ID = BeadIndex;
	  PolymerArray[PolymerIndex].BeadArray[Index].Location.x
	    = PolymerX + (Index - BeadShift + 0.2) * BeadSpacing;
	  PolymerArray[PolymerIndex].BeadArray[Index].Location.y
	    = PolymerY;
	  PolymerArray[PolymerIndex].BeadArray[Index].Location.z
	    = PolymerZ;

	  //	  if ((Index % 2) == 0)
	  if (Index < (PolymerLength / 2))
	    PolymerArray[PolymerIndex].BeadArray[Index].Type = 1;
	  else
	    PolymerArray[PolymerIndex].BeadArray[Index].Type = 2;
	}
      }

  configstream << endl;

  // now write out the bead coordinates and assign type

  configstream << "Atoms" << endl << endl;

  for (int PolymerIndex = 0; PolymerIndex < NPolymer; PolymerIndex++)
    for (int Index = 0; Index < PolymerLength; Index++)

    configstream << PolymerArray[PolymerIndex].BeadArray[Index].ID
	 << '\t'
	 << PolymerIndex + 1
	 << '\t'
	 << PolymerArray[PolymerIndex].BeadArray[Index].Type
	 << '\t'
	 << PolymerArray[PolymerIndex].BeadArray[Index].Location.x
	 << '\t'
	 << PolymerArray[PolymerIndex].BeadArray[Index].Location.y
	 << '\t'
	 << PolymerArray[PolymerIndex].BeadArray[Index].Location.z
	 << endl;

  // now write out velocities

  configstream << endl << "Velocities" << endl << endl;

  for (int PolymerIndex = 0; PolymerIndex < NPolymer; PolymerIndex++)
    for (int Index = 0; Index < PolymerLength;  Index++)

    configstream << PolymerArray[PolymerIndex].BeadArray[Index].ID
	 << '\t'
	 << 0
	 << '\t'
	 << 0
	 << '\t'
	 << 0
	 << '\t'
	 << endl;

  // now write out the bonds
  
  configstream << endl << "Bonds" << endl << endl;

  int BondNumber = 0;
  
  for (int PolymerIndex = 0; PolymerIndex < NPolymer; PolymerIndex++)
    for (int Index = 0; Index < PolymerLength - 1; Index++) {

      BondNumber++;
  
      configstream << BondNumber
      	 << '\t'
	 << 1
	 << '\t'
	 << PolymerArray[PolymerIndex].BeadArray[Index].ID
	 << '\t'
	 << PolymerArray[PolymerIndex].BeadArray[Index].ID + 1
	 << endl;
    }

  delete [] PolymerArray;
  configstream.close();
  return 0;
}
