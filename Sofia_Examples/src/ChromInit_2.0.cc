// -*- C++ -*-

// This is a setup code for a model of two chromsomes, which are
// initiallly configured into small loops before expansion via soft
// potentials.

// The model includes 1 centromere tethered to an SPB and 1 telomeric
// bead at each end of the chromosome.

// To build this code, you will need the libconfig++ and armadillo libraries.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

#include <armadillo>

#include <libconfig.h++>

using namespace std;
using namespace arma;
using namespace libconfig;

enum BEADTYPE {GENONE = 1, GENTWO = 2, TEL = 3, CEN = 4, SPB = 5, FAKE = 6};

struct Bead {
  vec Location;
  int Type;
  int ID;
  Bead()  { Location = vec(3); }
};

//  what we need:
// array of loops
// each loop is one chromosome, with its own length and centromere location
// rDNA has different parameters
//

class Chromosome {
protected:
  int ChrLen;
  int CenLoc;
  
public:
  Bead* BeadArray;
  Chromosome() {}
  void Initialize(int ChrID, int ChrLen_, int CenLoc_)
  {
    ChrLen = ChrLen_;
    CenLoc = CenLoc_;

    int rDNALoc = 564;
    rDNALoc--;
    int rDNALen = 938;
    
    BeadArray = new Bead[ChrLen];
    for (int j = 0; j < ChrLen; j++) {
      BeadArray[j].Type = 0;

      // 1 telomeric bead at each end (independent of scale)
      if (j == 0)
	BeadArray[j].Type = TEL;
      if (j  == ChrLen - 1)
	BeadArray[j].Type = TEL;
      if (j == CenLoc - 1)
	BeadArray[j].Type = CEN;
    }
      
    
  }
  ~Chromosome()
  { delete[] BeadArray; }
	
protected:
};

int main(int argc, char** argv)
{

  Config Configuration;

  ofstream headstream;
  headstream.open("Header");
  ofstream bodystream;
  bodystream.open("Body");
  int NChr = 2;
  double NuclearRadius = 30;
  double box_half_edge = 35;

  Chromosome* ChromosomeArray = new Chromosome[NChr];

  // customize bead structure here

  // chromosome lengths
  
  int ChrLen[2] = {700, 800};

  // centromere locations
  
  int CenLoc[2] = {300, 250};

  double dLat = 2 * M_PI / NChr;

  // initial bead spacing
  
  double BeadSpacing = 0.10;

  int Nphi;
  int Ntheta;
  
  double Xmax;
  double Ymax;
  double Zmax;

  double Mass1;
  double Mass2;
  double Mass3;

  try
    {
      Configuration.readFile("ConfigChrom.dat");
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
      BeadSpacing = Configuration.lookup("BeadSpacing");
      NuclearRadius = Configuration.lookup("NuclearRadius");
      box_half_edge = Configuration.lookup("box_half_edge");
      Nphi = Configuration.lookup("Nphi");
      Ntheta = Configuration.lookup("Ntheta");
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

  for (int ChrIndex = 0; ChrIndex < NChr; ChrIndex++) {
    ChromosomeArray[ChrIndex].Initialize(ChrIndex, ChrLen[ChrIndex], CenLoc[ChrIndex]);
  }

 // header

  bodystream << "Atoms" << endl << endl;

  // internal data

    int BeadCount = 0;
    
  // now set up array of chromosomes

  // loop over chromosome count

  for (int ChrIndex = 0; ChrIndex < NChr; ChrIndex++) {

    // set latitude of plane for chromosome loop
					   
    double Latitude = (2 * M_PI / NChr) * ChrIndex;

    double CL = cos(Latitude);
    double SL = sin(Latitude);

    mat LatRotate(3, 3);
    // LatRotate << CL << -SL << 0 << endr
    // 	      << SL << CL << 0 << endr
    // 	      << 0 << 0 << 1;
    LatRotate << 1 << 0 << 0 << endr
	      << 0 << CL << SL << endr
	      << 0 << -SL << CL;
    double ChrLength = ChrLen[ChrIndex] * BeadSpacing;
    double LoopRadius = ChrLength / (2 * M_PI);

    double BeadAngleSpacing = 2 * M_PI / ChrLen[ChrIndex];

    // run through chromosome, assigning locations

    for (int BeadIndex = 0; BeadIndex < ChrLen[ChrIndex]; BeadIndex++) {

      BeadCount++;
	
	ChromosomeArray[ChrIndex].BeadArray[BeadIndex].ID = BeadCount;
	
	vec BeadLoc(3);
	double Offset = 5 * BeadSpacing;
	double BeadAngle = BeadIndex * BeadAngleSpacing;
	double Z = Offset + LoopRadius * (1 + cos(BeadAngle));
	double X =  -LoopRadius * sin(BeadAngle);
	vec PlaneVec{X, 0, Z};
	vec RotatedVec = LatRotate * PlaneVec;
	ChromosomeArray[ChrIndex].BeadArray[BeadIndex].Location
	  = RotatedVec;

      // types already assigned

      }
    }    

  // now write out the bead coordinates

  double Rmax = 0;
  for (int ChrIndex = 0; ChrIndex < NChr; ChrIndex++)
    for (int BeadIndex = 0; BeadIndex < ChrLen[ChrIndex]; BeadIndex++) {

      int LocalType = ChromosomeArray[ChrIndex].BeadArray[BeadIndex].Type;
      if (LocalType == 0) {
	if (ChrIndex == 0) LocalType = GENONE;
	if (ChrIndex == 1) LocalType = GENTWO;
      }
      bodystream << ChromosomeArray[ChrIndex].BeadArray[BeadIndex].ID
	 << '\t'
	 << ChrIndex + 1
	 << '\t'
	 << LocalType
	 << '\t'
	 << ChromosomeArray[ChrIndex].BeadArray[BeadIndex].Location(0)
	 << '\t'
	 << ChromosomeArray[ChrIndex].BeadArray[BeadIndex].Location(1)
	 << '\t'
	 << ChromosomeArray[ChrIndex].BeadArray[BeadIndex].Location(2)
	 << endl;

    double Rtest
      = sqrt(pow(ChromosomeArray[ChrIndex].BeadArray[BeadIndex].Location(0), 2)
		 + pow(ChromosomeArray[ChrIndex].BeadArray[BeadIndex].Location(1), 2)
		 + pow(ChromosomeArray[ChrIndex].BeadArray[BeadIndex].Location(2), 2));
    if (Rtest > Rmax)
      Rmax = Rtest;
    }
  cerr << "Rmax: " << Rmax << endl;
  // write out the SPB separately

  {
    BeadCount++;
    bodystream << BeadCount
	 << '\t'
	 << NChr + 1
	 << '\t'
	 << SPB
	 << '\t'
	 << -NuclearRadius
	 << '\t'
	 << 0
	 << '\t'
	 << 0
	 << endl;
  }

  int SPBID = BeadCount;

  // create fake atoms for spherical boundary

  double dphi = 2 * M_PI / Nphi;
  double dtheta = M_PI /  (2 * Ntheta);

  for (int jphi = 0; jphi < Nphi; jphi++)
    {
      double phi = dphi * jphi;
      double Cp = cos(phi);
      double Sp = sin(phi);
      for (int jtheta = -Ntheta; jtheta <= Ntheta; jtheta++)
	{
	  double theta = dtheta * jtheta;
	  double Ct = cos (theta);
	  double St = sin(theta);
	  double X = NuclearRadius * Ct * Cp;
	  double Y = NuclearRadius * Ct * Sp;
	  double Z = NuclearRadius * St;
	  BeadCount++;
	  bodystream << BeadCount
	       << '\t'
	       << NChr + 2
	       << '\t'
	       << FAKE
	       << '\t'
	       << X
	       << '\t'
	       << Y
	       << '\t'
	       << Z
	       << endl;
	}
    }
  bodystream << endl;

  // now write out simple bond relationships

  bodystream << "Bonds" << endl << endl;
  
  int BondCount = 0;
  int BondType;
  for (int ChrIndex = 0; ChrIndex < NChr; ChrIndex++)
    for (int BeadIndex = 0; BeadIndex < ChrLen[ChrIndex] - 1; BeadIndex++) {
      BondCount++;
      BondType = 1;
      bodystream << BondCount
	   << '\t'
	   << BondType
	   << '\t'
	   << ChromosomeArray[ChrIndex].BeadArray[BeadIndex].ID
	   << '\t'
	   << ChromosomeArray[ChrIndex].BeadArray[BeadIndex].ID + 1
	   << endl;
    }

  // add the centromere-SPB bonds

  int SPB_Bond = 2;

  for (int ChrIndex = 0; ChrIndex < NChr; ChrIndex++) {
    {
      BondCount++;
      bodystream
	<< BondCount
	<< '\t'
	<< SPB_Bond
	<< '\t'
	<< SPBID
	<< '\t'
	<< ChromosomeArray[ChrIndex].BeadArray[(CenLoc[ChrIndex] - 1)].ID
	<< endl;
    }
  }

  bodystream << endl;

  // now make the header

  headstream << "LAMMPS datafile" << endl << endl;
  headstream << "#" << endl;
  headstream << "# nuclear radius: " << NuclearRadius << endl;
  headstream << "# box half edge: " << box_half_edge << endl;
  headstream << "# bead spacing: " << BeadSpacing << endl;
  headstream << "#" << endl;
  headstream << endl;
  
  headstream << BeadCount << '\t' << "atoms" << endl;
  headstream << BondCount << '\t' << "bonds" << endl;
  headstream << 0 << '\t' << "angles" << endl;
  headstream << 0 << '\t' << "dihedrals" << endl << endl;

  headstream << 6 << '\t' << "atom types" << endl;
  headstream << 2 << '\t' << "bond types" << endl;
  headstream << 0 << '\t' << "angle types" << endl;
  headstream << 0 << '\t' << "dihedral types" << endl << endl;

  headstream << -box_half_edge << '\t' << box_half_edge << '\t'
	     << "xlo" << '\t' << "xhi" << endl;
  headstream << -box_half_edge << '\t' << box_half_edge << '\t'
	     << "ylo" << '\t' << "yhi" << endl;
  headstream << -box_half_edge << '\t' << box_half_edge << '\t'
	     << "zlo" << '\t' << "zhi" << endl;
  headstream << endl;
  
  delete [] ChromosomeArray;
  bodystream.close();
  headstream.close();
  return 0;
}
