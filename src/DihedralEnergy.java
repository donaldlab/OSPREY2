/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University

	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.

	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.

	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129
			USA
			e-mail:   www.cs.duke.edu/brd/

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
 */

///////////////////////////////////////////////////////////////////////////////////////////////
//	DihedralEnergy.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.util.Arrays;


//This class evaluates the dihedral energy
//defined as in SimpleMinimizer (e.g. defined to be 0 at ideal rotamers)
//It is created and used by a ContSCObjFunction
//If terms based on rotamer probabilities (not 0 at ideal rotamers) are added later,
//the terms in here now can be added to them
public class DihedralEnergy extends EnergyFunction {

    Molecule m;


    int strDihedralAtNums[][][] = null;
		// array of dihedrals, this array is n by 4
		//  and has the moleculeAtomNumbers for each
		//  of the 4 atoms of the dihedral
    int numberOfStrands = 0;


    
    int strDihedPN[][] = new int[0][0];
    double strDihedWeight[][] = new double[0][0];
    int strDihedLocalNum[][] = new int[0][0];
    int strDihedNumTerms[] = new int[0];

    boolean includeDihedral[][];//which dihedrals to include

    double idealDihedVal[][] = null;//Ideal values for the dihedrals in their current rotamers
    //These will be initialized or updated when ContSCObjFunction.updateIdealDihedrals() is called


    //Based on SimpleMinimizer.setupDihedralTerms
    public DihedralEnergy(Molecule molec, int numStrands, int dihAtNums[][][],
            int numStrDihedrals[], Amber96ext a96ff){
    //The arguments are all things a ContSCObjFunction would know

        m = molec;
        strDihedralAtNums = dihAtNums;
        numberOfStrands = numStrands;


            // At this point we assume sysDihedralAtNums, ligDihedralAtNums,
            //  numSysDihedrals, and numLigDihedrals are all correct.

            int tmpStrSize = 1000;
            //int tmpLSize = 1000;
            strDihedPN = new int[numberOfStrands][tmpStrSize];
            //sDihedPN = new int[tmpSSize];
            //lDihedPN = new int[tmpLSize];
            strDihedWeight = new double[numberOfStrands][tmpStrSize];
            //sDihedWeight = new double[tmpSSize];
            //lDihedWeight = new double[tmpLSize];
            strDihedLocalNum = new int[numberOfStrands][tmpStrSize];
            //sDihedLocalNum = new int[tmpSSize];
            //lDihedLocalNum = new int[tmpLSize];
            strDihedNumTerms = new int[numberOfStrands];
            for(int i=0; i<numberOfStrands;i++)
                    strDihedNumTerms[i] = 0;
            //sDihedNumTerms = 0;
            //lDihedNumTerms = 0;
            int atom2 = 0, atom3 = 0;
            double fConst[] = new double[5];
            double eqAngle[] = new double[5];
            int terms[] = new int[5];
            int multiplic[] = new int[1];

            int newDihedPN[] = null;
            double newDihedWeight[] = null;
            int newDihedLocalNum[] = null;
            for(int str=0;str<numberOfStrands;str++){
            // SYSTEM DIHEDRALS
                    for(int i=0;i<numStrDihedrals[str];i++) {
                    // Get center two atoms
                            atom2 = strDihedralAtNums[str][i][1];
                            atom3 = strDihedralAtNums[str][i][2];
                    for(int q=0;q<m.connected[atom2][0];q++) {
                            if (m.connected[atom2][q+1] != atom3) {
                                    for(int w=0;w<m.connected[atom3][0];w++) {
                                            if (m.connected[atom3][w+1] != atom2) {
                                                    // At this point 'q',atom2,atom3,'w' is a dihedral
                                                    if (!(a96ff.getTorsionParameters(m.atom[m.connected[atom2][q+1]].type,m.atom[atom2].type,
                                                                                                                     m.atom[atom3].type,m.atom[m.connected[atom3][w+1]].type,
                                                                                                                     fConst,eqAngle,terms,multiplic))) {
                                                            System.out.println("WARNING: Could not find torsion parameters for " + m.connected[atom2][q+1] + "," + atom2 + "," + atom3 +"," + m.connected[atom3][w+1]);
                                                            System.exit(1);
                                                    }
                                                    for(int r=0;r<=multiplic[0];r++) {
                                                                    strDihedNumTerms[str]++;
                                                                    if (strDihedNumTerms[str] > tmpStrSize) {
                                                                    // increase the size of the arrays
                                                                            newDihedPN = new int[tmpStrSize + 500];
                                                                            newDihedWeight = new double[tmpStrSize + 500];
                                                                            newDihedLocalNum = new int[tmpStrSize + 500];
                                                                            System.arraycopy(strDihedPN[str],0,newDihedPN,0,tmpStrSize);
                                                                            System.arraycopy(strDihedWeight[str],0,newDihedWeight,0,tmpStrSize);
                                                                            System.arraycopy(strDihedLocalNum[str],0,newDihedLocalNum,0,tmpStrSize);
                                                                            strDihedPN[str] = newDihedPN;
                                                                            strDihedWeight[str] = newDihedWeight;
                                                                            strDihedLocalNum[str] = newDihedLocalNum;
                                                                            tmpStrSize += 500;
                                                            }
                                                                    strDihedPN[str][strDihedNumTerms[str]] = terms[r];
                                                                    strDihedWeight[str][strDihedNumTerms[str]] = fConst[r];
                                                                    strDihedLocalNum[str][strDihedNumTerms[str]] = i;
                                                    }
                                            }
                                    }
                            }
                    }
            }

            // Shrink the size of the array
                    newDihedPN = new int[strDihedNumTerms[str]];
                    newDihedWeight = new double[strDihedNumTerms[str]];
                    newDihedLocalNum = new int[strDihedNumTerms[str]];
                    System.arraycopy(strDihedPN[str],0,newDihedPN,0,strDihedNumTerms[str]);
                    System.arraycopy(strDihedWeight[str],0,newDihedWeight,0,strDihedNumTerms[str]);
                    System.arraycopy(strDihedLocalNum[str],0,newDihedLocalNum,0,strDihedNumTerms[str]);
                    strDihedPN[str] = newDihedPN;
                    strDihedWeight[str] = newDihedWeight;
                    strDihedLocalNum[str] = newDihedLocalNum;
            }

            includeDihedral = new boolean[numberOfStrands][];
            for(int str=0; str<numberOfStrands; str++){
                includeDihedral[str] = new boolean[strDihedralAtNums[str].length];
                Arrays.fill(includeDihedral[str],true);
            }
    }
    


    

    public double getEnergy(){

        double dihedE = 0.0;
        double d2r = 3.14159265 / 180.0;

        // System dihedral terms first
        for(int str=0; str<numberOfStrands;str++){
                for(int i=0;i<strDihedNumTerms[str];i++){

                        int j = strDihedLocalNum[str][i];//Indicates which dihedral this is

                        if(includeDihedral[str][j]){

                            double relDihedVal = m.getTorsion(strDihedralAtNums[str][j][0],strDihedralAtNums[str][j][1],
                                strDihedralAtNums[str][j][2],strDihedralAtNums[str][j][3]) - idealDihedVal[str][j];//Dihedral value minus ideal value for rotamer


                            dihedE += strDihedWeight[str][i] * (1 - Math.cos(strDihedPN[str][i]*relDihedVal*d2r));
                        }
                }
        }

        return(dihedE);

    }


    public void restrictToRes(int useStr, int useStrResNum){
        //Restrict the Dihedral Energy to dihedral in a particular residue
        for(int str=0; str<numberOfStrands; str++){
            for(int j=0; j<strDihedralAtNums[str].length; j++){

                if(str!=useStr)
                    includeDihedral[str][j] = false;
                else{
                    if( m.atom[strDihedralAtNums[str][j][0]].strandResidueNumber == useStrResNum )
                        includeDihedral[str][j] = true;
                    else
                        includeDihedral[str][j] = false;
                }
            }
        }
    }


    //Dihedral energies are fast so we are not going to provide partial computations
    public void setupPartialComputation(int res[][]){

    }

    public double getEnergy(int a){
        return getEnergy();
    }

    

}
