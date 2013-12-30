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
//	ReducedEnergyMatrix.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.io.Serializable;
import java.util.Arrays;

//This class represents the energy bounds for unpruned rotamers in a compact format
//for use in A*
public class ReducedEnergyMatrix implements Serializable {

    //Residues, AA types, and normal rotamer indices for each of the reduced rotamer indices
    int indicesEMatrixPos[];//Numbered among flexible residues
    int indicesEMatrixAA[];//Numbered among all possible AA types
    int indicesEMatrixRot[];//Numbered among all RCs for a given residue and AA type

    int shortRedIndices[][][];//Reverse lookup for reduced indices (per residue, as a function of pos, AA, rot)

    int numRes;//number of residues

    private double RedEmat[][];
    // the min energy matrix: the last column contains the intra-energy for each rotamer; the last row
    // contains the shell-residue energy for each rotamer

    
    int numRotTotal;//the total number of possible rotamers for the given mutation

    int resOffsets[];//Rotamer/RC index in pairwiseEnergyMatrix at which each new residue starts

    

    public ReducedEnergyMatrix(int numRotamers, double[][] a, int indPos[], int indAA[], int indRot[]){
        numRotTotal = numRotamers;
        RedEmat = a;
        indicesEMatrixPos = indPos;
        indicesEMatrixAA = indAA;
        indicesEMatrixRot = indRot;


        numRes = indicesEMatrixPos[indicesEMatrixPos.length-1]+1;

        //Calculate resOffsets
        resOffsets = new int[numRes];
        int place=0;

        for(int curRes=0; curRes<numRes-1; curRes++){
            resOffsets[curRes] = place;
            while( indicesEMatrixPos[place] == curRes ){
                place++;
            }
        }
        resOffsets[numRes-1] = place;


        //Compute index reverse-lookup array
        //Assuming indicesEMatrixPos, etc. are in ascending order
        shortRedIndices = new int[numRes][][];
        for(int res=0; res<numRes; res++){
            int lastRC;//Last (long reduced index) RC for this residue
            if(res==numRes-1)
                lastRC = indicesEMatrixPos.length-1;
            else
                lastRC = resOffsets[res+1]-1;

            int lastAA = indicesEMatrixAA[lastRC];//biggest AA index we'll need

            shortRedIndices[res] = new int[lastAA+1][];

            //Count how many rotamers are at each AA type
            int p=resOffsets[res];
            for(int curAA=0; curAA<=lastAA; curAA++){
                lastRC = -1;
                while(indicesEMatrixAA[p]==curAA && indicesEMatrixPos[p]==res){
                    lastRC = indicesEMatrixRot[p];//Largest RC number for this AA type
                    p++;
                    if(p>=indicesEMatrixPos.length)
                        break;
                }
                shortRedIndices[res][curAA] = new int[lastRC+1];
                Arrays.fill( shortRedIndices[res][curAA], -1 );//-1 for rotamers that aren't available for our A* tree
            }

            //Fill in the indices
            
            p=resOffsets[res];
            for(int curAA=0; curAA<=lastAA; curAA++){
                while( indicesEMatrixAA[p] == curAA && indicesEMatrixPos[p]==res){
                    shortRedIndices[res][curAA][indicesEMatrixRot[p]] = p-resOffsets[res];
                    p++;
                    if(p>=indicesEMatrixPos.length)
                        break;
                }
            }
        }
    }

    public double getPairwiseE(int index1, int index2){
        return RedEmat[index1][index2];
    }

    public double getIntraE(int index){//the intra-energy is in the last column
        return RedEmat[index][numRotTotal];
    }

    public double getShellRotE(int index){
        double shlRotE = 0.0f;
        //Skip the first energy which is the intra-rot energy still
        for(int i=numRotTotal; i<RedEmat.length;i++){
                shlRotE += RedEmat[i][index];
        }
        return shlRotE;
    }

    public double getIntraAndShellE(int index){
        return getIntraE(index) + getShellRotE(index);
    }


            


    //Getting energies based on (residue, short reduced index) pairs
    public double getPairwiseE(int res1, int res2, int RC1, int RC2){
        int offset1 = resOffsets[res1];
        int offset2 = resOffsets[res2];
        return RedEmat[offset1+RC1][offset2+RC2];
    }

    public double getIntraE(int res, int RC){//the intra-energy is in the last column
        int offset = resOffsets[res];
        return RedEmat[offset+RC][numRotTotal];
    }

    public double getShellRotE(int res, int RC){
        int offset = resOffsets[res];
        return getShellRotE(offset+RC);
    }

    public double getIntraAndShellE(int res, int RC){
        return getIntraE(res,RC) + getShellRotE(res,RC);
    }


}