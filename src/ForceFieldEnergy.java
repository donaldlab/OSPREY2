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
//	ForceFieldEnergy.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

//This is an energy function based on a force field
//It wraps an a96ff for use as an energy function (suitable for adding other terms,
//quasi-Newton minimization, etc.)
public class ForceFieldEnergy extends EnergyFunction {

    Molecule m;

    Amber96ext a96ff;

    public ForceFieldEnergy(Molecule molec, Amber96ext ff){
        m=molec;
        a96ff=ff;
    }

    public double getEnergy(){
        if( ! (m.validBB&&m.validSC) )//Invalid conformation
            return Double.POSITIVE_INFINITY;

        return a96ff.calculateTotalEnergy(m.actualCoordinates, -1)[0];
    }


    public void setupPartialComputation(int res[][]){

        //public void setupPartialArrays(int numRows, int maxNumColumns, int atomList[][],int numColumns[]){

        int numRows = res.length;

        int atomList[][] = new int[numRows][];
        int numColumns[] = new int[numRows];
        int maxNumColumns = 0;

        for(int row=0; row<numRows; row++){

            for(int molResNum : res[row])
                numColumns[row] += m.residue[molResNum].numberOfAtoms;

            maxNumColumns = Math.max(maxNumColumns, numColumns[row]);

            atomList[row] = new int[numColumns[row]];

            int atCount=0;

            for( int molResNum : res[row] ){
                Residue curRes = m.residue[molResNum];

                for( Atom at : curRes.atom ){
                    atomList[row][atCount] = at.moleculeAtomNumber;
                    atCount++;
                }
            }
        }


        a96ff.setupPartialArrays(numRows,maxNumColumns,atomList,numColumns);
    }


    //Partial energy
    public double getEnergy(int a){
        if( ! (m.validBB&&m.validSC) )//Invalid conformation
            return Double.POSITIVE_INFINITY;
        
        return a96ff.calculateTotalEnergy(m.actualCoordinates, a)[0];
    }


    
  @Override
    public Amber96ext getAmber96ext(){
        return a96ff;
    }

}
