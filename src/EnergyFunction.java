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
//	EnergyFunction.java
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


//An EnergyFunction is a scalar function of the molecular coordinates
//Examples: pairwise energy between two residues, total energy of molecule,
//dihedral energy, energy with some arithmetic done on it for convex lower-bound computation, etc.
public abstract class EnergyFunction implements Serializable {

    //Get energy
    public abstract double getEnergy();

    //During minimization we want to calculate only the part of the energy
    //involving certain residues without creating a separate EnergyFunction object
    //This function sets up this computation; residues[a] is the a'th set of residues
    //to calculate energies involving.  These are molecule residue numbers.  
    public abstract void setupPartialComputation(int residues[][]);
    //And this function calculates energies involving the a'th set of residues
    public abstract double getEnergy(int a);

    //NOTE: some subclasses may not support partial computation and will return the entire energy regardless


    //Add a term to the energy function by creating a MultiTermEnergyFunction
    //If this is already a MultiTermEnergyFunction then the term will be added directly
    //but this is in the overriding function MultiTermEnergyFunction.addTerm
    public EnergyFunction addTerm(EnergyFunction ef){
        EnergyFunction terms[] = { this, ef };
        double coeffs[] = {1,1};
        return new MultiTermEnergyFunction(terms, coeffs);
    }

    //This function is implemented by energy functions that have an Amber96ext, 
    //which override this definition.  Without such an implementation there needs to be an error
    public Amber96ext getAmber96ext(){
        throw new RuntimeException("ERROR: getAmber96ext() operation unsupported for "+
                getClass().getName());
    }


}
