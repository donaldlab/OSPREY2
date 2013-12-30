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
//	ObjectiveFunction.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import cern.colt.matrix.DoubleMatrix1D;

/**
 *
 * @author mhall44
 */
//Objective function for CCDMinimizer or another general minimizer
public interface ObjectiveFunction {

    //Return number of degrees of freedom
    public int getNumDOFs();


    //Get constraints on the degrees of freedom
    public DoubleMatrix1D[] getConstraints();

    //Set all DOFs to values in x (e.g. in the molecule)
    public void setDOFs(DoubleMatrix1D x);

    //Set just one degree of freedom
    public void setDOF(int dof, double val);

    //Value and gradient at a given point (specified as values for all DOFs)
    public double getValue(DoubleMatrix1D x);
    //public DoubleMatrix1D getGradient(DoubleMatrix1D x);

    //Value at a given value for a given DOF,
    //and, for efficiency, possibly omitting energy terms that don't depend on that DOF
    //Other DOFs kept as they are currently set
    //(omitted terms must be consistent for a given objective function and DOF)
    public double getValForDOF(int dof, double val);


    public double getInitStepSize(int dof);//Get a good initial step size for a DOF
    //(e.g. for first initial value checking in CCD)
    
    
    public boolean isDOFAngle(int dof);//Is the given degree of freedom an angle?
}
