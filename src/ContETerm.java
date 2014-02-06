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
//	ContETerm.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import java.io.Serializable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;


//This is a base class for an energy term represented in a way that includes continuous variations
//for example, representing a pairwise or intra+shell energy as a polynomial in internal coordinates (EPIC)

public abstract class ContETerm implements Serializable {

    //The continuously varying degrees of freedom the terms acts on
    int DOFNums[];
    int numDOFs;
    //bounds on the degrees of freedom
    DoubleMatrix1D DOFmax, DOFmin;
    
    
    //Store discrete DOF values specific to the RC(s) represented in this ContETerm
    ArrayList<Integer> discreteDOFNums = new ArrayList<>();
    ArrayList<Double> discreteDOFVals = new ArrayList<>();
    
    
    
    DoubleMatrix1D center;//the center (typically, minimum point) for relative coordinates for this term
    double minE;//minimum energy for voxel represented by this (value of
    //pairwise or intra+shell center at DOFs=center)
    //will be included in evaluate if includeMinE is called, else center will evalute to 0
    
    String fitDescription = "N/A";//description of fit used to get this ContETerm, if applicable
    
    public ContETerm(int[] DOFNums, int numDOFs, DoubleMatrix1D DOFmax, DoubleMatrix1D DOFmin, DoubleMatrix1D center, double minE) {
        this.DOFNums = DOFNums;
        this.numDOFs = numDOFs;
        this.DOFmax = DOFmax;
        this.DOFmin = DOFmin;
        this.center = center;
        this.minE = minE;
    }


    public abstract double evaluate(DoubleMatrix1D x, boolean includeMinE, Molecule m);
    //evaluate at coordinates x, including minE if indicated
    //m is only used if there is SVE intended to use a shared molec; m should already be set
    //to match DOF values x, and may be null if not wanted or no SVE
    

    DoubleMatrix1D toRelCoords(DoubleMatrix1D x){
        //convert coordinates to be relative to center
        DoubleMatrix1D ans = x.copy();
        ans.assign(center,Functions.minus);
        return ans;
    }
    
    
    //reverse conversion
    DoubleMatrix1D toAbsCoords(DoubleMatrix1D z){
        //convert coordinates to be relative to center
        DoubleMatrix1D ans = z.copy();
        ans.assign(center,Functions.plus);
        return ans;
    }
    


    
    
    double bisectionBeta(DoubleMatrix1D x, double r){
        //computes scaling of (x-center) to get the value (from evaluate()) equal to r (w/o minE)
        //x is assumed to be in range
        
        DoubleMatrix1D relx = add( x, scale(center,-1) );
        double bottom = 0;//Won't scale below 0
        //figure out how high we can get by scaling all the way to the edge of the voxel, use as upper bound
        double top = Double.POSITIVE_INFINITY;
        //we minimize over DOFs because we want min_dofs (max scaling possible for dof)
        for(int dof=0; dof<numDOFs; dof++){
            if(relx.get(dof)<-0.001)//relx points in negative direction for dof
                top = Math.min( top, (DOFmin.get(dof)-center.get(dof))/relx.get(dof) );
            else if(relx.get(dof)>0.001)//positive direction
                top = Math.min( top, (DOFmax.get(dof)-center.get(dof))/relx.get(dof) );
        }        
        
        if(Double.isInfinite(top)){//x is basically at center
            if(r==0)
                return 0;
            else
                return Double.NaN;//can't really scale up 0 to anything
        }
            
        
        if( evaluate(add(center,scale(relx,top)),false,null) < r )//at top, r is still not achieved
            return Double.NaN;
        
        //now bisect to find the right scaling
        double mid = 0;
        
        double bisecPrecision = 1e-10;
        
        while ( top-bottom > bisecPrecision ) {
            mid = (top+bottom)/2;

            double val = evaluate(add(center,scale(relx,mid)),false,null);

            if( val>=r )
                top = mid;
            else
                bottom = mid;
        }
        
        return mid;
    }
    

    //Vector operations for convenience
    DoubleMatrix1D oneComponentVector(int n, double a){
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(numDOFs);
        ans.set(n, a);
        return ans;
    }

    static DoubleMatrix1D add(DoubleMatrix1D... v){
        DoubleMatrix1D ans = v[0].copy();
        for(int n=1; n<v.length; n++)
            ans.assign(v[n], Functions.plus);
        return ans;
    }

    static DoubleMatrix1D scale(DoubleMatrix1D v, double d){
        DoubleMatrix1D ans = v.copy();
        ans.assign(Functions.mult(d));
        return ans;
    }

    
    
    public boolean valuesInRange(DoubleMatrix1D x){
        //Check if x is in range for this voxel
        for(int dof=0; dof<numDOFs; dof++){
            if(x.get(dof)>DOFmax.get(dof) || x.get(dof)<DOFmin.get(dof))
                return false;
        }
        return true;
    }
    


}
