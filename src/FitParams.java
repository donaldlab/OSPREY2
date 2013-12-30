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
//	FitParams.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

//Fit parameter class, used by NWisePolyFitter
public class FitParams {
    
    static FitParams CLASH = new FitParams();//markers for clashes
    
    
    int numDims;
    
    int order;//order of polynomials
    
    //stuff for principal components (activated if PCOrder>order)
    //boolean isPC[];
    double PCFac;//eigenvalue factor used to choose components
    int PCOrder;//max order of principal components
    
    
    boolean includeConst;
            
    double explicitVDWCutoff;
    
    int numPCParams = -1;//number of parameters from only-PC monomial degrees

    //set to -1 at first to indicate not yet computed; when we have an EPolyPC object
    //we can make this correct
    
    
    
    public FitParams(){
        
    }
    
    
    public FitParams(int numDims, int order, double PCFac, int PCOrder, boolean includeConst, double explicitVDWCutoff) {
        this.numDims = numDims;
        this.order = order;
        this.PCFac = PCFac;
        this.PCOrder = PCOrder;
        this.includeConst = includeConst;
        this.explicitVDWCutoff = explicitVDWCutoff;
    }
    
    
    
    
    static FitParams quadratic(int nd, boolean includeConst) {
        return new FitParams(nd,2,0,2,includeConst,0);
    }
    
    
    
    
    void increaseOrder(){
        //go up to the next level!!
        //starts with quadratic
        //let's keep it simple for starters
        if(order==2)
            order = 4;
        /*else if(order==4)
            order = 6;*/
        else
            System.err.println("ERROR: Cannot increase order further");
    }
    
    
    
    int numSamplesNeeded(){
        return 10*numParams();
    }
    
    
    int numParams(){
        int ans = SeriesFitter.getNumParams(numDims, includeConst, order);
        if(PCOrder>order){
            if(numPCParams==-1){
                System.err.println("ERROR: Tried to call fitParams.numParams() on PC fit without"
                        + "precomputing numPCParams!");
                System.exit(1);
            }
            ans += numPCParams;   
        }
        return ans;
    }
    
    
    String getDescription(){
        
        if(this==CLASH)
            return "CLASH";
        
        String description = "FIT ORDER: "+order;
        if(PCOrder>order)//PC fit
            description = description+" PCORDER="+PCOrder+" FAC="+PCFac;
        if(explicitVDWCutoff>0)//has sve
            description = description+" EXPLICITVDWCUTOFF="+explicitVDWCutoff;
        if(includeConst)
            description = description+" INCLUDECONST";
        
        return description;
    }
    
    
}
