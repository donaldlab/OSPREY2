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
//	EPICSettings.java
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
import java.io.Serializable;


//This class contains the settings for EPIC: whether to use it, thresholds, etc.
//It will be referenced during all CETMatrix and A* runs.  

public class EPICSettings implements Serializable {
    
    boolean useEPIC;
    boolean checkEPIC;//do analysis of EPIC vs. regular energy for A*-enumerated conformations in GMEC runs
    double EPICThresh1;
    double EPICThresh2;
    
    boolean enforceEPICThresh;//if true, then make sure validity conditions involving EPIC thresholds are followed
    //this means keeping thresh2 at least Ival+Ew in DEE calculations, and complaining and quitting
    //if thresh1 or K* validity conditions aren't meant
    
    double EPICGoalResid;//(default = 1e-4)

        
    
    boolean useSVE = true;
    boolean usePC = true;
    
    
    
    
    //We'll need to precompute the lowest bound without polynomial fits for A* GMEC calculations
    boolean gettingLowestBound = false;//temporarily set to true for the precomputation 
    double lowestBound;//set during precomputation
    
    
    public EPICSettings(){
        //by default, no EPIC
        //this is cool for operations like K* mutation list or perturbation selection
        useEPIC = false;
        checkEPIC = false;
    }

    
    public EPICSettings(ParamSet params){
        //initialize from input parameter set
        useEPIC = (new Boolean((String)params.getValue("USEEPIC", "false"))).booleanValue();
        checkEPIC = (new Boolean((String)params.getValue("CHECKEPIC", "false"))).booleanValue();
        EPICThresh1 = (new Double((String)params.getValue("EPICTHRESH1","10"))).doubleValue();
        EPICThresh2 = (new Double((String)params.getValue("EPICTHRESH2","25"))).doubleValue();
        enforceEPICThresh = (new Boolean((String)params.getValue("ENFORCEEPICTHRESH", "true"))).booleanValue();
        EPICGoalResid = (new Double((String)params.getValue("EPICGOALRESID","0.0001"))).doubleValue();
        
        useSVE = (new Boolean((String)params.getValue("EPICUSESVE", "true"))).booleanValue();
        usePC = (new Boolean((String)params.getValue("EPICUSEPC","true"))).booleanValue();

        
        if(EPICThresh2<EPICThresh1){
            System.err.println("ERROR: EPICThresh2 must be at least EPICThresh1!  EPICThresh2="+EPICThresh2+" EPICThresh1="+EPICThresh1);
            System.exit(1);
        }
    }
    
    
    
    void checkThresh1Validity(CETObjFunction lof, DoubleMatrix1D bestVals){
        //Given a CETObjFunction and its optimal values bestVal, if need to check threshold 1 (EPICThreshMode>0),
        //do so
        //basically no polynomial may exceed EPICThresh1
        //may consider more elegant error handling later
        
        
        if(enforceEPICThresh){
            lof.setDOFs(bestVals);
            double maxTerm = lof.ef.maxTerm();

            if(maxTerm>EPICThresh1){
                System.err.println("ERROR: Insufficient EPICThresh1="+EPICThresh1+".  Encountered greater polynomial value: "+maxTerm);
                System.exit(1);
            }
        }
        
    }
    
    /* The validity condition for EPICThresh2 is that 
      we can only trust conformational enumeration up to 
      an energy of lowestBound+EPICThresh2
      For non-K* A*, we can ensure this just by setting EPICThresh2 to be at least Ival+Ew
      for K*, we complain if we find ourselves enumerating a conformation with 
      energy > lowestBound + EPICThresh2 */
    
    
    
    
}
