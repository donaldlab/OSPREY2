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
//	CETObjFunction.java
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
import cern.colt.matrix.DoubleFactory1D;
import java.io.Serializable;


//Objective function built from continuous energy terms
public class CETObjFunction implements ObjectiveFunction, Serializable {

    CETEnergyFunction ef;
    //Values of the DOFs will be stored in ef.DOFVals rather than in this objective function
    int numDOFs;
    DoubleMatrix1D constraints[];//constraints[0] holds the minima for all the DOFs; constraints[1], the maxima

    DegreeOfFreedom[] globalDOFList;//The global DOF list that the ef.DOFNums refer to

    //we can do SVE more efficiently if the SVE terms all share a molecule
    ContSCObjFunction sveOF;//sveOF will set DOFs in this molecule
    int sveOFMap[];//map from DOFs here to DOFs of sveOF
    int sveOFMapRev[];//reverse mapping
    
    public CETObjFunction(int[] DOFNums, double[][] constr, DegreeOfFreedom[] allDOFs,
            CETEnergyFunction efun, ContSCObjFunction sveObjFcn){
        //sveObjFcn can be null if we aren't doing SVE with molecule sharing

        numDOFs = DOFNums.length;
        ef = efun;
        constraints = new DoubleMatrix1D[] { DoubleFactory1D.dense.make(constr[0]), DoubleFactory1D.dense.make(constr[1]) };
        globalDOFList = allDOFs;
        
        sveOF = sveObjFcn;
        if(sveOF!=null){
            ef.m = sveOF.m;
            sveOFMap = new int[sveOF.numDOFs];
            sveOFMapRev = new int[numDOFs];
            for(int dof=0; dof<sveOF.numDOFs; dof++){                
                sveOFMap[dof] = ef.revDOFNums.get(sveOF.DOFNums[dof]);
                sveOFMapRev[sveOFMap[dof]] = dof;
            }
        }
    }

    public int getNumDOFs(){
        return numDOFs;
    }


    public DoubleMatrix1D[] getConstraints(){
        //This is called at the beginning of minimization
        //at this point, if we're doing strand translations and rotations
        //we need to get those back to the "ideal" state for the current rotamer
        //which is stored in the Atom.coord arrays--also with strand trans/rot info
        
        if(sveOF==null){
            for(ContETerm cet : ef.terms){
                if(cet instanceof EPoly){//not null and may have sve
                    EPoly s = (EPoly)cet;
                    if(s.sve!=null){
                        s.sveOF.getConstraints();//zeroes the trans/rot info
                        //adaptation of m.updateCoordinates since sve may have some null atoms to save space
                        for(Atom a : s.sve.m.atom){
                            if(a!=null)
                                s.sve.m.updateCoordinates(a);
                        }
                    }
                }
            }
        }
        else//going to use single sveOF instead of those in the individual terms
            sveOF.getConstraints();
        
        
        
        
        //To use if we're comparing SVE.getEnergy with and without argument (called from EPoly.evaluate)
        //DEBUG!!!
        /*for(ContETerm cet : ef.terms){
                if(cet instanceof EPoly){//not null and may have sve
                    EPoly s = (EPoly)cet;
                    if(s.sve!=null){
                        s.sveOF.getConstraints();//zeroes the trans/rot info
                        //adaptation of m.updateCoordinates since sve may have some null atoms to save space
                        for(Atom a : s.sve.m.atom){
                            if(a!=null)
                                s.sve.m.updateCoordinates(a);
                        }
                    }
                }
            }*/
        
        
               
        //the actual constraints are stored in the constraints field though
        return constraints;
    }
    

    //Set all DOFs to values in x (e.g. in the molecule)
    public void setDOFs(DoubleMatrix1D x){
        
        if(sveOF!=null){
            DoubleMatrix1D z = DoubleFactory1D.dense.make(sveOFMap.length);
            for(int dof=0; dof<sveOFMap.length; dof++)
                z.set(dof, x.get(sveOFMap[dof]));
            
            sveOF.setDOFs(z);
        }
        
        
        for(int a=0; a<numDOFs; a++)
            ef.DOFVals.set(a, x.get(a));
    }

    //Value and gradient at a given point (specified as values for all DOFs)
    public double getValue(DoubleMatrix1D x){
     
        setDOFs(x);
        double val = ef.getEnergy();

        if(Double.isNaN(val))
            System.err.println("NaN value returned by CETObjFunction!!");
        
        return val;
    }

    //Set just one degree of freedom
    public void setDOF(int dof, double val){
        
        if(sveOF!=null){
            sveOF.setDOF(sveOFMapRev[dof], val);
        }
        
        ef.DOFVals.set(dof, val);
    }

    public double getValForDOF(int dof, double val){
        setDOF(dof,val);
        return ef.getEnergyForDOF(dof);
    }


    public double getInitStepSize(int dof){
        //initial step size to use for dof
        int DOFNum = ef.DOFNums[dof];
        if(DOFNum<0)//sin or cosine...
            return 0.05;
        else {
            DegreeOfFreedom curDOF = globalDOFList[DOFNum];
            return curDOF.getMeshWidth();
        }
    }

    
    @Override
    public boolean isDOFAngle(int dof){
        int DOFNum = ef.DOFNums[dof];
        if(DOFNum<0)//sine or cosine can't be an angle
            return false;
        else
            return globalDOFList[DOFNum].isDOFAngle();
    }
     
    
    public void printTerms(DoubleMatrix1D x){
        //print the values of the individual LSB terms, all on a line
        setDOFs(x);
        ef.printTerms();
    }


}
