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
//	EPoly.java
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
import cern.jet.math.Functions;
import java.io.FileInputStream;
import java.util.ArrayList;


public class EPoly extends ContETerm {
    //energy term represented by a polynomial fit


    double coeffs[];//coefficients for series (expanded in original, centered coordinates)
    int order;//order of polynomial (can be 3 or 4)

    
    //Sparse VDW terms
    SparseVDWEnergy sve = null;
    ContSCObjFunction sveOF = null;//To set DOFs for SVE
    double baseSVE = 0;//value of SVE terms at center (SVE will be evaluated relative to this)
    
    
    

    public EPoly(int[] DOFNums, int numDOFs, DoubleMatrix1D DOFmax, DoubleMatrix1D DOFmin, 
            DoubleMatrix1D center, double minE, double[] series, int ord ) {
        
        super(DOFNums,numDOFs,DOFmax,DOFmin,center,minE);
        coeffs = series;
        order = ord;
    }

    
    
    @Override
    public double evaluate(DoubleMatrix1D x, boolean includeMinE, Molecule m) {
        //evaluate this EPoly as a function of internal coordinates x
         
        DoubleMatrix1D z = toRelCoords(x);
        double serVal = evalSeries(z);
        
        if(includeMinE)
            serVal += minE;
                 
        if(sve!=null){
            if(m==null){//do SVE with the molecule stored in SVE
                sveOF.setDOFs(x); 
                return serVal + sve.getEnergy() - baseSVE;
            }
            else {//m provided, and should have DOFs already set to x    
                
                
                //DEBUG!!!  To make sure results using an external molecule or the stored molecule agree
                //use with the corresponding DEBUG block in CETObjFunction.getConstraints
                /*
                sveOF.setDOFs(x);
                double check1 = sve.getEnergy();
                double check2 = sve.getEnergy(m);
                
                if(Math.abs(check1-check2)/Math.max(Math.abs(check1),Math.abs(check2))>0.01 && Math.abs(check1-check2)>0.01){
s
                    sve.getEnergy();
                    sve.getEnergy(m);
                    sve.getEnergy();
                    sve.getEnergy(m);
                    sve.getEnergy();
                    sve.getEnergy(m);
                    sve.getEnergy();
                    sve.getEnergy(m);
                    sve.getEnergy();
                    sve.getEnergy(m);
                    sve.getEnergy();
                    sve.getEnergy(m);
                    sve.getEnergy();
                    sve.getEnergy(m);
                    
                    
                    
                    
                     //DEBUG!!!
                    new KSParser().saveMolecule(m, "genm.pdb", 0);
                    for(int res=0; res<m.numberOfResidues; res++){//can't output sve.m directly but can do this
                        if(m.residue[res].name.equalsIgnoreCase(sve.m.residue[res].name)){
                            if(sve.m.residue[res].atom[0]!=null){
                                for(int at=0; at<m.residue[res].numberOfAtoms; at++){
                                    m.residue[res].atom[at].coord = sve.m.getActualCoord( sve.m.residue[res].getAtomNameToMolnum(m.residue[res].atom[at].name) );
                                    m.updateCoordinates(m.residue[res].atom[at]);
                                }
                            }
                        }
                    }
                    new KSParser().saveMolecule(m, "termm.pdb", 0);
                    System.exit(0);
                }*/
                
                //DEBUG!!!!
                //new KSParser().saveMolecule(m, "checkm"+System.currentTimeMillis()+".pdb",0);
          
                return serVal + sve.getEnergy(m) - baseSVE;
            }
        }
        else 
            return serVal;
    }
    
    
    
    double evalSeries(DoubleMatrix1D z){
        //evaluate the actual series
        //(function of relative coordinates)
        return SeriesFitter.evalSeries(coeffs, z, numDOFs, false, order);
    }
    
    
}
