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
//	CETEnergyFunction.java
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
import java.util.ArrayList;
import java.util.HashMap;


public class CETEnergyFunction extends EnergyFunction {

    DoubleMatrix1D DOFVals;//Values of the degrees of freedom
    //This object can be shared with the minimizer

    int DOFNums[];//DOFNums of DOFs in DOFVals

    HashMap<Integer,Integer> revDOFNums = new HashMap<Integer,Integer>();//revDOFNums[DOFNum] is the index of DOFNum in DOFNums
    HashMap<Integer,Double> discreteDOFs = new HashMap<Integer,Double>();//discrete DOFs (DOFNum mapped to value): fixed for all terms
    
    int numDOFs;

    double constTerm = 0;//MAYBE TO BE USED LATER

    ContETerm terms[];
    //These are terms dependent on values in DOFVals
    //that are summed together to get the total energy
    //each term depends on a subset of those values (specified by term.DOFNums);
    int numTerms;

    boolean includeMinE;//include the voxel minimum

    double termCoeffs[] = null;//coefficients for terms (default = all 1, if termCoeffs stays null)
    
    boolean termHasDOF[][];//which terms (numbered as in terms) have which DOFs (numbered as in DOFVals)
    
    
    public CETEnergyFunction(DoubleMatrix1D x, int DOFList[], ContETerm termArr[], boolean ime){

        DOFVals = x;
        DOFNums = DOFList;
        numDOFs = DOFVals.size();
        terms = termArr;
        numTerms = terms.length;

        includeMinE = ime;
        
        for( int a=0; a<DOFNums.length; a++ )
            revDOFNums.put(DOFNums[a], a);
        
        
        //store discrete DOFs (fixed for all the terms)
        for(ContETerm term : termArr){
            if(term!=null){
                if(term.discreteDOFNums!=null){
                    for(int ddof=0; ddof<term.discreteDOFNums.size(); ddof++){
                        int DOFNum = term.discreteDOFNums.get(ddof);
                        double val = term.discreteDOFVals.get(ddof);
                        if(discreteDOFs.containsKey(DOFNum)){
                            if(discreteDOFs.get(DOFNum)!=val)
                                System.err.println("ERROR: discrete DOF # "+DOFNum+" has >1 different values for different terms!  "
                                    +val+" "+discreteDOFs.get(DOFNum));
                        }
                        else
                            discreteDOFs.put(DOFNum, val);
                    }
                }
            }
        }
        
        
        termHasDOF = new boolean[numTerms][numDOFs];
        for(int t=0; t<numTerms; t++){
            if(terms[t]!=null){
                for(int DOFNum : terms[t].DOFNums){
                    termHasDOF[t][revDOFNums.get(DOFNum)] = true;
                }
            }
        }
        
    }


    //get total energy
    public double getEnergy(){

        return getEnergyForDOF(-1);
    }
    
    
    
    //Get the values of DOFs needed for a particular term
    //We need the term's degrees of freedom; the term itself can be null
    //unless we might have to convert between angles and sines/cosines
    DoubleMatrix1D getTermDOFVals(DoubleMatrix1D overallVals, int termDOFs[], ContETerm term){
        
        
        DoubleMatrix1D termDOFVals = DoubleFactory1D.dense.make(termDOFs.length);
        for(int dof=0; dof<termDOFs.length; dof++){
            
            if (revDOFNums.containsKey(termDOFs[dof])) {
                int valIndex = revDOFNums.get(termDOFs[dof]);//revDOFNums[termDOFs[dof]];
                termDOFVals.set(dof,overallVals.get(valIndex));
            }
            else if(discreteDOFs.containsKey(termDOFs[dof])){//it might be a fixed discrete DOF...
                termDOFVals.set(dof,discreteDOFs.get(termDOFs[dof]));
            }
            else {
                System.err.println("ERROR: can't find DOF # "+termDOFs[dof]);
                new Exception().printStackTrace();
                System.exit(1);
            }
        }
        
        return termDOFVals;
    }
    
    
    
    
    //Get energy terms involving a particular DOF (numbered as in DOFVals, etc.)
    //if -1, then get total energy
    //note that this is different from other energy functions because we do it by res and we don't setupPartialComputation
    //hence why we aren't overriding getEnergy(int a)
    public double getEnergyForDOF(int dof){
        double E = constTerm;

        for(int t=0; t<terms.length; t++){
            
            if(dof!=-1){//if not total energy, only take terms involving this DOF
                if(!termHasDOF[t][dof])
                    continue;
            }
            
            ContETerm term  = terms[t];
            
            if(term!=null){
                DoubleMatrix1D termDOFVals = getTermDOFVals(DOFVals,term.DOFNums,term);
                double termVal = term.evaluate(termDOFVals,includeMinE);
                
                                
                if(termCoeffs==null){
                    E += termVal;
                }
                else{
                    E += termVal*termCoeffs[t];
                }
                
                
                if(Double.isNaN(E)){
                    System.err.println("ERROR: CETEnergyFunction returning NaN.  Writing ContETerm to badTerm.dat");
                    System.err.println("termVal: "+termVal+"  termDOFVals: "+termDOFVals);
                    System.err.println("DOFVals: "+DOFVals+" DOFNums: ");
                    for(int DOFNum : DOFNums)
                        System.err.println(DOFNum);
                    KSParser.outputObject(term, "badTerm.dat");
                    System.exit(0);
                }   
            }
        }

        return E;
    }
    
    
    public double getEnergy(int a){
        return getEnergy();
    }

    public void setupPartialComputation(int[][] res){

    }
    
    
    public void printTerms(){
        checkTerms(true);
    }
    
    public double maxTerm(){
        return checkTerms(false);
    }
    
    
    public double checkTerms(boolean print){
        //go through the values of the individual CETs, printing all on a line (if print)
        //constant terms not included
        //return maximum term
        //this is mostly to check out what bCutoffs we might need, etc.
        //based on getEnergy()

        double maxTerm = 0;

        for(int t=0; t<terms.length; t++){
            ContETerm term  = terms[t];
            
            if(term==null){
                if(print)
                    System.out.print("null ");
            }
            else {
                DoubleMatrix1D termDOFVals = getTermDOFVals(DOFVals,term.DOFNums,term);
                double termVal = term.evaluate(termDOFVals,false);//minE constant terms not included
                
                
                maxTerm = Math.max(maxTerm,termVal);
                
                if(print)
                    System.out.print(termVal+" ");
            }
        }
        
        return maxTerm;
    }

}
