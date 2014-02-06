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
//	SparseVDWEnergy.java
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
import java.util.HashMap;
import java.io.Serializable;


public class SparseVDWEnergy extends EnergyFunction implements Serializable {
    //This class represents a couple of VDW energies for a molecule
    
    Molecule m;
    int atomNums[];
    double coeffs[];//the corresponding coefficients
    
    //if we need to call this on other molecules, we'll need to call by residue # and atom # within the residue
    //since we can only count on the residues we call to be of the right AA type
    int resNums[];//molecule-based residue numbers
    //int resAtomNums[];//residue-based atom numbers
    String atomNames[];
    
    //each term corresponds to two indices in atomNums and coeffs:
    //the two moleculeAtomNumbers in atomNums; the r^-12 and then r^-6 coefficients in coeffs 
    int numTerms;
    
    DoubleMatrix1D DOFVals;//PROBABLY SHARED WITH AN OBJECTIVE FUNCTION
    double distCutoff;
    
    double V6, V12;//VDW multipliers for attractive and repulsive terms
        
    //TODO: EXPLICITLY GET ATOM LOCATIONS FROM DOFS INSTEAD OF HAVING TO MESS WITH MOLECULE/
    //MOVE ALL ATOMS AFFECTED BY THEM
    
    static boolean includeElec = true;//include the electric along with the VDW energy term
    double elecCoeffs[];//coefficients for this
    boolean distDepDielec = false;
    double coulombFactor = 0;
    
    
    //SparseVDWEnergy sve = new SparseVDWEnergy(of.m,of.efunc,explicitVDWCutoff,of.curDOFVals);
   public SparseVDWEnergy(ContSCObjFunction of,double explicitVDWCutoff,DoubleMatrix1D[] sampAbs){
       init(of.m,of.efunc,explicitVDWCutoff,of.curDOFVals,of,sampAbs);
   }
    
    
    public SparseVDWEnergy(Molecule molec, EnergyFunction ef, double cutoff, DoubleMatrix1D dv){
        init(molec,ef,cutoff,dv,null,null);
    }
    
    void init(Molecule molec, EnergyFunction ef, double cutoff, DoubleMatrix1D dv,
            ContSCObjFunction of, DoubleMatrix1D[] samples){
        //Take the VDW interactions with distances below distCutoff
        //if of and samples aren't null then count distances below distCutoff in any of the samples
        //(from the DOFs of of)
        //if they are null then just check distances in the current geometry
        
        DOFVals = dv;
        m = molec;
        distCutoff = cutoff;
        
        double VDWMultiplier = ef.getAmber96ext().vdwMultiplier;
        V6 = VDWMultiplier*VDWMultiplier;//square then cube
        V6=V6*V6*V6;
        V12 = V6*V6;
        
        if(includeElec){
            Amber96ext a96ff = ef.getAmber96ext();//to get general parameters
            distDepDielec = a96ff.distDepDielect;
            coulombFactor = a96ff.constCoulomb / a96ff.dielectric;
        }
        
                
        //now traverse ef to get all VDW terms with specified distance, and their coefficients in ef
        HashMap<Integer,double[]> terms = new HashMap<Integer,double[]>();
        //will record pair of molecule atom numbers (as (lower number)*m.numberOfAtoms+(higherNumber) )
        //mapped to the pair (coefficient for r^-12, coefficient for r^-6)
        
        if(of==null){
            getTerms(ef,terms,1);
        }
        else{
            for(DoubleMatrix1D samp : samples){
                of.setDOFs(samp);
                HashMap<Integer,double[]> sampTerms = new HashMap<Integer,double[]>();
                getTerms(ef,sampTerms,1);
                //merge sampTerms into terms
                for( int pair : sampTerms.keySet() ){
                    terms.put(pair, sampTerms.get(pair));
                    //coefficients will be the same for each of the samples
                    //since these come from the energy function, not geometry
                }
            }
        }
        
        //now fill in atomNums and coeffs accordingly
        numTerms = terms.size();
        atomNums = new int[2*numTerms];
        coeffs = new double[2*numTerms];
        
        resNums = new int[2*numTerms];
        
        //resAtomNums = new int[2*numTerms];
        atomNames = new String[2*numTerms];
        
        if(includeElec)
            elecCoeffs = new double[numTerms];
            
        int count = 0;
        for( int index : terms.keySet() ){
            int atom1 = index / m.numberOfAtoms;
            int atom2 = index % m.numberOfAtoms;
            double termCoeffs[] = terms.get(index);
            
            atomNums[2*count] = atom1;
            atomNums[2*count+1] = atom2;
            coeffs[2*count] = termCoeffs[0];
            coeffs[2*count+1] = termCoeffs[1];
            
            if(includeElec)
                elecCoeffs[count] = termCoeffs[2];
            
            getResBasedNums(2*count);
            getResBasedNums(2*count+1);
            
            count++;
        }
        
        System.out.println("Created SVE with "+numTerms+" terms");
    }
    
    
    void getResBasedNums(int index){
        //fill in the residue # and residue-based atom number 
        //once the molecule-based atom number has been put into atomNums
        Atom at = m.atom[atomNums[index]];
        resNums[index] = at.moleculeResidueNumber;
        //resAtomNums[index] = at.residueAtomNumber;
        atomNames[index] = at.name;
    }
    
    
    void getTerms(EnergyFunction ef,HashMap<Integer,double[]> termMap, double coeff){
        //get the appropriate VDW terms from ef, which has coefficient coeff in the energy,
        //and put them in termMap
        
        if(ef instanceof ForceFieldEnergy){//The terms are only in ForceFieldEnergies
            
            double nbTerms[] = ((ForceFieldEnergy)ef).a96ff.nonBondedTerms;
            RotMatrix r = new RotMatrix();
            
            for(int t=0; t<nbTerms.length/4; t++){
                
                int atom1 = (int)nbTerms[4*t];
                int atom2 = (int)nbTerms[4*t+1];
                
                double dist = r.norm( r.subtract( m.getActualCoord(atom1), m.getActualCoord(atom2) ) );
                
                if( dist<distCutoff){
                    
                    int index = Math.min(atom1,atom2)*m.numberOfAtoms + Math.max(atom1,atom2);
                    double coeff12 = coeff*V12*nbTerms[4*t+2];//coefficient for this term for r^-12
                    double coeff6 = coeff*V6*nbTerms[4*t+3];
                                        
                    
                    double coeffElec = coeff*m.atom[atom1].charge*m.atom[atom2].charge*coulombFactor;
                    
                    if(termMap.containsKey(index)){
                        double[] termCoeffs = termMap.get(index);
                        termCoeffs[0] += coeff12;
                        termCoeffs[1] += coeff6;
                        termCoeffs[2] += coeffElec;
                        if(termCoeffs[0]==0 && termCoeffs[1]==0 && termCoeffs[2]==0)//terms cancelled out
                            termMap.remove(index);
                    }
                    else{
                        double termCoeffs[] = new double[] {coeff12,coeff6,coeffElec};
                        termMap.put(index, termCoeffs);
                    }
                } 
            }
        }
        else if(ef instanceof MultiTermEnergyFunction){//recursively look for ForceFieldEnergies
            MultiTermEnergyFunction mtef = (MultiTermEnergyFunction)ef;
            for(int t=0; t<mtef.numTerms; t++)
                getTerms(mtef.terms[t],termMap,mtef.coeffs[t]*coeff);
        }
        //else we have an EnergyFunction not containing VDW, so do nothing
    }
    
    
    @Override
    public double getEnergy(){
        
        RotMatrix r = new RotMatrix();
        double ans = 0;
    
        for(int t=0; t<numTerms; t++){
            double dist2 = r.normsq( r.subtract( m.getActualCoord(atomNums[2*t]), m.getActualCoord(atomNums[2*t+1]) ) );
            ans += termValue(t,dist2);
        }
        
        return ans;
    }
    
    
    private double termValue(int t, double dist2){
        //evaluate term t, given the squared distance between the atoms involved
        double dist6 = dist2*dist2*dist2;
        double dist12 = dist6*dist6;
        double ans = coeffs[2*t]/dist12 - coeffs[2*t+1]/dist6;

        if(includeElec){
            if(distDepDielec)
                ans += elecCoeffs[t]/dist2;
            else
                ans += elecCoeffs[t]/Math.sqrt(dist2);
        }
        
        return ans;
    }
    
    
    public double getEnergy(Molecule otherMolec){
        //like getEnergy() but evaluate using a different molecule
        
        
        RotMatrix r = new RotMatrix();
        double ans = 0;
    
        for(int t=0; t<numTerms; t++){
            //int atomNum1 = otherMolec.residue[resNums[2*t]].atom[resAtomNums[2*t]].moleculeAtomNumber;
            //int atomNum2 = otherMolec.residue[resNums[2*t+1]].atom[resAtomNums[2*t+1]].moleculeAtomNumber;
            int atomNum1 = otherMolec.residue[resNums[2*t]].getAtomNameToMolnum(atomNames[2*t]);
            int atomNum2 = otherMolec.residue[resNums[2*t+1]].getAtomNameToMolnum(atomNames[2*t+1]);
            
            double dist2 = r.normsq( r.subtract( otherMolec.getActualCoord(atomNum1), otherMolec.getActualCoord(atomNum2) ) );
            ans += termValue(t,dist2);
        }
        
        return ans;
    }

    
    //Energy evaluated is already pretty minimal...not expected to save anything with partial computation
    public double getEnergy(int a){
        return getEnergy();
    }

    public void setupPartialComputation(int[][] res){

    }
    
    
    
    
}
