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
//	MultiTermEnergyFunction.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
import java.util.Arrays;


//This an energy function consisting of terms that are other energy functions
//For example, it can be AMBER energy + dihedral energy
//or ( energy of a pair of residues ) - (energy of residue 1) - (energy of residue 2)
//Total energy = sum_i (coeff[i] * energy i)

public class MultiTermEnergyFunction extends EnergyFunction {
    
    EnergyFunction terms[];
    double coeffs[];

    int numTerms;

    //Constructor with no terms
    public MultiTermEnergyFunction(){
        numTerms = 0;
        terms = new EnergyFunction[0];
        coeffs = new double[0];
    }

    //General constructor
    public MultiTermEnergyFunction( EnergyFunction terms[], double coeffs[] ){
        this.terms = terms;
        this.coeffs = coeffs;
        numTerms = terms.length;
    }


    //Energy = unweighted sum of terms
    public MultiTermEnergyFunction( EnergyFunction posTerms[] ){
        terms = posTerms;
        numTerms = terms.length;
        coeffs = new double[numTerms];
        Arrays.fill(coeffs, 1);
    }

    //Energy = (sum of terms in plus) - (sum of terms in minus)
    public MultiTermEnergyFunction(EnergyFunction[] plus, EnergyFunction[] minus){

        numTerms = plus.length + minus.length;

        terms = new EnergyFunction[numTerms];

        System.arraycopy(plus, 0, terms, 0, plus.length);
        System.arraycopy(minus, 0, terms, plus.length, minus.length);

        coeffs = new double[numTerms];
        Arrays.fill(coeffs, 0, plus.length, 1);//Positive coefficients
        Arrays.fill(coeffs, plus.length, numTerms, -1);//Negative coefficients
    }


  @Override
    public EnergyFunction addTerm(EnergyFunction ef){
        EnergyFunction newTerms[] = new EnergyFunction[numTerms+1];
        double newCoeffs[] = new double[numTerms+1];
        System.arraycopy(terms,0,newTerms,0,numTerms);
        System.arraycopy(coeffs,0,newCoeffs,0,numTerms);
        newTerms[numTerms] = ef;
        newCoeffs[numTerms] = 1;
        return new MultiTermEnergyFunction(newTerms, newCoeffs);
    }



    //add a term to this with a given coefficient
    public void addTermWithCoeff(EnergyFunction ef, double coeff){
        
        EnergyFunction newTerms[] = new EnergyFunction[numTerms+1];
        double newCoeffs[] = new double[numTerms+1];
        System.arraycopy(terms,0,newTerms,0,numTerms);
        System.arraycopy(coeffs,0,newCoeffs,0,numTerms);
        newTerms[numTerms] = ef;
        newCoeffs[numTerms] = coeff;

        terms = newTerms;
        coeffs = newCoeffs;
        numTerms++;
    }
    

    public double getEnergy(){
        
        double E = 0;
        
        for(int a=0;a<numTerms; a++)
            E += coeffs[a] * terms[a].getEnergy();

        if(Double.isNaN(E) || Double.isInfinite(E))//This can happen if there are positive and negative terms
            //with infinite energy...we assume this to be an impossible conformation
            //and thus return inifinity
            return Double.POSITIVE_INFINITY;

        return E;
    }



    //Partial computations
    public void setupPartialComputation(int residues[][]){
        for( EnergyFunction ef : terms )
            ef.setupPartialComputation(residues);
    }

    public double getEnergy(int part){
        double E = 0;

        for(int a=0;a<numTerms; a++)
            E += coeffs[a] * terms[a].getEnergy(part);

        if(Double.isNaN(E))
            return Double.POSITIVE_INFINITY;

        return E;
    }


    //This returns the Amber96ext from some ForceFieldEnergy term
    //for purposes of getting parameters
  @Override
    public Amber96ext getAmber96ext(){
        
        for(EnergyFunction ef : terms){
            if(ef instanceof ForceFieldEnergy)
                return ((ForceFieldEnergy)ef).getAmber96ext();
        }

        System.err.println("ERROR: No Amber96ext found in multi-term energy function, but"+
                "getAmber96ext() was called.");
        System.exit(1);
        return null;
    }
  
  
  
  public void removeLastTerm(){
        EnergyFunction newTerms[] = new EnergyFunction[numTerms-1];
        double newCoeffs[] = new double[numTerms-1];
        System.arraycopy(terms,0,newTerms,0,numTerms-1);
        System.arraycopy(coeffs,0,newCoeffs,0,numTerms-1);
        terms = newTerms;
        coeffs = newCoeffs;
        numTerms--;
  }
    

}
