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
//	LoopClosureAdjustment.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/* This perturbation involves finding alternate ways to close loops
 * A plain "loop closure adjustment" does this by finding alternate sets of dihedrals
 * for a three-residue segment without altering the chain elsewhere
 * This is a discrete type of perturbation

 * There are also some specialized versions:
 * 1. A helix or sheet reduction is a loop closure adjustment containing one or more non-loop residues
 * (which will be removed from their helix or sheet)
 * 2. A helix or sheet expansion is formed by adding a loop residue to an adjacent helix or sheet
 * and then trying to close the loop
 * 3. A partial structure switch uses experimental data (maybe its own class??)
 *
 */
import java.util.HashMap;
import java.util.ArrayList;


public class LoopClosureAdjustment extends Perturbation {

    

    double predParams[][] = new double[0][];//Perturbation parameters for predecessors.
    //Indices: 1. Which set of predecessor parameters this is (state index) 2. Which perturbation

    double[][][][][] rotMatrices = new double[0][][][][];//Rotation matrices
    //indices: 1. State index (like predParams) 2. # of solution 3. number of matrix
    //(0 for the middle sidechain & CA & everything before it; 1 for everything after) 4-5. Matrix indices

    //Should the first and/or last residue copy the dihedrals (i.e. the secondary structure) of the one
    //before or after it respectively?  
    boolean addSSStart = false;
    boolean addSSEnd = false;
    //Possibly also add the option of some out-of-plane motion before and/or after?
    //(Added the same way)

    TripeptideClosure tc;//This object stores the bond lengths & angles and omega dihedrals
    //and will be used to compute the alternate conformations

//Partial structure switch: here or separate?  Not sure

    public LoopClosureAdjustment(String t, Molecule molec, int resList[]){

        type=t;
        m=molec;
        resDirectlyAffected=resList;

        if(type.equalsIgnoreCase("SSNE"))//Secondary structure N-terminal expansion
            addSSEnd = true;
        else if(type.equalsIgnoreCase("SSCE"))//Secondary structure C-terminal expansion
            addSSStart = true;

        explicitUnperturbed = !(addSSStart||addSSEnd);

        tc = new TripeptideClosure(m, getTripepRes());
        //This grabs the bond lengths, angles, and omega dihedrals before we start changing them
    }

    public void setDefaultParams(){//There can be up to 16 solutions so we create that many states
        //They're all discrete though
        minParams = new double[16];
        maxParams = new double[16];
        for(int a=0;a<16;a++)
            minParams[a] = maxParams[a] = a;
    }



    public boolean doPerturbationMotion(double param){

        int predState = -1;//If the predecessors have been in this state before this will indicate which set of known parameters this is

        for( int a=0; ( a<predParams.length ) && ( predState == -1 ); a++ ){//Look for a match in predParams to the current predecessor parameters

            predState = a;

            for( int b=0; b<predecessors.length; b++ ){
                if( predParams[predState][b] != m.perts[predecessors[b]].curParam ){//Predecessor parameters
                    predState = -1;
                    break;
                }
            }

        }

        if(predState == -1){//We need to calculate the solutions
            calcSolns();
            predState = predParams.length - 1;
        }


        if( (param < 0) || param >= ( rotMatrices[predState].length) ){
            invalidState = true;//The parameter value is invalid for this state of the predecessors
            m.validBB = false;//So flag the molecule to have a bad backbone
            return false;//So we fail to apply the perturbation
        }

        applySSChanges();

        double[][][] rm = rotMatrices[predState][(int)param];

        if(addSSStart){
            int tripepRes[] = {resDirectlyAffected[1], resDirectlyAffected[2], resDirectlyAffected[3]};
            applyTripeptideRot(rm, tripepRes);
        }
        else
            applyTripeptideRot(rm, resDirectlyAffected);

        invalidState = false;
        m.updateValidBB();//See if we can lift the bad BB flag
        return true;
    }





    public void applyTripeptideRot(double rm[][][], int res[]){
        //Perform the loop closure adjustment rotations on three consecutive residues,
        //given the rotation matrices and the molecule residue numbers


        //All the rotation are about anchor-point CAs, so order is not important here like for the shear
        //The second peptide plane is rotated about the ending anchor CA; the other stuff is rotated about the starting one
        int CA0AtNum = m.residue[res[0]].getAtomNameToMolnum("CA");
        int CA2AtNum = m.residue[res[2]].getAtomNameToMolnum("CA");

        
        //We rotate the middle residue except the carbonyl, plus the carbonyl group of the first, with the 1st matrix
        int pep1sc[] = concatenateArrays( m.residue[res[0]].getAtomList(false, false, false, true),
                m.residue[res[1]].getAtomList(true, true, true, false) );
        m.rotateAtomList( pep1sc, rm[0], m.actualCoordinates[3*CA0AtNum],
                m.actualCoordinates[3*CA0AtNum+1], m.actualCoordinates[3*CA0AtNum+2], false);


        //Rotate the amide group of the last residue and carbonyl of the middle residue
        //with the second matrix
        int pep2[] = concatenateArrays( m.residue[res[1]].getAtomList(false, false, false, true),
                m.residue[res[2]].getAtomList(true,false,false,false) );
        m.rotateAtomList( pep2, rm[1], m.actualCoordinates[3*CA2AtNum],
                m.actualCoordinates[3*CA2AtNum+1], m.actualCoordinates[3*CA2AtNum+2], false);

    }



    public int[] getTripepRes(){
        //Identify the residues in the tripeptide to be closed
        if(addSSStart || addSSEnd){
            int[] tripepRes = new int[3];

            if(addSSStart)
                System.arraycopy(resDirectlyAffected, 1, tripepRes, 0, 3);
            else
                System.arraycopy(resDirectlyAffected, 0, tripepRes, 0, 3);

            return tripepRes;
        }
        else
            return resDirectlyAffected;
    }



    public int calcSolns(){//Solve the loop closure equations given the current state of predecessor perturbations
        //Store conformational change info and return the number of solutions

        int tripepRes[] = getTripepRes();//Residues in the tripeptide to be closed



        //Get necessary lengths, angles, dihedrals

        //MAYBE TEMPORARY CODE!!!
        //DON'T CHANGE TC IF GOT FREESC PERTURBED STATE UNDER IT
        boolean updateTC = true;


        if(updateTC)
            tc = new TripeptideClosure(m, getTripepRes());


        RotMatrix rm = new RotMatrix();
        double[][][][] matrices = null;


        double[] genChi1 = null;

        //Backup coordinates for any residues to be moved (for secondary-structure adjustment purposes)
        //Would need to do the other same for other motions outside the closed tripeptide
        //We use oldCoords, like in undo; since we have not applied the perturbation yet this does not interfere with undoing
        if(addSSStart || addSSEnd){
            oldCoords = new ArrayList<HashMap<String, double[]>>(resDirectlyAffected.length);
            genChi1 = getGenChi1();//Note the generalized chi1 for each residue so we can restore it after the changes

            for(int a=0;a<resDirectlyAffected.length;a++)
                oldCoords.add(a, storeResBB(m.residue[resDirectlyAffected[a]]) );
        }

        applySSChanges();

        

        double r_soln_n[][][] = new double[16][3][3];
        double r_soln_a[][][] = new double[16][3][3];
        double r_soln_c[][][] = new double[16][3][3];

        double firstN[] = m.getActualCoord( m.residue[ tripepRes[0] ].getAtomNameToMolnum("N") );//Coordinates of the first residue's N
        double firstCA[] = m.getActualCoord( m.residue[ tripepRes[0] ].getAtomNameToMolnum("CA") );//Its CA
        double firstC[] = m.getActualCoord( m.residue[ tripepRes[0] ].getAtomNameToMolnum("C") );

        double midN[] = m.getActualCoord( m.residue[ tripepRes[1] ].getAtomNameToMolnum("N") );//Starting coordinates of the middle N
        double midCA[] = m.getActualCoord( m.residue[ tripepRes[1] ].getAtomNameToMolnum("CA") );
        double midC[] = m.getActualCoord( m.residue[ tripepRes[1] ].getAtomNameToMolnum("C") );

        double lastN[] = m.getActualCoord( m.residue[ tripepRes[2] ].getAtomNameToMolnum("N") );
        double lastCA[] = m.getActualCoord( m.residue[ tripepRes[2] ].getAtomNameToMolnum("CA") );
        double lastC[] = m.getActualCoord( m.residue[ tripepRes[2] ].getAtomNameToMolnum("C") );


        int numSoln = tc.solve_3pep_poly( firstN, firstCA, lastCA, lastC, r_soln_n, r_soln_a, r_soln_c);


        int shift = 0;
        if( addSSStart || addSSEnd )
            shift = 1;
        //This is to make room for the unperturbed state at matrices[0]

        matrices = new double[numSoln+shift][2][][];

        int unperturbed = -1;//The unperturbed state: will be found by least-squares comparison to the original state
        double lowestSum = Double.POSITIVE_INFINITY;

        
        for(int s=0; s<numSoln; s++){

            //First rotation: Based on midCA and midN
            matrices[s+shift][0] = rm.getSuperposingRotMatrix( rm.subtract( midCA, firstCA ), rm.subtract( r_soln_a[s][1], firstCA),
                    rm.subtract( midN, firstCA ), rm.subtract( r_soln_n[s][1], firstCA) );

            //This rotation is about the last CA instead of the first one
            matrices[s+shift][1] = rm.getSuperposingRotMatrix( rm.subtract( midC, lastCA ), rm.subtract( r_soln_c[s][1], lastCA),
                    rm.subtract( lastN, lastCA ), rm.subtract( r_soln_n[s][2], lastCA) );


            if( ! (addSSStart || addSSEnd) ){//See if this might be the unperturbed state
                double checkSum = rm.normsq( rm.subtract(firstC, r_soln_c[s][0]) )
                        + rm.normsq( rm.subtract(midN, r_soln_n[s][1]) )
                        + rm.normsq( rm.subtract(midCA, r_soln_a[s][1]) )
                        + rm.normsq( rm.subtract(midC, r_soln_c[s][1]) )
                        + rm.normsq( rm.subtract(lastN, r_soln_n[s][2]) );

                if(checkSum < lowestSum){
                    lowestSum = checkSum;
                    unperturbed = s;
                }
            }

        }

        if(addSSStart || addSSEnd){//This is like undo.  In this case we do not expect the unperturbed state to have been generated
            for(int a=0;a<resDirectlyAffected.length;a++)
                restoreResBB( m.residue[resDirectlyAffected[a]], oldCoords.get(a) );

            setGenChi1(genChi1);
        }
        else{

            /*if(numSoln == 0){//At least the unperturbed state should have been generated
                System.err.println("Error applying tripeptide closure to residues " + resDirectlyAffected[0] + " through " + resDirectlyAffected[2]);
                System.exit(1);
            }

            //The unperturbed state needs to have parameter 0, so switch it to that position
            //We don't need matrices for the unperturbed state so those can be overwritten
            matrices[unperturbed] = matrices[0];
             *
             */

            //With peptide translations/rotations preceding this LCA there might be cases where there is no unperturbed state...
            //this will simply lead to invalidState being set true as long as the predecessors are leading to this
            if(numSoln > 0){//Switch the unperturbed state into parameter 0
                double unpertMatrices[][][] = matrices[unperturbed];
                matrices[unperturbed] = matrices[0];
                matrices[0] = unpertMatrices;//We keep matrices for the unperturbed state since we might need them to close the BB chain
            }


        }


        //Store the matrices and current perturbation parameters at the end of rotMatrices and predParams respectively
        double newPredParams[][] = new double[predParams.length+1][];
        double newRotMatrices[][][][][] = new double[rotMatrices.length+1][][][][];
        System.arraycopy(predParams,0,newPredParams,0,predParams.length);
        System.arraycopy(rotMatrices,0,newRotMatrices,0,rotMatrices.length);
        predParams = newPredParams;
        rotMatrices = newRotMatrices;

        predParams[predParams.length-1] = new double[predecessors.length];
        for(int b=0;b<predecessors.length;b++)
            predParams[predParams.length-1][b] = m.perts[predecessors[b]].curParam;

        rotMatrices[rotMatrices.length-1] = matrices;

        return numSoln;

    }



    public void applySSChanges(){
        //Add loop residues to adjacent secondary structures by copying the dihedrals of the nearest residues in those structures

        if( ! (addSSStart || addSSEnd) )
            return;//Nothing to do

        RamachandranChecker rcheck = RamachandranChecker.getInstance();
        RotMatrix rm = new RotMatrix();

        int closeStart=0;//Start of tripeptide to close

        double anchor1Old[] = null, anchor2Old[] = null;//Positions of anchor CAs before SS changes

        if(addSSStart){//Add the first residue

            Residue res0 = m.residue[resDirectlyAffected[0]];
            if(!m.checkNBonded(resDirectlyAffected[0])){
                System.err.println("ERROR: Trying to add " + res0.fullName +
                        " to a preceding secondary structure, but it is not bonded to any preceding residue");
                System.exit(1);
            }

            Residue res1 = m.residue[resDirectlyAffected[1]];

            anchor1Old = m.getActualCoord( res1.getAtomNameToMolnum("CA") );

            double phiPsi[] = rcheck.getPhiPsi(m, resDirectlyAffected[0]-1);//phi, psi for residue before perturbation starts

            //First set phi
            int atomList1[] = res0.getAtomList(false, true, true, true);
            int atomList2[] = res1.getAtomList(true, true, true, false);
            int[] atomList = concatenateArrays(atomList1, atomList2);

            m.setTorsion( m.residue[resDirectlyAffected[0]-1].getAtomNameToMolnum("C") ,
                    res0.getAtomNameToMolnum("N"), res0.getAtomNameToMolnum("CA"),
                    res0.getAtomNameToMolnum("C"), phiPsi[0], atomList, atomList.length, false);
            
            
            //Now psi
            atomList1 = res0.getAtomList(false, false, false, true);
            atomList2 = res1.getAtomList(true, true, true, false);
            atomList = concatenateArrays(atomList1, atomList2);

            m.setTorsion( res0.getAtomNameToMolnum("N"), res0.getAtomNameToMolnum("CA"),
                    res0.getAtomNameToMolnum("C"), res1.getAtomNameToMolnum("N"),
                    phiPsi[1], atomList, atomList.length, false);

            closeStart=1;
        }

        if(addSSEnd){
            
            Residue resEnd = m.residue[resDirectlyAffected[resDirectlyAffected.length-1]];
            if(!m.checkCBonded(resDirectlyAffected[resDirectlyAffected.length-1])){
                System.err.println("ERROR: Trying to add " + resEnd.fullName +
                        " to an ensuing secondary structure, but it is not bonded to any ensuing residue");
                System.exit(1);
            }

            Residue resPrior = m.residue[resDirectlyAffected[resDirectlyAffected.length-2]];

            anchor2Old = m.getActualCoord( resPrior.getAtomNameToMolnum("CA") );

            double phiPsi[] = rcheck.getPhiPsi(m, resDirectlyAffected[resDirectlyAffected.length-1]+1 );
            //phi, psi for residue after perturbation ends

            //First set psi
            int atomList1[] = resEnd.getAtomList(true, true, true, false);
            int atomList2[] = resPrior.getAtomList(false, true, true, true);
            int[] atomList = concatenateArrays(atomList1, atomList2);

            m.setTorsion( m.residue[resDirectlyAffected[resDirectlyAffected.length-1]+1].getAtomNameToMolnum("N") ,
                    resEnd.getAtomNameToMolnum("C"), resEnd.getAtomNameToMolnum("CA"),
                    resEnd.getAtomNameToMolnum("N"), phiPsi[1], atomList, atomList.length, false);


            //Now phi
            atomList1 = resEnd.getAtomList(true, false, false, false);
            atomList2 = resPrior.getAtomList(false, true, true, true);
            atomList = concatenateArrays(atomList1, atomList2);

            m.setTorsion( resEnd.getAtomNameToMolnum("C"), resEnd.getAtomNameToMolnum("CA"),
                    resEnd.getAtomNameToMolnum("N"), resPrior.getAtomNameToMolnum("C"),
                    phiPsi[0], atomList, atomList.length, false);
        }
        else
            anchor2Old = m.getActualCoord( m.residue[resDirectlyAffected.length-1].getAtomNameToMolnum("CA") );


        Residue anchorRes1 = m.residue[resDirectlyAffected[closeStart]];
        Residue middleRes = m.residue[resDirectlyAffected[closeStart+1]];
        Residue anchorRes2 = m.residue[resDirectlyAffected[closeStart+2]];

        double anchor1New[] = m.getActualCoord( anchorRes1.getAtomNameToMolnum("CA") );
        double anchor2New[] = m.getActualCoord( anchorRes2.getAtomNameToMolnum("CA") );



        //Finally translate everything in between--the stuff to be closed--into place

        //The middle CA and everything before it is translated to avoid a discontinuity at the first CA of the tripeptide
        //Everything after is translated to avoid a discontinuity at the last CA of the tripeptide
        //But each translation only needs to happen if that CA has been messed with during secondary structure application


        if(addSSStart){
            int atomList[] = concatenateArrays( anchorRes1.getAtomList(false, false, false, true),
                middleRes.getAtomList(true, true, true, false) );

            double tr[] = rm.subtract(anchor1New, anchor1Old);
            m.translateAtomList(atomList, tr, false, false);
        }

        if(addSSEnd){
            int atomList[] = concatenateArrays( middleRes.getAtomList(false, false, false, true),
                anchorRes2.getAtomList(true, false, false, false) );

            double tr[] = rm.subtract(anchor2New, anchor2Old);
            m.translateAtomList(atomList, tr, false, false);
        }
        

        //Now the two rotations used to apply the loop closure adjustment will put the backbone in right
        //and sidechain idealization will fix the sidechains
    }
    
    
    private int[] concatenateArrays(int[]... toMerge){
        //Concatenates integer arrays (e.g. atom lists)

        int fullSize = 0;

        for ( int[] arr : toMerge )
            fullSize += arr.length;

        int[] ans = new int[fullSize];

        int cursor = 0;

        for( int[] arr : toMerge ){
            System.arraycopy(arr, 0, ans, cursor, arr.length);
            cursor += arr.length;
        }

        return ans;
    }


    
    @Override
    public boolean isParamAngle(){
        return false;
    }




}
