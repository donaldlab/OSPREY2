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
//	CETMatrix.java
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
import cern.colt.matrix.DoubleFactory2D;
import java.io.Serializable;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;


//Matrix of ContETerms for the various intra+shell and pairwise energies

public class CETMatrix implements Serializable {

    //these arrays are indexed by (residue #'d within active site, AA index, rotamer/RC #)
    //simlar to PairwiseEnergyMatrix

    ContETerm intraAndShellBounds[][][] = null;
    ContETerm pairwiseBounds[][][][][][] = null;
    
    
    int numRes;
    double ivalCutoff;//RCs unpruned at this ival are guaranteed to be included in this matrix
    //lets us know if we need to recompute the matrix 

    DegreeOfFreedom[] DOFList;
    
    String resAATypes[][];//Names of AA types for each residue



    public CETMatrix(int numMutable, int[][] strandMut, RotamerSearch rs){

        numRes = numMutable;
        DOFList = rs.m.DOFs;

        intraAndShellBounds = new ContETerm[numRes][][];
        pairwiseBounds = new ContETerm[numRes][][][][][];
        
        resAATypes = new String[numRes][];

        int p1 = 0;//First residue number within active site (not a loop variable because we loop over strands str1, then positions i within a strand)

        for (int str1=0; str1<strandMut.length; str1++){
            for (int i=0; i<strandMut[str1].length; i++){

                intraAndShellBounds[p1] = new ContETerm[rs.strandRot[str1].rl.getNumAAallowed()][];
                pairwiseBounds[p1] = new ContETerm[rs.strandRot[str1].rl.getNumAAallowed()][][][][];
                
                resAATypes[p1] = rs.strandRot[str1].rl.getAAtypesAllowed();

                for (int a1=0; a1<rs.strandRot[str1].getNumAllowable(strandMut[str1][i]); a1++){
                    int curAAind1 = rs.strandRot[str1].getIndexOfNthAllowable(strandMut[str1][i],a1);
                    int numRot1 = rs.getNumRot(str1, strandMut[str1][i], curAAind1);

                    intraAndShellBounds[p1][curAAind1] = new ContETerm[numRot1];
                    pairwiseBounds[p1][curAAind1] = new ContETerm[numRot1][][][];

                    for (int r1=0; r1<numRot1; r1++){
                        pairwiseBounds[p1][curAAind1][r1] = new ContETerm[numRes][][];

                        int p2=0;
                        for (int str2=0;str2<strandMut.length;str2++){
                                for (int j=0; j<strandMut[str2].length; j++){

                                        pairwiseBounds[p1][curAAind1][r1][p2] = new ContETerm[rs.strandRot[str2].rl.getNumAAallowed()][];
                                        for (int a2=0; a2<rs.strandRot[str2].getNumAllowable(strandMut[str2][j]); a2++){

                                            int curAAind2 = rs.strandRot[str2].getIndexOfNthAllowable(strandMut[str2][j],a2);
                                            int numRot2 = rs.getNumRot( str2, strandMut[str2][j], curAAind2 );
                                            pairwiseBounds[p1][curAAind1][r1][p2][curAAind2] = new ContETerm[numRot2];
                                        }
                                        p2++;
                                }
                        }
                    }
                }

                p1++;
            }
        }
    }




    

    public void mergeIn(CETMatrix M, int residueMutatable[]){
        //merge in bounds from a matrix M from a slave node,
        //which has the pairwise or intra+shell bounds for
        //the residue(s) indicated by 1s in residueMutatable


        int res1 = -1, res2 = -1;
        for(int i=numRes-1;i>=0;i--){//Go through the residues backwards to make sure res1>res2 (matching the matrix organization)
            //actually not needed now for LSBMatrix
                if (residueMutatable[i]==1) {
                        if (res1 == -1)
                                res1 = i;
                        else
                                res2 = i;
                }
        }


        if(res1 == -1)//Template run: nothing to merge in
            return;
        else if(res2 == -1) {
            intraAndShellBounds[res1] = M.intraAndShellBounds[res1];
        }
        else{

            for(int AAind=0; AAind<pairwiseBounds[res1].length; AAind++){
                if( pairwiseBounds[res1][AAind] != null ){
                    for(int rot=0; rot<pairwiseBounds[res1][AAind].length; rot++){
                        pairwiseBounds[res1][AAind][rot][res2] = M.pairwiseBounds[res1][AAind][rot][res2];
                    }

                }
            }

            for(int AAind=0; AAind<pairwiseBounds[res2].length; AAind++){
                if( pairwiseBounds[res2][AAind] != null ){
                    for(int rot=0; rot<pairwiseBounds[res2][AAind].length; rot++){
                        pairwiseBounds[res2][AAind][rot][res1] = M.pairwiseBounds[res2][AAind][rot][res1];
                    }

                }
            }
        }
    }


    public void setShellRotE(int res, int AAind, int rot, ContETerm b){
        tryMolecCompression(res,AAind,rot,b);
        intraAndShellBounds[res][AAind][rot] = b;
    }

    public void setPairwiseE(int res1, int AAind1, int rot1, int res2, int AAind2,
            int rot2, ContETerm b){
        
        tryMolecCompression(res1,AAind1,rot1,res2,AAind2,rot2,b);

        pairwiseBounds[res1][AAind1][rot1][res2][AAind2][rot2] = b;
        pairwiseBounds[res2][AAind2][rot2][res1][AAind1][rot1] = b;
    }
    
    
    //tryMolecCompression: these methods are used when adding a new series b w/ SVE to the mtarix
    //we try to substitute in a molecule that is already in the matrix
    //to avoid storing a new one in b
    //we can do this iff there is another bound with SVE at the same residue with the same 
    //discreteDOFVals (including amino acid type)
    //since then setting continuous DOFs (during ContETerm evaluation) will make the relevant residues
    //of the molecules identical
    void tryMolecCompression(int res, int AAind, int rot, ContETerm b){
        if(b instanceof EPoly){
            EPoly ep = (EPoly)b;
            
            if(ep.sve != null){
                
                for(ContETerm b2 : intraAndShellBounds[res][AAind]){//can only share m with other residues at the same (res, AA type)
                    if(b2 instanceof EPoly){//and thus not null...
                        EPoly ep2 = (EPoly)b2;
                        
                        if(ep2.sve!=null){
                            
                            boolean match = true;
                            //since res and AA type match, b and b2 will have the same discrete DOFs
                            //check that the values match
                            for(int dd=0; dd<b.discreteDOFVals.size(); dd++){
                                
                                double dd1 = b.discreteDOFVals.get(dd);//want to compare doubles, not Doubles in the array list...need numerical comparison
                                double dd2 = b2.discreteDOFVals.get(dd);
                                
                                if(dd1!=dd2){
                                    match = false;
                                    break;
                                }
                            }
                            
                            if(match){
                                ep.sve.m = ep2.sve.m;
                                ep.sveOF.m = ep2.sve.m;
                                return;//succesfully compressed
                            }
                        }
                    }
                }
            }
        }
    }


    //pairwise version
    void tryMolecCompression(int res1, int AAind1, int rot1, int res2, int AAind2,
            int rot2, ContETerm b){
        
        if(b instanceof EPoly){
            EPoly ep = (EPoly)b;
            
            if(ep.sve != null){
                
                for(ContETerm c[][][] : pairwiseBounds[res1][AAind1]){//can only share m with other residues at the same (res, AA type)
                    if(c!=null){
                        for(ContETerm b2 : c[res2][AAind2]){
                            
                            if(b2 instanceof EPoly){//and thus not null...
                                EPoly ep2 = (EPoly)b2;

                                if(ep2.sve!=null){
                                    boolean match = true;
                                    //since res and AA type match, b and b2 will have the same discrete DOFs
                                    //check that the values match
                                    for(int dd=0; dd<b.discreteDOFVals.size(); dd++){
                                        double dd1 = b.discreteDOFVals.get(dd);//want to compare doubles, not Doubles in the array list...need numerical comparison
                                        double dd2 = b2.discreteDOFVals.get(dd);

                                        if(dd1!=dd2){
                                            match = false;
                                            break;
                                        }
                                    }
                                    if(match){
                                        ep.sve.m = ep2.sve.m;
                                        ep.sveOF.m = ep2.sve.m;
                                        return;//succesfully compressed
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    public ContETerm[] getCETList(int[] AAIndices, int rot[]){
        //given the AAs and rotamers assigned so far, get the list of LSBs
        //we'll use to compute the A* lower bound

        int numTerms = (AAIndices.length+1)*AAIndices.length/2;

        ContETerm terms[] = new ContETerm[numTerms];

        int termCount = 0;

        for(int res=0; res<AAIndices.length; res++){
            //Note: this works if all residues are assigned
            //or just the first several

            terms[termCount] = intraAndShellBounds[res][AAIndices[res]][rot[res]];
            termCount++;

            for(int res2=0; res2<res; res2++){
                terms[termCount] = pairwiseBounds[res][AAIndices[res]][rot[res]][res2][AAIndices[res2]][rot[res2]];
                termCount++;
            }
        }


        return terms;
    }
    
    
    
    
    
    
    
    public CETObjFunction getObjFunc( int[] AAIndices, int rot[], boolean includeMinE, boolean rcs, ContSCObjFunction sveOF){
        //Get the objective function for the given rotamers
        //if !rcs: given as (AA#, rot#) assignments for each residue (can be -1)
        //if rcs: given as RC set assignments, in AStarAxe.RCSets numbering, stored in rot (AAIndices can be null)
        //if sveOF isn't null we can use it to assign DOFs for any SVEs this objective function may call
        
        ContETerm terms[] = getCETList(AAIndices, rot);

        //we'll minimize over all DOFs that the terms depend on
        //constraints will be taken from DOFmin, DOFmax in the bounds
        //when getting constraints all bounds should agree (either using
        //single-voxel or all-unpruned voxels range of values for each DOF...)
        HashMap<Integer,double[]> DOFMap = new HashMap<Integer,double[]>();//map DOFNum-->constraints
        for(ContETerm b : terms){

            if(b!=null){

                for(int dof=0; dof<b.numDOFs; dof++){
                    double constr[] = new double[] {b.DOFmin.get(dof),b.DOFmax.get(dof)};
                    int DOFNum = b.DOFNums[dof];

                    //check against any previous bounds' constraints on same DOF...
                    if(DOFMap.containsKey(DOFNum)){
                        double prevConstr[] = DOFMap.get(DOFNum);
                        prevConstr[0] = Math.max(constr[0],prevConstr[0]);
                        prevConstr[1] = Math.min(constr[1],prevConstr[1]);
                        //Each assigned RC restricts the relevant DOFs to
                        //the range for that RC, and removes RCs at other positions
                        //in the case of pruned pairs.  So all the DOF upper and lower bounds
                        //in terms apply to the conformations set we're minimizing over--
                        //thus we use the intersection of the ranges
                        //that different terms may provide for our DOF
                        
                        if(prevConstr[0] > prevConstr[1]){
                            System.err.println("ERROR: inconsistent DOF constraints for objective function! DOFNum: "
                                    +DOFNum+" prevConstr: "+prevConstr[0]+" "+prevConstr[1]+" new constr: "+
                                    constr[0]+" "+constr[1]);
                        }
                    }
                    else
                        DOFMap.put(DOFNum, constr);
                }
            }
        }

        int numDOFs = DOFMap.size();

        int DOFNums[] = new int[numDOFs];
        double constraints[][] = new double[2][numDOFs];
        int count = 0;

        for(int DOFNum : DOFMap.keySet()){
            DOFNums[count] = DOFNum;
            double constr[] = DOFMap.get(DOFNum);
            constraints[0][count] = constr[0];
            constraints[1][count] = constr[1];
            count++;
        }

        DoubleMatrix1D DOFVals = DoubleFactory1D.dense.make(numDOFs);
        CETEnergyFunction ef = new CETEnergyFunction( DOFVals, DOFNums, terms, includeMinE );
        return new CETObjFunction(DOFNums,constraints,DOFList,ef,sveOF);
    }

    
    
    public void cleanupSVE(){
        //Whenever two ContETerms have the same amino-acid types,
        //AND any discreteDOFVals agree,
        //we can use the same molecule for sve (since the relevant res of the molecule will be exactly the same
        //once the continuous DOFs are set)
        //so let's link all the redundant molecules
        for(ContETerm[][] b : intraAndShellBounds){
            if(b!=null)
            for(ContETerm[] bb : b){
                if(bb!=null){
                    
                    //Molecule mutMolec = null;//molecule with given mutations
                    //should be the same for all bbb in bb
                    //WAIT have to account for discreteDOFVals
                    
                    ArrayList<Molecule> mutMolec = new ArrayList<>();
                    ArrayList<double[]> mutMolecDD = new ArrayList<>();//discrete DOFVals for each mutMolec
                    
                    for(ContETerm bbb : bb){
                        if(bbb instanceof EPoly){//not null and may have sve
                            EPoly ep = (EPoly)bbb;
                            if(ep.sve!=null){
                                
                                /*if(mutMolec==null)
                                    mutMolec = ep.sve.m;
                                else{
                                    ep.sve.m = mutMolec;
                                    ep.sveOF.m = mutMolec;
                                }*/
                                boolean alreadyStored = false;
                                for(int m=0; m<mutMolec.size(); m++){
                                    double[] ddVal = mutMolecDD.get(m);
                                    boolean match = true;
                                    for(int dd=0; dd<ddVal.length; dd++)//discreteDOFVals must match
                                        //if AA types, res match then there'll be the same number of discrete DOFs
                                        match = match && (ep.discreteDOFVals.get(dd)==ddVal[dd]);

                                    if(match){
                                        ep.sve.m = mutMolec.get(m);
                                        ep.sveOF.m = ep.sve.m;
                                        alreadyStored = true;
                                        break;
                                    }
                                }

                                if(!alreadyStored){
                                    mutMolec.add(ep.sve.m);
                                    double ddVals[] = new double[ep.discreteDOFVals.size()];
                                    for(int dd=0; dd<ddVals.length; dd++)
                                        ddVals[dd] = ep.discreteDOFVals.get(dd);
                                    mutMolecDD.add(ddVals);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        for(ContETerm[][][][][] b : pairwiseBounds){
            if(b!=null)
            for(ContETerm[][][][] bb : b){
                if(bb!=null){
                    
                    //to reuse molecule, first & second res & AA must match
                    //HashMap<Integer,Molecule> mutMolec = new HashMap<>();
                    //key for mutMolec is res2 + numRes*curAA2
                    //which is OK because res2 can't exceed numRes
                    
                    //We need a mutMolec for each (res2, curAA2, any discreteDOFVals)
                    ArrayList<Molecule> mutMolec = new ArrayList<>();
                    ArrayList<double[]> mutMolecInfo = new ArrayList<>();
                    //for each molecule in mutMolec, (res2, curAA2, any discreteDOFVals)
                    
                    for(ContETerm bbb[][][] : bb){
                        if(bbb!=null)
                        for(int res2=0; res2<bbb.length; res2++){
                            ContETerm bbbb[][] = bbb[res2];
                            if(bbbb!=null)
                            for(int curAA2=0; curAA2<bbbb.length; curAA2++){
                                ContETerm bbbbb[] = bbbb[curAA2];
                                if(bbbbb!=null)
                                for(ContETerm bbbbbb: bbbbb){
                                    if(bbbbbb instanceof EPoly){//not null and may have sve
                                        EPoly ep = (EPoly)bbbbbb;
                                        if(ep.sve!=null){
                                            
                                            /*if(mutMolec.get(res2 + numRes*curAA2)==null)
                                                mutMolec.put(res2 + numRes*curAA2, ep.sve.m);
                                            else{
                                                ep.sve.m = mutMolec.get(res2 + numRes*curAA2);
                                                ep.sveOF.m = ep.sve.m;
                                            }*/
                                            boolean alreadyStored = false;
                                            for(int m=0; m<mutMolec.size(); m++){
                                                double[] info = mutMolecInfo.get(m);
                                                boolean match = (info[0]==res2) && (info[1]==curAA2);
                                                for(int dd=2; dd<info.length; dd++)//discreteDOFVals must match
                                                    //if AA types, res match then there'll be the same number of discrete DOFs
                                                    match = match && (ep.discreteDOFVals.get(dd-2)==info[dd]);
                                                
                                                if(match){
                                                    ep.sve.m = mutMolec.get(m);
                                                    ep.sveOF.m = ep.sve.m;
                                                    alreadyStored = true;
                                                    break;
                                                }
                                            }
                                            
                                            if(!alreadyStored){
                                                mutMolec.add(ep.sve.m);
                                                double info[] = new double[ep.discreteDOFVals.size()+2];
                                                info[0] = res2;
                                                info[1] = curAA2;
                                                for(int dd=2; dd<info.length; dd++)
                                                    info[dd] = ep.discreteDOFVals.get(dd-2);
                                                mutMolecInfo.add(info);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    

    
    static double squareDist(DoubleMatrix1D a, DoubleMatrix1D b){
        //square of distance between two vectors
        return a.zDotProduct(a) - 2*a.zDotProduct(b) + b.zDotProduct(b);
    }
    
    
    
    void analyzeFitTypes(){
        //Look at what types of fits are in this matrix
        System.out.println("Fit types in continuous energy term matrix (# of terms, description): ");
        
        HashMap<String,Integer> typeCounts = new HashMap<>();
        
        for(ContETerm c[][] : intraAndShellBounds){
            for(ContETerm cc[] : c){
                if(cc!=null){
                    for(ContETerm ccc : cc){
                        if(ccc!=null)
                            countFitType(ccc.fitDescription,typeCounts);
                    }
                }
            }
        }
        
        for(int res=0; res<numRes; res++){
            for(ContETerm cc[][][][] : pairwiseBounds[res]){
                if(cc!=null){
                    for(ContETerm ccc[][][] : cc){
                        if(ccc!=null){
                            for(int res2=0; res2<res; res2++){
                                for(ContETerm cccc[] : ccc[res2]){
                                    if(cccc!=null){
                                        for(ContETerm ccccc : cccc){
                                            if(ccccc!=null)
                                                countFitType(ccccc.fitDescription,typeCounts);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        for(String descr : typeCounts.keySet()){
            System.out.println(typeCounts.get(descr)+" "+descr);
        }
    }
    
    
    void countFitType(String descr, HashMap<String,Integer> typeCounts){
        //given the description of a fit, add it to the counts of descriptions in typeCounts
        if(typeCounts.containsKey(descr))
            typeCounts.put(descr, typeCounts.get(descr)+1);
        else
            typeCounts.put(descr, 1);
    }


    
}