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
//	DegreeOfFreedom.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.util.BitSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.io.Serializable;
/**
 *
 * @author mhall44
 */
public class DegreeOfFreedom implements RyanComparable, Serializable {//DOF for short

    //Possible types
    //Listed in the usual order for A*
    static int AATYPE = 0;
    static int STRRIGIDMOTION = 1;//Strand translation or rotation
    //This type could probably be treated as a perturbation
    static int PERTURBATION = 2;
    static int SCDIHEDRAL = 3;//Side-chain dihedral
    
    
    //new type?
    static int MUTATION = 4;

    

    int type;
    int DOFNum;//The number of this DOF, numbered among all DOFs in the system
    //boolean continuous;//Continuous DOFs allow minimization & can be used in A* node potential terms


    //Type-specific information
    int strResNum;//Strand residue number of affected residue: for AA types and dihedrals
    int dihedNum;//For dihedrals
    Perturbation pert;//For perturbations
    int pertNum;//index of pert in m.perts
    int rigidDOFNum;//Which rigid-body motion this is (must be in range 0-5)
    int strandNum;//Affected strand: for everything but perturbations 

    
    int flexResAffected[];//Residues affected by the DOF
    //indicated by numbers among flexible residues (e.g. those handled by A*)
    //Note: this includes all effects, whether direct or indirect (the distinction is for perturbations)
    //and also "flexible residues" means flexible in the design, not just in any given pairwise minimization


    //Interval info.  For now, intervals corresponding to RCs are used
    //Fused intervals can be attempted later (to reduce tree branching factor)

    int numIntervals;

    double intervals[][] = null;//Intervals indicated as bounds (the default way)
    //AA types are inconvenient to represent this way though
    
    //The following three fields are only for AA-type DOFs
    String WTAAType = null;//The wildtype AA type
    int WTAANum = -1;//its corresponding number
    String AATypeSets[] = null;//sets of amino-acid types take the place of intervals
    boolean WTCompatible[] = null;//Indicates which intervals are compatible with WT


    static boolean makeMutDOFs = true;//signals us to create mutation DOFs
    //associated with mutation to a particular AA type (mutAANum)
    //redundant with AATYPE, but suitable for making ellipse bounds
    int mutAANum;
    
    //Compatibility with RCs

    BitSet compatibleRCs[][][];//Indicates which RCs at each affected residue are compatible with each interval
    //First index: Which interval this is
    //Second index: Which of the affected residues this is
    //Third index: Which AA type this is
    //BitSet index: Which of the RCs at that residue && AA this is


    //This field is not created by constructors but by reduceCompatibleRCs
    BitSet reducedCompatibleRCs[][];//Indicates which unpruned RCs at each affected residue are compatible with each interval
    //First index: Which interval this is
    //Second index: Which of the affected residues this is
    //BitSet index: Which of the unpruned RCs at that residue this is (i.e. reduced numbering from A*)


    //Constructor for AA-type, mutation, or dihedral DOFs
    //if mutDOF==true we give mutAANum instead of dihedNumber
    public DegreeOfFreedom(int DOFNumber, int strNumber, int strResNumber, int flexResNumber, 
            String WT, int dihedNumber, StrandRotamers strandRot, boolean mutDOF){

        DOFNum = DOFNumber;
        strandNum = strNumber;
        strResNum = strResNumber;
        flexResAffected = new int[] { flexResNumber };

        if(dihedNumber == -1)
            initAATypeDOF(flexResNumber, WT, strandRot);
        else if(mutDOF)
            initMutDOF(flexResNumber, dihedNumber, strandRot);
        else
            initSCDihedral(flexResNumber, dihedNumber, strandRot);
    }
    

    public final void initAATypeDOF(int flexResNumber, String WT, StrandRotamers strandRot ){
        
            type = AATYPE;
            //continuous = false;

            WTAAType = WT;
            WTAANum = strandRot.rl.getAARotamerIndex(WT);

            numIntervals = strandRot.getNumAllowable(strResNum);
            WTCompatible = new boolean[numIntervals];
            AATypeSets = new String[numIntervals];

            compatibleRCs = initializeBitSets(numIntervals, 1, strandRot.rl.getNumAAallowed());
            
            for(int AA=0; AA<numIntervals; AA++){
                int AAIndex = strandRot.getIndexOfNthAllowable(strResNum, AA);
                AATypeSets[AA] = strandRot.rl.getAAName(AAIndex);
                if( AATypeSets[AA].equalsIgnoreCase(WTAAType) )
                    WTCompatible[AA] = true;

                int numRCs = getNumRCs(strandRot, AAIndex);
                compatibleRCs[AA][0][AAIndex].set(0,numRCs);//Set all RCs with the current AA type to be compatible with it
            }
    }
    
    
    public final void initMutDOF(int flexResNumber, int mutAAIndex, StrandRotamers strandRot ){
        
            type = MUTATION;
            mutAANum = mutAAIndex;
            numIntervals = 2;//we can either have this mutation (param==1) or not (param==0)
            compatibleRCs = initializeBitSets(numIntervals, 1, strandRot.rl.getNumAAallowed());
            
            intervals = new double[][] { {0,0}, {1,1} };
            
            int numAllowable = strandRot.getNumAllowable(strResNum);
            
            for(int AA=0; AA<numAllowable; AA++){
                int AAIndex = strandRot.getIndexOfNthAllowable(strResNum, AA);
                int numRCs = getNumRCs(strandRot, AAIndex);
                
                if(AAIndex==mutAAIndex)//param==1
                    compatibleRCs[1][0][AAIndex].set(0,numRCs);//Set all RCs with the current AA type to be compatible with it
                else
                    compatibleRCs[0][0][AAIndex].set(0,numRCs);
            }
    }


    public final void initSCDihedral(int flexResNumber,
            int dihedNumber, StrandRotamers strandRot){
        
        type = SCDIHEDRAL;
        //continuous = true;
        dihedNum = dihedNumber;

        int numAllowable = strandRot.getNumAllowable(strResNum);

        //We collect all the ideal rotameric values for this dihedral
        //for all possible AA types in this HashSet.  They will be used to make intervals
        HashSet<Integer> idealDihed = new HashSet<Integer>();

        for(int AA=0; AA<numAllowable; AA++){
            int AAIndex = strandRot.getIndexOfNthAllowable(strResNum, AA);
            if( dihedNum < strandRot.rl.getNumDihedrals(AAIndex) ){//If this dihedral is present for this residue
                int numRot = strandRot.rl.getNumRotForAAtype( AAIndex );
                for(int rot=0; rot<numRot; rot++)//Actual sidechain rotamers, not RCs
                    idealDihed.add( strandRot.rl.getRotamerValues(AAIndex, rot, dihedNum) );
            }
        }

        //Now make intervals for each of the ideal dihedral values
        //Plus an extra null interval, corresponding to amino-acid types for which
        //this dihedral is not defined

        double maxMovement = (double)ContSCObjFunction.maxMovement;//Maximum continuous sidechain movement from ideal value

        numIntervals = idealDihed.size() + 1;
        intervals = new double[numIntervals][];
        compatibleRCs = initializeBitSets(numIntervals, 1, strandRot.rl.getNumAAallowed());

        int interv=0;

        for( int dval : idealDihed ){
            intervals[interv] = new double[] { dval-maxMovement, dval+maxMovement };

            //Fill in compatibility
            for(int AA=0; AA<numAllowable; AA++){
                int AAIndex = strandRot.getIndexOfNthAllowable(strResNum, AA);

                if( dihedNum < strandRot.rl.getNumDihedrals(AAIndex) ){//If this dihedral is present for this residue
                    int numRot = strandRot.rl.getNumRotForAAtype( AAIndex );

                    for(int rot=0; rot<numRot; rot++){//Actual sidechain rotamers, not RCs
                        int rotDihed = strandRot.rl.getRotamerValues(AAIndex, rot, dihedNum);

                        if(rotDihed == dval){//Rotamer rot corresponds to the current intervals

                            if( strandRot instanceof StrandRCs ){//Set all RCs containing the current rotamer to be compatible with the interval
                                int[] curRCRots = ((StrandRCs)strandRot).RCRots[strResNum][AAIndex];
                                for(int RC=0; RC<curRCRots.length; RC++){
                                    if(curRCRots[RC] == rot)
                                        compatibleRCs[interv][0][AAIndex].set(RC);
                                }
                            }
                            else//Just set the current rotamer to be compatible with the interval
                                compatibleRCs[interv][0][AAIndex].set(rot);
                        }
                    }
                }
            }

            interv++;
        }

        //Leave interval[numIntervals-1] null but fill in compatibility
        for(int AA=0; AA<numAllowable; AA++){
                int AAIndex = strandRot.getIndexOfNthAllowable(strResNum, AA);

                if( dihedNum >= strandRot.rl.getNumDihedrals(AAIndex) ){//If this dihedral is absent for this residue
                    int numRCs = getNumRCs(strandRot, AAIndex);
                    compatibleRCs[interv][0][AAIndex].set(0,numRCs);
                }
        }

    }


    //Constructor for perturbations
    public DegreeOfFreedom(int DOFNumber, int pertNumber, Perturbation pe, int flexToMolResNumMap[], StrandRotamers strandRot[]){

        type=PERTURBATION;
        DOFNum = DOFNumber;
        //continuous = pert.
        pert = pe;
        pertNum = pertNumber;


        Molecule molec = pert.m;

        //Make a list of affected residues, and a corresponding list of flexible residue numbers
        ArrayList<Integer> affectedFlexRes = new ArrayList<Integer>();
        ArrayList<Residue> affectedRes = new ArrayList<Residue>();

        for(int flexRes=0; flexRes<flexToMolResNumMap.length; flexRes++){
            Residue res = molec.residue[flexToMolResNumMap[flexRes]];
            for(int p=0; p<res.perts.length; p++){
                if(res.perts[p]==pertNum){
                    affectedRes.add(res);
                    affectedFlexRes.add(flexRes);
                    break;
                }
            }
        }

        int numResAffected = affectedRes.size();
        flexResAffected = new int[numResAffected];
        for(int a=0; a<numResAffected; a++)
            flexResAffected[a] = affectedFlexRes.get(a);
        

        numIntervals = pert.minParams.length;
        intervals = new double[numIntervals][2];

        compatibleRCs = new BitSet[numIntervals][numResAffected][];

        for(int a=0; a<numIntervals; a++){
            
            intervals[a][0] = pert.minParams[a];
            intervals[a][1] = pert.maxParams[a];

            //Fill in RC compatibility
            for(int ares=0; ares<numResAffected; ares++){//Loop through affected residues
                Residue res = affectedRes.get(ares);
                int strandNumber = res.strandNumber;
                int strandResNumber = res.strandResidueNumber;

                StrandRCs sRC = (StrandRCs)strandRot[strandNumber];

                int numAATypes = sRC.getNumAATypes();

                compatibleRCs[a][ares] = initializeBitSets(numAATypes);

                for(int AAIndex=0; AAIndex<numAATypes; AAIndex++ ){
                    for(int RC=0; RC<sRC.getNumRCs(strandResNumber, AAIndex); RC++){
                        int resPertState = sRC.RCPertStates[strandResNumber][AAIndex][RC];

                        //Figure out which of the residues affecting this perturbations is pert
                        for( int resPertNum=0; resPertNum<res.perts.length; resPertNum++ ){
                            if( res.perts[resPertNum] == pertNum ){
                                //resPertNum is the number of perts among perturbations affecting res
                                int curInterval = res.pertStates[resPertState][resPertNum];

                                if( curInterval == a )
                                    compatibleRCs[a][ares][AAIndex].set(RC);
                                
                                break;
                            }
                        }
                    }
                }
            }
        }
    }


    
    //Constructor for strand translations/rotations
    public DegreeOfFreedom(int DOFNumber, int rigidDOFNumber, int strandNumber, int mutRes2Strand[], StrandRotamers strandRot){

        DOFNum = DOFNumber;
        strandNum = strandNumber;
        type = STRRIGIDMOTION;
        rigidDOFNum = rigidDOFNumber;


        //Figure out which residues are affected by this (it will be all the residues on this strand)
        int numResAffected = 0;
        for(int str : mutRes2Strand){
            if(str == strandNum)
                numResAffected++;
        }

        flexResAffected = new int[numResAffected];

        int affRes=0;

        for(int flexRes=0; flexRes<mutRes2Strand.length; flexRes++){
            if( mutRes2Strand[flexRes] == strandNum ){
                flexResAffected[affRes] = flexRes;
                affRes++;
            }
        }

        //For now, we will only have one interval for these motions
        //Translation up to 1.2 angstroms in either direction, as in SimpleMinimizer
        //Rotation up to 5 degrees in either direction (not specified in SimpleMinimizer,
        //but this is consistent with prior behavior and probably about right for a local minimizer)
        numIntervals = 1;

        if(rigidDOFNumber<3)//translation
            intervals = new double[][] { new double[] { -1.2f, 1.2f } };
        else//rotation
            intervals = new double[][] { new double[] { -5f, 5f } };

        //All RCs at residues in the strand are compatible with this interval
        compatibleRCs = initializeBitSets(numIntervals, numResAffected, strandRot.rl.getNumAAallowed());

        for(int ares=0; ares<numResAffected; ares++){//Loop through affected residues

            int numAllowable = strandRot.getNumAllowable(strResNum);

            for(int AA=0; AA<numAllowable; AA++){
                int AAIndex = strandRot.getIndexOfNthAllowable(strResNum, AA);
                int numRCs = getNumRCs(strandRot, AAIndex);
                compatibleRCs[0][ares][AAIndex].set(0,numRCs);
            }
        }
    }


    
    //Make an array of all the degrees of freedom in a system
    //The intervals will be sufficient to cover all the RCs in strandRot
    public static DegreeOfFreedom[] makeDOFArray( StrandRotamers strandRot[], int strandMut[][], int mutRes2Strand[],
            int mutRes2StrandMutIndex[], Molecule m, String strandDefault[][] ){

        int numFlexRes = mutRes2Strand.length;

        //ORDER:
        //AA types and dihedrals by residue
        //then all perturbations
        //then trans/rot by strand
        

        ArrayList<DegreeOfFreedom> DOFList = new ArrayList<DegreeOfFreedom>();
        int DOFNum=0;

        for(int flexResNum=0; flexResNum<numFlexRes; flexResNum++){//Loop through residues
            //to make AA type and dihedral DOFs

            int str = mutRes2Strand[flexResNum];
            int strandResNum = strandMut[str][mutRes2StrandMutIndex[flexResNum]];
            
            String WT = m.strand[str].residue[strandResNum].name;//by default WT is current AA type
            if(strandDefault!=null)
                WT = strandDefault[str][mutRes2StrandMutIndex[flexResNum]];//wild-type amino acid
            //Note: as per KSParser.checkWT, if there is only one AA type allowed at a position it is treated as WT
          

            //Figure out how many dihedral the residue has
            int numDihedrals = 0;
            int numAllowable = strandRot[str].getNumAllowable(strandResNum);
            for(int AA=0; AA<numAllowable; AA++){
                int AAIndex = strandRot[str].getIndexOfNthAllowable(strandResNum, AA);
                numDihedrals = Math.max( numDihedrals, strandRot[str].rl.getNumDihedrals(AAIndex) );
                
                if(makeMutDOFs){
                    //make mutation DOFs
                    if( !WT.equalsIgnoreCase(strandRot[str].rl.getAAName(AAIndex) ) ){//not wild-type
                        DOFList.add( new DegreeOfFreedom(DOFNum,str,strandResNum,flexResNum,WT,AAIndex,strandRot[str],true) );
                        DOFNum++;
                    }
                }
            
            }

            for(int d=-1; d<numDihedrals; d++){//Start at -1 for AA type DOF, then make dihedral DOFs
                DOFList.add( new DegreeOfFreedom(DOFNum,str,strandResNum,flexResNum,WT,d,strandRot[str],false) );
                DOFNum++;
            }
        }


        if( m.perts.length > 0 ){//If there are perturbations
            //we need to create a map from molecule residue numbers to residue numbers among flexible residues
            int flexToMolResNumMap[] = new int[m.residue.length];
            for(int flexResNum=0; flexResNum<numFlexRes; flexResNum++){
                int str = mutRes2Strand[flexResNum];
                int strandResNum = strandMut[str][mutRes2StrandMutIndex[flexResNum]];
                int molResNum = m.strand[str].residue[strandResNum].moleculeResidueNumber;
                flexToMolResNumMap[flexResNum] = molResNum;
            }

            for(int p=0; p<m.perts.length; p++){
            //Inherit the intervals from the perturbations
                DOFList.add( new DegreeOfFreedom(DOFNum, p, m.perts[p], flexToMolResNumMap, strandRot) );
                DOFNum++;
            }
        }

        

        //Strand translations & rotations
        for( int str=0; str<m.strand.length; str++ ){
            if(m.strand[str].rotTrans){
                for(int rigidDOFNum=0; rigidDOFNum<6; rigidDOFNum++){
                    DOFList.add( new DegreeOfFreedom(DOFNum, rigidDOFNum, str, mutRes2Strand, strandRot[str]) );
                    DOFNum ++;
                }
            }
        }

        DegreeOfFreedom[] ans = new DegreeOfFreedom[DOFList.size()];
        ans = DOFList.toArray(ans);

        return ans;
    }


    
    public void reduceCompatibleRCs(int numFlexRes, ReducedEnergyMatrix rem){
        //Reduce the compatible RCs
        //using the RC indexing from rem

        reducedCompatibleRCs = initializeBitSets( numIntervals, flexResAffected.length );

        //Build a map from flexible residues to affected residues for this DOF
        int flexToAffectedResMap[] = new int[numFlexRes];
        Arrays.fill(flexToAffectedResMap, -1);
        for(int a=0; a<flexResAffected.length; a++)
            flexToAffectedResMap[flexResAffected[a]] = a;


        for(int i=0; i<numIntervals; i++ ){

            int numUnprunedRCs = rem.indicesEMatrixPos.length;


            //RC indices: redRCLong is among all RCs in rem (all unpruned RCs), as in rem
            //redRCShort is among all unpruned RCs at the given position
            //unredRC is among all RCs at the given position and AA type, as in strandRCs or strandRotamers
            int redRCShort = 0;
            int oldPos = -1;

            for(int redRCLong=0; redRCLong<numUnprunedRCs; redRCLong++){//Go through unpruned RCs

                int pos = rem.indicesEMatrixPos[redRCLong];
                if(pos!=oldPos){
                    redRCShort = 0;
                    oldPos = pos;
                }

                int affPos = flexToAffectedResMap[pos];

                if(affPos != -1){//If this RC is at a position affected by this DOF

                    int curAA = rem.indicesEMatrixAA[redRCLong];
                    int unredRC = rem.indicesEMatrixRot[redRCLong];

                    boolean compatible = compatibleRCs[i][affPos][curAA].get(unredRC);
                    reducedCompatibleRCs[i][affPos].set(redRCShort,compatible);
                }

                redRCShort++;
            }
        }
    }



    public double getMeshWidth(){//Get the mesh width to use for this DOF in a FuncLBTerm that depends on it
        
        if(type==SCDIHEDRAL)
            return 1.;
        else if(type==PERTURBATION)
            return pert.getMeshWidth();
        else if(type==STRRIGIDMOTION){
            if(rigidDOFNum<3)//translation
                return 0.1;
            else//rotation
                return 0.25;
        }
        else if(type==MUTATION)
            return 0.2;
        else{
            System.out.println("Warning: No mesh width defined for DOF type "+type);
            return Double.NaN;//not defined
        }
    }




    public int compareTo(Object other){
        
        DegreeOfFreedom DOF2 = (DegreeOfFreedom)other;
        //Returns 1 is this should come after other in the A* ordering,
        //-1 if before, 0 if it doesn't matter


        //For now let's use this simple ordering
        //Can refine later by using gap in energies for intervals,
        //or size of quadratic lower-bound terms, to rank
        int typeCompare = new Integer(type).compareTo(DOF2.type);

        if(typeCompare != 0)//Different types
            return typeCompare;

        //Look for useful type-specific orderings
        if(type == PERTURBATION)
            return new Integer(DOFNum).compareTo(DOF2.DOFNum);//Use order of application

        if(type == SCDIHEDRAL)
            return new Integer(dihedNum).compareTo(DOF2.dihedNum);

        //We haven't found a reason to specify ordering so we return 0
        return 0;
    }


    public boolean isContinuous(){
        //Check if the DOF has at least one interval of nonzero width,
        //allowing continuous minimization
        if( intervals == null )
            return false;

        for( double[] interval : intervals ){
            if(interval != null){
                if( interval[1] != interval[0] )
                    return true;
            }
        }
        return false;
    }


    public boolean isContinuousForRC(int curRes, int curAA, int curRC){
        //Check if the DOF can be continuously minimized in the given RC,
        //i.e. there is an interval of nonzero width compatible with the RC
        if( intervals == null )
            return false;

        //Convert curRes (numbered among all flexible residues) to numbering among
        //affected residue
        int ares = -1;
        for(int a=0; a<flexResAffected.length; a++){
            if(flexResAffected[a] == curRes){
                ares = a;
                break;
            }
        }

        if(ares==-1)//This DOF doesn't affect curRes so the effect can't be continuous
            return false;


        for( int intNum=0; intNum<intervals.length; intNum++ ){
            double[] interval = intervals[intNum];
            if(interval != null){
                if( interval[1] != interval[0] ){
                    if(compatibleRCs[intNum][ares][curAA].get(curRC))
                        return true;
                }
            }
        }
        return false;
    }



    public boolean hasOverlapRes(DegreeOfFreedom... otherDOFs){
        //Is there a residue affected by this and all the otherDOFs?

        for(int res1 : flexResAffected){

            boolean allAffect = true;//Do they all affect res1?

            for(DegreeOfFreedom dof : otherDOFs){
                boolean affected = false;//Does dof affect res1?
                for(int res2 : dof.flexResAffected){
                    if(res1==res2)
                        affected = true;
                }
                allAffect = allAffect && affected;
            }

            if(allAffect)
                return true;
        }

        return false;
    }


    //Initializing BitSet arrays with BitSet objects instead of null objects
    public static BitSet[] initializeBitSets(int num){
        BitSet ans[] = new BitSet[num];
        for(int a=0; a<num; a++)
            ans[a] = new BitSet();
        return ans;
    }

    public static BitSet[][] initializeBitSets(int num1, int num2){
        BitSet ans[][] = new BitSet[num1][];
        for(int a=0; a<num1; a++)
            ans[a] = initializeBitSets(num2);

        return ans;
    }
    
    public static BitSet[][][] initializeBitSets(int num1, int num2, int num3){
        BitSet ans[][][] = new BitSet[num1][][];
        for(int a=0; a<num1; a++)
            ans[a] = initializeBitSets(num2,num3);

        return ans;
    }



    //For a single-residue DOF, get the number of RCs at the current residue
    //for the given AAIndex
    private int getNumRCs(StrandRotamers strandRot, int AAIndex){

        if( strandRot instanceof StrandRCs )
            return ((StrandRCs)strandRot).getNumRCs(strResNum, AAIndex);
        else{
            int numRot = strandRot.rl.getNumRotForAAtype(AAIndex);
            if(numRot == 0)
                return 1;
            else
                return numRot;
        }
    }
    
    
    public boolean isDOFAngle(){
        //is this DOF an angle?
        if(type==SCDIHEDRAL)
            return true;
        else if(type==STRRIGIDMOTION){
            if(rigidDOFNum>2)//rotation angle
                return true;
        }
        else if(type==PERTURBATION){
            if(pert.isParamAngle())
                return true;
        }
        
        return false;
    }


       /* public int getNumIntervals(){
        if (type == AATYPE)
            return AATypeSets.length;
        else
            return intervals.length;
    }*/


        //Dihedrals for a residue have a "hard dependency" on AA type
    //because they are likely not well defined unless the AA type is defined
    //So for a dihedral hardDependency is the DOFNum of the AA type for that residue
    //int hardDependency = -1;
    //"Soft dependencies" are DOFs (indicated by their DOFNums) that should probably come before
    //the current DOF in the A* tree because they move the same atoms by a larger amount
    //(e.g. earlier dihedrals at a residue)
    //This list will contain both hard and soft dependencies
    //ArrayList<Integer> dependencies = new ArrayList<Integer>();
    //This can be handled in compareTo instead:
            //For dependency-based ordering: We will consider the ordering irrelevant if there is no overlap in affected residues
       /* if( hasOverlap( flexResAffected, DOF2.flexResAffected ) ){

            //The DOF affecting more residues comes first
            if( Math.max(flexResAffected.length, DOF2.flexResAffected.length) > 1 )
                return new Integer(DOF2.flexResAffected.length).compareTo(flexResAffected.length);

            //If they both just affect 1 residue, the AA type comes first
            if( type == AATYPE )
                return -1;
            if(DOF2.type == AATYPE)
                return 1;

            //If they're both SC dihedrals the one with the lower
            if( type == SCDIHEDRAL && DOF2.type == SCDIHEDRAL  )


            }

        }*/


}
