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
//	ContSCObjFunction.java
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
import java.util.HashSet;
import java.io.Serializable;

/**
 *
 * @author mhall44
 */
public class ContSCObjFunction implements ObjectiveFunction, Serializable {
    //This is an objective function for minimizing with respect to sidechain dihedrals
    //can have perturbation parameters and strand rigid-body motions as DOFs too



        //A lot of this is copied from SimpleMinimizer

        Molecule m = null;//Assumed to contain a DOF array (m.DOFs)
	EnergyFunction efunc = null;

        DihedralEnergy de = null;//If we're using dihedral energies, this
        //dihedral term will be allocated and added to efunc

        DoubleMatrix1D curDOFVals;//Value of degrees of freedom we're minimizing over.
        //Some components of the energy function may refer to this
        //(a DOFTermsEnergyFunction will do this)

        int DOFNums[];//DOF numbers for each degree of freedom being minimized over
        //(i.e. map from indexes in curDOFVals to indices in m.DOFs)
        int numDOFs;//Number of degrees of freedom

        StrandRotamers[] strandRot = null;
        int numberOfStrands = 0;

	int numFlexRes[] = null;
		// the number of flexible residues
        int totalFlexRes = 0;

        int flexToMolResMap[];//Map from molecule residue number to index among flexible residues


        double specialBoxConstraints[][];//Special constraints on individual degrees of freedom
        //If specialBoxConstraints[a] is null will use default values (e.g. +/- 9 degrees around ideal dihedrals)
        //Otherwise the bounds on the DOF will be specialBoxConstraints[a][0] and specialBoxConstraints[a][1]



        //DIHEDRAL STUFF
	int strDihedralAtNums[][][] = null;
		// array of dihedrals, this array is n by 4
		//  and has the moleculeAtomNumbers for each
		//  of the 4 atoms of the dihedral
	int strDihedralDistal[][][] = null;
		// array of dihedrals, the first index is over all
		//  dihedrals, each row consists of an atomList
		//  that can be passed to the setTorsion
	int strNumAtomsDistal[][] = null;
		// the number of atoms distal for each dihedral,
		//  ie. the first x entries in a row of
		//  sysDihedralDistal are valid atoms

        static double maxMovement = 9.0;
		// maximum degrees by which a torsion can
		//  cumulatively change

        int strDihedToResNum[][] = null;
		// mapping from dihedral number to residue
		//  number

        int[] numStrDihedrals = null;
        double idealDihedVal[][] = null;//Ideal values for the dihedrals in their current rotamers
	int MAX_NUM_ATOMS_DISTAL = 30;

	boolean doDihedEnergy = false;
	// If true then dihedral energies are computed and
	//  used in minimization


        boolean minimizePerturbations = true;

        int totalTransRotStrands = 0;
        boolean transRotStrands[];//which strands can rotate and translate

        
        
        
    public ContSCObjFunction(Molecule theM, int numStrands, EnergyFunction ef,
		StrandRotamers strandRotamers[], int curAANum[], boolean doDihedral, boolean[] rotTransStrands){



                // snag the local variables
		m = theM;
		efunc = ef;
		doDihedEnergy = doDihedral;

		strandRot = strandRotamers;

		numberOfStrands = numStrands;

                
                if(rotTransStrands==null){
                    //by default, all strands that can translate & rotate
                    //will be allowed to do so during minimization
                    
                    transRotStrands = new boolean[numberOfStrands];
                    for(int i=0; i<numberOfStrands;i++){
			if(m.strand[i].rotTrans){
				transRotStrands[i] = true;
                        }
                    }
                }
                else
                    transRotStrands = rotTransStrands;


                //Set up dihedral stuff
		numStrDihedrals = new int[numberOfStrands];

		for(int i=0;i<numberOfStrands;i++){
			for(int j=0;j<m.strand[i].numberOfResidues;j++){
				if(m.strand[i].residue[j].flexible)
					MAX_NUM_ATOMS_DISTAL = Math.max(MAX_NUM_ATOMS_DISTAL, m.strand[i].residue[j].numberOfAtoms);
			}
		}

		for(int i=0;i<numberOfStrands;i++){
			for(int j=0;j<m.strand[i].numberOfResidues;j++){
				if(m.strand[i].residue[j].flexible)
					numStrDihedrals[i] += strandRot[i].rl.getNumDihedrals(curAANum[m.strand[i].residue[j].moleculeResidueNumber]);
			}
		}

		strDihedralAtNums = new int[numberOfStrands][][];
		strDihedralDistal = new int[numberOfStrands][][];
		strNumAtomsDistal = new int[numberOfStrands][];
		strDihedToResNum  = new int[numberOfStrands][];

		for(int i=0;i<numberOfStrands;i++){
			strDihedralAtNums[i] = new int[numStrDihedrals[i]][4];
			strDihedralDistal[i] = new int[numStrDihedrals[i]][MAX_NUM_ATOMS_DISTAL];
			strNumAtomsDistal[i] = new int[numStrDihedrals[i]];
			strDihedToResNum[i]  = new int[numStrDihedrals[i]];
		}

		// Count number of flexible residues
		numFlexRes = new int[numberOfStrands];
		for(int str=0; str<numberOfStrands;str++){
			for(int i=0;i<m.strand[str].numberOfResidues;i++){
				if(m.strand[str].residue[i].flexible)
					numFlexRes[str]++;
                        }
		}


		totalFlexRes = 0;
		for(int i=0; i<numberOfStrands;i++)
			totalFlexRes += numFlexRes[i];

		totalTransRotStrands = 0;
		for(int i=0; i<numberOfStrands;i++){
			//if(m.strand[i].rotTrans){
                        if(transRotStrands[i]){
				totalTransRotStrands++;
                        }
                }

		/*flexResAtomList = new int[totalFlexRes+totalTransRotStrands][MAX_NUM_ATOMS_DISTAL];
		flexResListSize = new int[totalFlexRes+totalTransRotStrands];*/

		// Now determine which atoms are involved with each system dihedral
		//int curTransRotInd = -1;
		//int strResNumPP = -1;
		int curNumFlex = 0;

		for(int str=0; str<numberOfStrands;str++){
		int atoms[] = new int[4];
		int tmpAtomList[] = new int[MAX_NUM_ATOMS_DISTAL];
		int at3num = -1, at2num = -1;
		int numDihed = 0;


			//if(m.strand[str].rotTrans){
                        /*if(transRotStrands[str]){
				curTransRotInd++;
				strResNumPP = totalFlexRes+curTransRotInd;*/
				/*flexResListSize[strResNumPP] = theM.strand[str].numberOfAtoms;
				flexResAtomList[strResNumPP] = new int[flexResListSize[strResNumPP]];*/
			//}
			Residue localRes = null;
			//StrandRotamers localRH = null;
			//int prevNumAtoms = 0;
			for(int i=0;i<m.strand[str].numberOfResidues;i++){
				localRes = m.strand[str].residue[i];
				/*if(m.strand[str].rotTrans){
					for(int k=0;k<localRes.numberOfAtoms;k++){
						flexResAtomList[strResNumPP][k+prevNumAtoms] = localRes.atom[k].moleculeAtomNumber;
					}
				}*/
				//prevNumAtoms += localRes.numberOfAtoms;
				if(localRes.flexible){

					for(int j=0;j<strandRot[str].rl.getNumDihedrals(curAANum[theM.strand[str].residue[i].moleculeResidueNumber]);j++){
						strDihedToResNum[str][numDihed] = curNumFlex;
						atoms = strandRot[str].rl.getDihedralInfo(m,str,i,j);
					// note: atoms are residue relative numbering
						strDihedralAtNums[str][numDihed][0] = localRes.atom[atoms[0]].moleculeAtomNumber;
						strDihedralAtNums[str][numDihed][1] = localRes.atom[atoms[1]].moleculeAtomNumber;
						strDihedralAtNums[str][numDihed][2] = localRes.atom[atoms[2]].moleculeAtomNumber;
						strDihedralAtNums[str][numDihed][3] = localRes.atom[atoms[3]].moleculeAtomNumber;

					at2num = localRes.atom[atoms[2]].moleculeAtomNumber;
					at3num = localRes.atom[atoms[3]].moleculeAtomNumber;
					for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
						tmpAtomList[k]=1;
					}
					tmpAtomList[atoms[1]]=0;
					tmpAtomList[atoms[2]]=0;
						strandRot[str].getAtomsMoreDistal(m,localRes.moleculeResidueNumber,m.atom[at2num],tmpAtomList);
						strNumAtomsDistal[str][numDihed]=0;
					for(int k=0;k<MAX_NUM_ATOMS_DISTAL;k++){
						if ((tmpAtomList[k]==2) && (k != atoms[3])){
								strDihedralDistal[str][numDihed][strNumAtomsDistal[str][numDihed]]=localRes.atom[k].moleculeAtomNumber;
								strNumAtomsDistal[str][numDihed] += 1;
						}
					}
					/*if (j==0){
						// If this is the first dihedral for this residue snag information
						flexResListSize[curNumFlex] = strNumAtomsDistal[str][numDihed]+1;
						flexResAtomList[curNumFlex][0] = at3num;
							for(int k=0;k<strNumAtomsDistal[str][numDihed];k++){
								flexResAtomList[curNumFlex][k+1] = strDihedralDistal[str][numDihed][k];
						}
					}*/
					numDihed++;
				}
				curNumFlex++;
				}
			}
		}



		// If computing dihedral energies initialize them
		if (doDihedEnergy){//This assumes ef.getAmber96ext() is well defined (since dihedral energies are made to be used with force field energies).  There will be an error otherwise
                    de = new DihedralEnergy(m, numberOfStrands, strDihedralAtNums, numStrDihedrals, efunc.getAmber96ext());
                    efunc = efunc.addTerm(de);
		}


                calcDOFNums(m.DOFs);//Fill in DOFNums, numDOFs, flexToMolResMap
                setupPartialComp();//Set up partial computations for energy function
                
                specialBoxConstraints = new double[numDOFs][];//Initialize all special constraints to null

                curDOFVals = DoubleFactory1D.dense.make(numDOFs);
                updateIdealDihedrals();
    }


    private void calcDOFNums(DegreeOfFreedom[] DOFList){
        //Calculate DOFNums, numDOFs, and flexToMolResMap

        flexToMolResMap = new int[totalFlexRes];

        //Figure out which perturbations to minimize over
        HashSet<Integer> pertSet=new HashSet<Integer>();

        int flexResCount = 0;


        if(minimizePerturbations){
            //Find the perturbations affecting flexible residues
            for(int i=0;i<numberOfStrands;i++){
                    for(int j=0;j<m.strand[i].numberOfResidues;j++){

                            Residue localResidue = m.strand[i].residue[j];

                            if(localResidue.flexible){

                                for(int a=0;a<localResidue.perts.length;a++){//Add this residue's perturbations to the perturbation set, if they're not already in it

                                    if( !pertSet.contains(localResidue.perts[a]) )
                                        pertSet.add(localResidue.perts[a]);
                                }

                                flexToMolResMap[flexResCount] = localResidue.moleculeResidueNumber;
                                flexResCount++;
                            }
                    }
            }
        }
        

        //calculate numDOFs
        numDOFs = getNumTotDihed() + pertSet.size() + 6*totalTransRotStrands;
        DOFNums = new int[numDOFs];

        int DOFCount = 0;
        int dihedNum = 0;
        int oldRes = -1;

        //Get DOF numbers for dihedrals
        for(int str=0; str<numberOfStrands;str++){
                for(int j=0;j<numStrDihedrals[str];j++) {

                    int curRes = strDihedToResNum[str][j];
                    int strResNum = m.atom[strDihedralAtNums[str][j][0]].strandResidueNumber;

                    if(curRes != oldRes)
                        dihedNum = 0;

                    for( DegreeOfFreedom DOF : DOFList ){
                        if( DOF.type == DegreeOfFreedom.SCDIHEDRAL ){
                            if(DOF.strResNum == strResNum){
                                if( DOF.strandNum == str ){
                                    if( DOF.dihedNum == dihedNum ){//If this DOF is the current dihedral
                                        DOFNums[DOFCount] = DOF.DOFNum;
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    oldRes = curRes;
                    dihedNum++;
                    DOFCount++;
                }
        }


        //perturbations
        for( DegreeOfFreedom curDOF : DOFList ){
            if(curDOF.type==DegreeOfFreedom.PERTURBATION){
                if(pertSet.contains(curDOF.pertNum) && curDOF.isContinuous()){//exclude discrete perturbations
                    DOFNums[DOFCount] = curDOF.DOFNum;
                    DOFCount++;
                }
            }
            else if(curDOF.type == DegreeOfFreedom.STRRIGIDMOTION){//Any strand rotation/translation can be minimized over...
                if(transRotStrands[curDOF.strandNum]){
                    DOFNums[DOFCount] = curDOF.DOFNum;
                    DOFCount++;
                }
            }
        }

        if(DOFCount<numDOFs){//Some of the perturbations weren't continuous, so we can reduce numDOFs accordingly
            numDOFs = DOFCount;
            int newDOFNums[] = new int[numDOFs];
            System.arraycopy(DOFNums,0,newDOFNums,0,numDOFs);
            DOFNums = newDOFNums;
        }
    }



    private void setupPartialComp(){
        //Setup partial computations for energy function
        //partialSets[0] will be the set of residues affected by any DOF minimized over
        //partialSets[a] for a>0 will be the set of residues affected by DOF #(a-1)
        //being minimized over
        int partialSets[][] = new int[numDOFs+1][];

        HashSet<Integer> allAffectedRes = new HashSet<Integer>();
        //molecule-based numbers of all residues moved by continuous DOFs

        for(int dof=0; dof<numDOFs; dof++){
            DegreeOfFreedom curDOF = m.DOFs[DOFNums[dof]];
            if(curDOF.type==DegreeOfFreedom.STRRIGIDMOTION){
                //Moves all the residues in its strand, whether flexible or not
                int str = curDOF.strandNum;
                int numRes = m.strand[str].numberOfResidues;
                partialSets[dof+1] = new int[numRes];
                for(int res=0; res<numRes; res++){
                    int molResNum = m.strand[str].residue[res].moleculeResidueNumber;
                    partialSets[dof+1][res] = molResNum;
                    allAffectedRes.add(molResNum);
                }
            }
            /*else{//Other DOF types only affect flexible residues
                int numRes = curDOF.flexResAffected.length;
                partialSets[dof+1] = new int[numRes];
                for(int res=0; res<numRes; res++){
                    int molResNum = flexToMolResMap[curDOF.flexResAffected[res]];
                    partialSets[dof+1][res] = molResNum;
                    allAffectedRes.add(molResNum);
                }
            }*
             * //
             * //Can't use curDOF.flexResAffected because those are all residues flexible in the design
            //and we might only be minimizing w.r.t. a pair or a single residue
             */
            else if(curDOF.type==DegreeOfFreedom.SCDIHEDRAL){
                int molResNum = m.strand[curDOF.strandNum].residue[curDOF.strResNum].moleculeResidueNumber;
                partialSets[dof+1] = new int[] { molResNum };
                allAffectedRes.add(molResNum);
            }
            else if(curDOF.type==DegreeOfFreedom.PERTURBATION){
                HashSet<Integer> pertAffectedRes = new HashSet<Integer>();

                for(Residue res : m.residue){
                    //if(res.flexible){//This actually shouldn't change anything since only flexible or template (so unperturbed) residues have energies counted
                        for(int p : res.perts){
                            if(p==curDOF.pertNum){
                                pertAffectedRes.add(res.moleculeResidueNumber);
                                allAffectedRes.add(res.moleculeResidueNumber);
                            }
                        }
                    //}
                }

                partialSets[dof+1] = new int[pertAffectedRes.size()];
                int count = 0;
                for( int molResNum : pertAffectedRes ){
                    partialSets[dof+1][count] = molResNum;
                    count++;
                }
            }
        }

        partialSets[0] = new int[allAffectedRes.size()];
        int count = 0;
        for( int molResNum : allAffectedRes ){
            partialSets[0][count] = molResNum;
            count++;
        }

        /*
        //The first row contains all flexible residues
        //Then we have each of the flexible residues separately
        int partialSets[][] = new int[totalFlexRes+1][];

        partialSets[0] = flexToMolResMap;
        for(int a=0; a<totalFlexRes; a++)
            partialSets[a+1] = new int[] { flexToMolResMap[a] };*/


        efunc.setupPartialComputation(partialSets);
        //a96ff.setupPartialArrays(totalFlexRes+totalTransRotStrands,MAX_NUM_ATOMS_DISTAL,flexResAtomList,
        //			flexResListSize);
    }


    public int getNumTotDihed(){
		int numDihedrals = 0;
		for(int str=0;str<numberOfStrands;str++)
			numDihedrals+=numStrDihedrals[str];

		return numDihedrals;
    }


    public int getNumDOFs(){
        return numDOFs;
    }


/*    public int getNumDOFs(){
        return getNumTotDihed();
        //Add rot, trans, perturbations later!!!
    }*/

    //This function calculates the idealDihedVal based on the current dihedrals
    //they are updated here and in de
    public void updateIdealDihedrals(){

        if(idealDihedVal == null){//allocate
            idealDihedVal = new double[numberOfStrands][];
            for(int str=0; str<numberOfStrands;str++)
                idealDihedVal[str] = new double[numStrDihedrals[str]];
        }

        //Get current values of the dihedrals
        for(int str=0; str<numberOfStrands;str++){

                for(int j=0;j<numStrDihedrals[str];j++) {

                    double dihedVal = m.getTorsion(strDihedralAtNums[str][j][0],strDihedralAtNums[str][j][1],
                            strDihedralAtNums[str][j][2],strDihedralAtNums[str][j][3]);

                    dihedVal = enforceAngleRangeConvention(dihedVal);//Make sure dihedVal follows the angle range convention ((-180,180])

                    idealDihedVal[str][j] = dihedVal;
                }
        }

        if(doDihedEnergy)
            de.idealDihedVal = idealDihedVal;
    }
    

    //Get constraints on the degrees of freedom
    //This should be called before minimization
    public DoubleMatrix1D[] getConstraints(){

        //Re-initialize strand translation/rotation info
        //This needs to be done whenever we are setting up the constraints
        for(int i=0; i<numberOfStrands;i++){
                //if(m.strand[i].rotTrans){
                if(transRotStrands[i]){
                        //totalTransRotStrands++;
                        //Initialize translation/rotation to 0
                        m.strand[i].curTrans = new double[3];
                        m.strand[i].curRotAngles = new double[3];
                        m.strand[i].curRotMatrix = new RotMatrix().identity();
                        m.strand[i].strStartCOM = m.getStrandCOM(i);
                }
        }



        DoubleMatrix1D DOFmin = DoubleFactory1D.dense.make( numDOFs );
        DoubleMatrix1D DOFmax = DoubleFactory1D.dense.make( numDOFs );


        for(int dof=0; dof<numDOFs; dof++){//Loop over minimizing DOFs

            if( specialBoxConstraints[dof] == null ){//Default constraints

                DegreeOfFreedom curDOF = m.DOFs[DOFNums[dof]];

                if(curDOF.type == DegreeOfFreedom.SCDIHEDRAL){
                    int str = getDihedStrNum(dof);
                    int j = getDihedStrBasedNum(dof,str);
                    DOFmin.set(dof, idealDihedVal[str][j] - maxMovement );
                    DOFmax.set(dof, idealDihedVal[str][j] + maxMovement );
                }
                else if(curDOF.type == DegreeOfFreedom.PERTURBATION){
                    Perturbation pert = curDOF.pert;
                    int curState = pert.curState;//Current state of perturbation
                    DOFmin.set(dof, pert.minParams[curState]);
                    DOFmax.set(dof, pert.maxParams[curState]);
                }
                else if(curDOF.type == DegreeOfFreedom.STRRIGIDMOTION){
                    //curDOF just has one interval so we bound its values in that interval
                    DOFmin.set(dof, curDOF.intervals[0][0]);
                    DOFmax.set(dof, curDOF.intervals[0][1]);
                }
                else{
                    System.err.println("ERROR: Cannot minimize over DOF type "+curDOF.type);
                    System.exit(1);
                }
            }
            else{
                DOFmin.set(dof, specialBoxConstraints[dof][0]);
                DOFmax.set(dof, specialBoxConstraints[dof][1]);
            }

            

        }


        /*int DOFNum = 0;

        for(int str=0; str<numberOfStrands;str++){
                for(int j=0;j<numStrDihedrals[str];j++) {

                    if( specialBoxConstraints[DOFNum] == null ){//Default constraints
                        DOFmin.set(DOFNum, idealDihedVal[str][j] - maxMovement );
                        DOFmax.set(DOFNum, idealDihedVal[str][j] + maxMovement );
                    }
                    else{
                        DOFmin.set(DOFNum, specialBoxConstraints[DOFNum][0]);
                        DOFmax.set(DOFNum, specialBoxConstraints[DOFNum][1]);
                    }
                        

                    DOFNum++;
                }
        }*/


        //This stuff was from the beginning of SimpleMinimizer.minimize()...
        //could use for str/rot
        /*
          		double strTorque[] = new double[3];
		double strTrans[] = new double[3];

		strStartCOM = new double[numberOfStrands][3];
		strCurCOM = new double[numberOfStrands][3];
		for(int str=0;str<numberOfStrands;str++){
			strStartCOM[str] = m.getStrandCOM(str);
			strCurCOM[str][0] = strStartCOM[str][0];
			strCurCOM[str][1] = strStartCOM[str][1];
			strCurCOM[str][2] = strStartCOM[str][2];
		}
         */



        return new DoubleMatrix1D[] { DOFmin, DOFmax };
    }



    //Value and gradient at a given point (specified as values for all DOFs)
    public double getValue(DoubleMatrix1D x){

        setDOFs(x);

        /*
        //COULD ACCELERATE THIS WITH PARTIAL ARRAYS TOO
        if(totalTransRotStrands>0)//even the template might move...
            return efunc.getEnergy();*/
        
        double energy = efunc.getEnergy(0);//This should account for dihedral energies
        //The 0 indicates energy terms involving all residues affected by DOFs will be used
        //secondEnergy[0] += computeOneDihedEnergyDiff(strNumber,dihedNumForCur,cumulStep+stepSize);
        
        
        return energy;
    }

   

    //Set the DOFs for the molecule in the state indicated by x
    public void setDOFs(DoubleMatrix1D x){

        for(int dof=0; dof<numDOFs; dof++)
            setDOF(dof,x.get(dof));
    }


    public void setDOF(int dof, double val){

        curDOFVals.set(dof, val);
        DegreeOfFreedom curDOF = m.DOFs[DOFNums[dof]];

        if(curDOF.type == DegreeOfFreedom.SCDIHEDRAL){
            int str = getDihedStrNum(dof);
            int j = getDihedStrBasedNum(dof,str);

            //Set the dihedal value appropriately
            m.setTorsion(strDihedralAtNums[str][j][0],strDihedralAtNums[str][j][1],
                                strDihedralAtNums[str][j][2],strDihedralAtNums[str][j][3],
                                val,strDihedralDistal[str][j],strNumAtomsDistal[str][j],false);
        }
        else if(curDOF.type == DegreeOfFreedom.PERTURBATION){
            Perturbation pert = curDOF.pert;
            pert.changePerturbationParameter((double)val);
        }
        else if(curDOF.type == DegreeOfFreedom.STRRIGIDMOTION){
            //Translation: val = translation relative to starting point in given dimension
            //Rotation: val = rotation in degrees in given dimension

            int str = curDOF.strandNum;
            Strand curStrand = m.strand[str];
            

            if(curDOF.rigidDOFNum<3){//translation

                double delta = (double) ( val - curStrand.curTrans[curDOF.rigidDOFNum] );
                //new translation that needs to happen in the specified dimension

                curStrand.curTrans[curDOF.rigidDOFNum] = val;

                switch(curDOF.rigidDOFNum){
                    case 0:
                        m.translateStrand(str, delta, 0, 0, false);
                        break;
                    case 1:
                        m.translateStrand(str, 0, delta, 0, false);
                        break;
                    case 2:
                        m.translateStrand(str, 0, 0, delta, false);
                        break;
                }
            }
            else{//rotation

                RotMatrix r = new RotMatrix();
                double curMatrix[][] = curStrand.curRotMatrix;//Current rotation matrix

                double angles[] = curStrand.curRotAngles;
                angles[curDOF.rigidDOFNum-3] = (double)val;//Update angles
                

                double newRotMatrix[][] = new double[3][3];
                r.getRotMatrix(1,0,0,angles[0],newRotMatrix);

                double dimRot[][] = new double[3][3];//Rotation matrix in each additional dimension
                r.getRotMatrix(0,1,0,angles[1],dimRot);
                newRotMatrix = r.multiplyMatrices(dimRot, newRotMatrix);
                r.getRotMatrix(0,0,1,angles[2],dimRot);
                newRotMatrix = r.multiplyMatrices(dimRot, newRotMatrix);

                curStrand.curRotMatrix = newRotMatrix;

                double changeMatrix[][] = r.multiplyMatrices( newRotMatrix, r.transpose(curMatrix) );//apply the transpose (i.e. inverse) of the current rotation matrix, then the new rotation matrix
                double curStrCenter[] = r.add( curStrand.strStartCOM , curStrand.curTrans );//center to rotate around
                m.rotateStrand(str, changeMatrix, curStrCenter[0], curStrCenter[1], curStrCenter[2], false);
            }
        }
        else{
            System.err.println("ERROR: Cannot minimize over DOF type "+curDOF.type);
            System.exit(1);
        }
    }
    
    
    
    public double getValForDOF(int dof,double val){
        //Computes energy terms after setting dof at val and keeping all other DOFs the same
        //(for example, they can have been set by setDOFs)
        //the energy returned must include all terms that depend on dof

        //DegreeOfFreedom curDOF = m.DOFs[DOFNums[dof]];

        setDOF(dof,val);

        
        double energy = efunc.getEnergy(dof+1);
        
        return energy;
        /*
        //THESE ENERGIES COULD BE MUCH IMPROVED BY HAVING SPECIAL PARTIAL ARRAYS FOR EACH DOF!
        if(curDOF.type == DegreeOfFreedom.SCDIHEDRAL){
            int str = getDihedStrNum(dof);
            int j = getDihedStrBasedNum(dof,str);

            //Return the energy terms involving the residue whose dihedral this is
            return efunc.getEnergy(strDihedToResNum[str][j]+1);
        }
        else if(curDOF.type == DegreeOfFreedom.PERTURBATION){
            return efunc.getEnergy(0);
        }
        else{//might move even the template?
            return efunc.getEnergy();
        }*/
    }


    public double getInitStepSize(int dof){
        //initial step size to use for dof

        DegreeOfFreedom curDOF = m.DOFs[DOFNums[dof]];
        
        //We'll peg this at a quarter of the mesh width, so even for a 4-component
        //diagonal DOF we can fit a step inside a subvoxel
        return curDOF.getMeshWidth()/4;

    }
    
    
        @Override
    public boolean isDOFAngle(int dof){
        DegreeOfFreedom curDOF = m.DOFs[DOFNums[dof]];
        return curDOF.isDOFAngle();
    }


    public static double enforceAngleRangeConvention(double theta){
        //Enforces angle range convention for ideal dihedrals
        //There can be a problem when it's supposed to be 180 but ends up being measured as -179.9999 or something
        //(A* will assume 180 and FuncLBTerm nodes can be off by 360 degrees)
        //Since the ideal dihedral comes from a standard rotamer, it is expected to be an integer
        //so we will make it be in the interval [-179.5,180.5)

        double shiftCycles = Math.floor((theta+179.5)/360);
        double shift = -shiftCycles*360;//We need to shift the nodes up by shift*360

        return theta+shift;
    }



    private int getDihedStrNum(int dihedNum){//Figure out which strand a dihedral is on
        //given a molecule-based number
        int base=0;
        for(int str=0; str<numberOfStrands; str++ ){
            base += numStrDihedrals[str];
            if( base > dihedNum )
                return str;
        }
        System.err.println("ERROR: Can't find a strand for dihedral # "+dihedNum);
        System.exit(1);
        return -1;
    }


    private int getDihedStrBasedNum(int dihedNum, int str){
        //Get a strand-based dihedral number from a molecule-based one and a strand number
        for(int str2=0; str2<str; str2++ ){//Subtract off numbers of dihedrals on earlier strands
            dihedNum -= numStrDihedrals[str2];
        }
        return dihedNum;
    }


}
