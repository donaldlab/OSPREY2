
import cern.colt.matrix.DoubleMatrix1D;
import java.util.ArrayList;
import java.util.HashMap;

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
//	MSAStar.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//        PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
//	  MAH        Mark A. Hallen	   Duke University         mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ivelin Georgiev (2004-2009)
 * 
 * The single-mutation-sequence A* uses the steric filter object; the steric filter cannot be used with multiple-mutation-sequence A*.
 * 
 * The multiple-mutation-sequence A* can allow the generation only of sequences with up to numMaxChanges mutations from the wildtype.
 * 
 * The parameters from the calling program need to be modified to fit the reduced matrices here:
 * 		only the entries for the possible mutations (AA assignments) are stored in the matrices;
 * 
 * The expansion queue is returned after each call, so that the next call can start with the
 * 		saved expansion queue: this allows to return the next-lowest-energy conformation
 * 
 * If A* is not able to find a non-pruned node for a given level (all of the possible nodes for that level
 * 		have been pruned by DEE), then it returns immediately, with the node for that level set to -1
*/

/**
 * Uses A* search for single or multiple mutation sequences simultaneously to return the minimum-energy conformation; 
 * 		each consecutive run returns the next lowest-energy conformation, in order.
 * 
 */
public class MSAStar {

	//number of residues under consideration
	private int numTreeLevels;
	
	//number of rotamers possible for each residue (given by the assigned AA type)
	private int numNodesForLevel[] = null;
		
	//the offset in the array index for each level
	private int nodeIndexOffset[] = null;
		
	//eliminated rotamers at residue i, for all residues
	//private boolean eliminatedNodesAtLevel [] = null;
	
	//the current sequence: the number of the corrsponding rotamer for each level assigned so far 
	private int curConf[] = null;

	//the reduced min pairwise energy matrix
	private ReducedEnergyMatrix pairwiseMinEnergyMatrix = null;
	
	//the leaf nodes visible from the expansions
	private ExpansionQueue curExpansion;
	
	int numExpanded = 0;
        int numFS = 0;
	
	int topL = 0;
	int numTopL = 0;
	
	//the steric check filter (*only* for single-sequence A*)
	StericCheck stericF = null;

        boolean[][] splitFlags = null;//Reduced-size array indicating which pairs of rotamers are pruned
        boolean[][][] tripleFlags = null;//Which triples are pruned

        //Use DEEPer
        boolean doPerturbations;

        EPICSettings es;
        
        //if using EPIC we'll need these:
        CETMatrix CETM = null;//the series
        DegreeOfFreedom DOFList[] = null;//the degrees of freedom for the system

        
        
        
        
        //stuff for SVE minimization
        //can be null if no EPIC
        Molecule m;
        StrandRotamers[] strandRot;
        int mutRes2Strand[];
        int strandMut[][];
        int mutRes2StrandMutIndex[];
        
        

	//constructor
	/*
	 * We assume that the parameters supplied (energy and DEE information) have already been modified
	 * 		to consider only the residues and rotamers for the possible mutations; i.e., the matrices are of
	 * 		reduced size
	 */
	MSAStar (int treeLevels, int numRotForRes[], ReducedEnergyMatrix arpMatrixRed, StericCheck stF, boolean[][] spFlags,
                boolean[][][] tripFlags, boolean doPerts, EPICSettings eset, 
                Molecule molec, StrandRotamers[] sr, int mr2s[], int sm[][], int mr2smi[]){
	
                es = eset;
                
		numTreeLevels = treeLevels;
                pairwiseMinEnergyMatrix = arpMatrixRed;

		
		numNodesForLevel = new int [numTreeLevels];
		nodeIndexOffset = new int [numTreeLevels];
		int offset = 0;

		for (int i=0; i<numTreeLevels; i++){
			nodeIndexOffset[i] = offset;
			numNodesForLevel[i] = numRotForRes[i];
			offset += numNodesForLevel[i];
		}

		
		//the current expansion list
		curExpansion = new ExpansionQueue();
		
		//the current conformation
		curConf = new int [numTreeLevels];
		for (int i=0; i<numTreeLevels; i++){
			curConf[i] = -1;
		}
		
		stericF = stF;

                splitFlags = spFlags;
                tripleFlags = tripFlags;

                doPerturbations = doPerts;
                
                
                try{
                    m = (Molecule)KSParser.deepCopy(molec);
                    strandRot = (StrandRotamers[])KSParser.deepCopy(sr);
                    for(Residue r : m.residue)
                        r.flexible = false;//we'll make only the assigned residues flexible
                }
                catch(Exception e){
                    throw new RuntimeException("ERROR: molecule or StrandRotamers[] deepCopy failed in A* initialization!");
                }
                mutRes2Strand = mr2s;
                strandMut = sm;
                mutRes2StrandMutIndex = mr2smi;
	}
	
	//Find the lowest-energy conformation and return an array with the number of the corresponding
	//		chosen rotamer for each residue;
	//		the mapping to the rotamer information is done by the calling procedure;
	//		the chosen conformation should be marked so that the next call to AStar will return the
	//			conformation with the next lowest energy value, and so on
	/*
	 * Look at the minimum value in the expansion list; determine the node level corresponding to this value;
	 * 		expand the node into the next level; update the expansion list (adding the new nodes and deleting
	 * 		the expanded node) and determine the f(n)=g(n)+h(n) scores for the new nodes
	 * 
	 * To get the next lowest conformation, the state of the expansion queue is saved after each complete run
	 * 		of A*: after the lowest-energy conformation is found, the queue is saved, and A* returns this
	 * 		conformation. To find the second conformation, A* runs on the saved queue, and this is repeated
	 * 		for all subsequent conformations
	 */
	public int[] doAStar (boolean run1, int numMaxChanges, int nodesDefault[], boolean prunedNodes[],
			StrandRotamers strandRot[], String strandDefault[][], int numForRes[], int strandMut[][], boolean singleSeq, 
			int mutRes2Strand[], int mutRes2MutIndex[]){

		int curLevelNum = 0;
		double hScore;
		double gScore;
		double fScore;
		QueueNode expNode = null;
		QueueNode newNode = null;
		
		for (int i=0; i<numTreeLevels; i++){ //initialize for this run
			curConf[i] = -1;
		}

		if (run1) {//if this is the first run of A*, then we need to set-up the empty expansion queue
			
                            //initially, we are at level zero; all nodes at that level are visible;
                            //	compute their f(n) score and add to the expansion queue
                            for (int curNode=0; curNode<numNodesForLevel[0]; curNode++){

                                    //if ((stericF==null)||(stericF.checkAllowedSteric(0,curConf,curNode))){//do not do a steric check if backbone minimization

                                            curConf[0] = curNode;//this is the only node in the current conformation

                                            //compute f for the current node
                                            gScore = gCompute (curLevelNum, curConf);
                                            hScore = hCompute (curLevelNum, curConf);
                                            fScore = gScore + hScore;

                                            //create a queueNode with the corresponding information
                                            newNode = new QueueNode (curNode, 0, curConf, fScore);

                                            
                                            //DEBUG!!!!  To use when checking minimization for a particular conformation
                                            /*if(CETM==null)
                                                CETM = (CETMatrix)KSParser.readObject("1l7mCETM_COM.dat");
                                            newNode.level = 6;
                                            newNode.confSoFar = new int[] {0,0,0,0,0,0,0};
                                            double E = FSTerm(newNode);
                                            */
                                            
                                            
                                            
                                            //insert in the expansion list
                                            curExpansion.insert(newNode);
                                    //}
                            }
                            if (curConf[0]==-1) //no sterically allowed nodes at the first residue, so no possible conformations				
                                    return curConf;	
		}

		boolean done = false;		
		//start forming the conformation by expanding the lowest-valued node and updating the expansion queue
		/*
		 * While not at the last level
		 * 	For the current minimum node in the queue
		 * 		For all possible nodes at the next level
		 * 			If the steric is allowed
		 * 				Compute the f(n) scores; Add to the queue
		 * 		Delete the expanded node from the queue
		 */
		while (!done) {	
			
			if(!run1){
				for (int i=0; i<numTreeLevels; i++) //reinitialize for each consecutive node to be expanded
					curConf[i] = -1;
			}
			
			expNode = curExpansion.curFront;//get the current min node
			
			if (expNode==null){//the queue is empty
				return curConf; //so return a sequence of -1's to flag the stop of the search
			}

                        if(es.useEPIC){
                            //We will only expand terms that have the FS term included
                            //until our lowest lower bound contains an FS term, we will
                            //compute this term for the lowest-bounded nodes we get
                            //and then throw them back in the queue (this is to minimize the number of
                            //times we have to compute the FS term)
                            while(!expNode.FSTermIncluded){
                                curExpansion.delete(expNode);
                                expNode.fScore += FSTerm(expNode);
                                expNode.FSTermIncluded = true;
                                curExpansion.insert(expNode);
                                expNode = curExpansion.curFront;
                            }
                            //Now expNode has an FS term so we can expand it
                        }
				
                        printState(expNode);
                        
                            for (int i=0; i<=expNode.level; i++){//get the corresponding conformation leading to this node
                                    curConf[i] = expNode.confSoFar[i];
                            }
                            curLevelNum = expNode.level;//get the corresponding level

                            //if the current min node is at the last level, we have found a full conformation
                            if (curLevelNum==numTreeLevels-1){
                                    done = true;
                            }
                            else {//we are not at the last level, so we can expand the current node to the next level

                                    curLevelNum++;	//since the new nodes are at the next level, increment the current level number


                                    int numChanges = 0;
                                    int curPruningInd = 0;
                                    if (!singleSeq) { //multiple sequences
                                            //Check the number of differences from the default node sequence: should not exceed numMaxChanges
                                            for (int i=0; i<curLevelNum; i++){ //first, get the actual conf for the assigned curConf
                                                    int curNodeInd = 0;
                                                    int actualNode = -1;
                                                    for (int curRot=0; curRot<numForRes[i]; curRot++){
                                                            if (!prunedNodes[curPruningInd]){
                                                                    if (curNodeInd==curConf[i])
                                                                            actualNode = curRot;
                                                                    curNodeInd++;
                                                            }
                                                            curPruningInd++;
                                                    }
                                                    int index = getAAIndex(actualNode,i,strandDefault,strandRot,strandMut,mutRes2Strand,mutRes2MutIndex);
                                                    if (index!=nodesDefault[i])
                                                            numChanges++;
                                            }
                                    }

                                    for (int curNode=0;curNode<numNodesForLevel[curLevelNum];curNode++){

                                            int numChanges2 = numChanges;
                                            if (!singleSeq) //multiple sequences, so compute the difference from WT
                                                    numChanges2 = getNewNumChanges(curPruningInd,curLevelNum,numForRes,prunedNodes,curNode,strandDefault,strandRot,strandMut,nodesDefault,numChanges,mutRes2Strand, mutRes2MutIndex);


                                            if ((singleSeq)||(numChanges2<=numMaxChanges)){ //continue only if (singleSeq) or (num differences does not exceed the max one)

                                                    //if ((stericF==null)||(stericF.checkAllowedSteric(curLevelNum,curConf,curNode))){//do not do a steric check if backbone minimization
                                                    if( ! isPruned(nodeIndexOffset[curLevelNum] + curNode, curConf, curLevelNum-1) ){//Don't make this new node if its conformation contains pruned pairs or triples


                                                            curConf[curLevelNum] = curNode;//add curNode to the conformation so far

                                                            //compute f for the current node
                                                            gScore = gCompute (curLevelNum, curConf);
                                                            hScore = hCompute (curLevelNum, curConf);
                                                            fScore = gScore + hScore;

                                                            //create a queueNode with the corresponding information
                                                            newNode = new QueueNode (curNode, curLevelNum, curConf, fScore);

                                                            if(es.useEPIC && expNode.optFSPoint!=null)//inherit optimal DOF values from parent, to use as initial values if we need to minimize
                                                                newNode.optFSPoint = expNode.optFSPoint.copy();

                                                            //insert in the expansion list
                                                            curExpansion.insert(newNode);
                                                    }
                                                    //}
                                            }
                                    }
                            }
                            

                        //delete the expanded node from the queue, since it has already contributed
                        curExpansion.delete(expNode);
			//}
		}

                System.out.println("A* returning conformation; lower bound = "+expNode.fScore+" nodes expanded: "+numExpanded+" FS terms evaluated: "+numFS);
		
                return curConf;
	}
	
	private int getAAIndex(int rotIndex, int curRes, String strandDefault[][], StrandRotamers strandRot[], int strandMut[][], int mutRes2Strand[], int mutRes2MutIndex[]){
		
		int str = mutRes2Strand[curRes];
		int strResNum = strandMut[str][mutRes2MutIndex[curRes]];
		int rotSum = 0;
		for (int i=0; i<strandRot[str].getNumAllowable(strResNum); i++){

                        int curRot;

                        if(doPerturbations)
                            curRot = ((StrandRCs)strandRot[str]).getNumRCs(strResNum, strandRot[str].getIndexOfNthAllowable(strResNum,i));
                        else{
                            curRot = strandRot[str].rl.getNumRotamers(strandRot[str].rl.getAAName(strandRot[str].getIndexOfNthAllowable(strResNum,i)));
                            if (curRot==0) //GLY or ALA
                               	curRot = 1;
                        }
			rotSum += curRot;			
			if (rotSum>rotIndex)
				return (strandRot[str].getIndexOfNthAllowable(strResNum,i));
		}
		return -1; //the AA in the last position
	}
	
	private int getNewNumChanges(int curPruningInd, int curLevelNum, int numForRes[], boolean prunedNodes[], int curNode, String strandDefault[][],
			StrandRotamers strandRot[], int strandMut[][], int nodesDefault[], int numCurChanges,int mutRes2Strand[], int mutRes2MutIndex[]){
		
		//Get the actual node number for curNode
		int actualNode = -1;
		int index = -1;
		int curPruningInd2 = curPruningInd;
		//if ( ! ((residueMap.length!=numTreeLevels)&&(curLevelNum==(numTreeLevels-1))) ) {//not (ligand present and at the ligand level)	
		int curNodeInd = 0;
		for (int curRot=0; curRot<numForRes[curLevelNum]; curRot++){
			if (!prunedNodes[curPruningInd2]){								
				if (curNodeInd==curNode)
					actualNode = curRot;
				curNodeInd++;
			}
			curPruningInd2++;
		}							
		
		index = getAAIndex(actualNode,curLevelNum,strandDefault,strandRot,strandMut,mutRes2Strand, mutRes2MutIndex);
		if (index!=nodesDefault[curLevelNum])
			numCurChanges++;
		//}
		
		return numCurChanges;
	}
	
	//Updates and prints the state of the queue
	private void printState(QueueNode expNode){
		
		numExpanded++;
		if (expNode.level+1>topL){
			topL = expNode.level+1;
			numTopL = 1;
		}
		else if (expNode.level+1==topL)
			numTopL++;
		
		if((numExpanded%1000)==0){
			System.out.print(curExpansion.numNodes+" "+expNode.fScore+" "+expNode.level+" ");
			for (int i=0;i<expNode.level;i++)System.out.print(expNode.confSoFar[i]+" ");
			System.out.println(topL+" "+numTopL);
		}
	}
//////////////////////////////////////////////////////////////////////////
	
//////////////////////////////////////////////////////////////////////////
	
	//Compute the h(n) score for the new node expanded by expNode
	//		called by doAStar(.)
	private double hCompute (int dLevel, int conf[]){
		
		double hn = 0.0f;

		for (int curLevel=dLevel+1;curLevel<numTreeLevels;curLevel++){
			hn += EnergyAtLevel(dLevel, curLevel, conf);
		}
			
		return hn;
	}
	
	//Called by hCompute(.)
	private double EnergyAtLevel(int topLevel, int curLevel, int conf[]){
		
		double minE = (double)Math.pow(10,30);
		double curE;		
		int index1;
		
		double minShellIntraE;			//formula term 1
		double sumMinPairE;			//formula term 2
		double sumMinMinPairE;		//formula term 3
		
		for (int i1=0; i1<numNodesForLevel[curLevel];i1++){		//the rotamers at j

			index1 = nodeIndexOffset[curLevel]+i1;	//the index of s at j

                        if( ! isPruned( index1, conf, topLevel ) ){
				
                            minShellIntraE = pairwiseMinEnergyMatrix.getIntraAndShellE(index1);
                            sumMinPairE = hSumMinPVE (topLevel, index1, conf);
                            sumMinMinPairE = sumMinMinPVE(curLevel, index1, conf, topLevel);

                            curE = minShellIntraE + sumMinPairE + sumMinMinPairE;
                            if (curE<minE)		//compare to the min energy found so far
                                    minE = curE;
                        }
		}
		
		return minE;
		
	}
	
	//Called by EnergyAtLevel(.)
	private double hSumMinPVE (int topLevel, int index1, int conf[]){
		
		double sum = 0.0f;
		int index2;
		
		for (int level=0; level<=topLevel; level++){
			
			index2 = nodeIndexOffset[level] + conf[level]; //the index of r at i
			
			sum += pairwiseMinEnergyMatrix.getPairwiseE( index2, index1 ); //the pairwise energy between the two nodes
		}
		
		return sum;
	}
	
	//Called by EnergyAtLevel(.)
	private double sumMinMinPVE(int jLevel, int firstIndex, int conf[], int topLevel){
		
		double sum = 0.0f;
		for (int level=jLevel+1; level<numTreeLevels; level++){
			sum += indMinMinPVE(level, firstIndex, conf, topLevel);
		}
		
		return sum;
	}
	
	//Called by sumMinMinPVE(.)
	private double indMinMinPVE (int kLevel, int firstIndex, int conf[], int topLevel){
		
		double minEn = (double)Math.pow(10,30);
		double curEn;
		int secondIndex;
		
		for (int i2=0; i2<numNodesForLevel[kLevel]; i2++){ //u at k
			
			secondIndex = nodeIndexOffset[kLevel]+i2;

                        if( ! isPruned( secondIndex, firstIndex, conf, topLevel ) ){
                            
                            curEn = pairwiseMinEnergyMatrix.getPairwiseE( firstIndex, secondIndex );
                            if (curEn<minEn)
                                	minEn = curEn;
                        }
		}
		
		return minEn;
	}
//////////////////////////////////////////////////////////////////////////
	
//////////////////////////////////////////////////////////////////////////
	
	//Compute the g(n) score for the new node expanded by expNode;
	//		called by doAStar(.)
	private double gCompute (int dLevel,int conf[]){
		
		double gn = 0.0f;
		int index1;
		
		double minShellIntraE;		//formula term 1
		double sumMinPairE;		//formula term 2
		
		for (int curLevel=0; curLevel<=dLevel; curLevel++){//copute using the formula
			
			index1 = nodeIndexOffset[curLevel] + conf[curLevel];//index of r at i
			
			minShellIntraE = pairwiseMinEnergyMatrix.getIntraAndShellE(index1);
			
			sumMinPairE = gSumMinPVE(dLevel, curLevel+1, index1, conf);
			
			gn += (minShellIntraE + sumMinPairE);
		}

		return gn;
	}


        private double FSTerm(QueueNode expNode){
            //Get contribution to g-score from fit series
            //include h-score too if useHSer
            //optDOFVals is initially from the parent, to use an initial values; update to be optimal here


            //if only want FS terms for full confs
            if( (expNode.level<numTreeLevels-1) && !(es.minPartialConfs) )
                return 0;

            CETObjFunction of = getNodeObjFunc(expNode);
            
            
            CCDMinimizer ccdMin = new CCDMinimizer(of,false);
            
            //TEMPORARILY DE-ACTIVATING OPTFSPOINT
            /*if(expNode.optFSPoint!=null && useHSer ){//we have initial values to use...
                //OUTFIT FOR CHANGING NUMBERS OF DOFS IN !USEHSER THOUGH...
                
                //if we're not getting much improvement by updating the FSTerm
                //even with the old initial values,
                //then the minimized value will not offer significant improvement over the current score
                //so don't spend time on it
                //dof_inh_lazy
                double initE = of.getValue( expNode.optFSPoint );
                if(useHSer && initE < expNode.fScore+0.1)
                    return expNode.fScore;
                
                
                ccdMin.singleInitVal = expNode.optFSPoint;
            }*/
                

            DoubleMatrix1D optDOFs = ccdMin.minimize();
            if(optDOFs==null)//for ellipse bounds, if the highest ellipses don't all intersect together,
                //we can exclude the conformations involving them (set bound to inf)
                return Double.POSITIVE_INFINITY;
            
            double LSBE = of.getValue( optDOFs );
            
            
            if(es.useSVE){
                //cof minimized m...revert to unminimized state
                m.updateCoordinates();
                if(doPerturbations)
                    m.revertPertParamsToCurState();
            }
            
            
            //store initial values for next time
            //TEMPORARILY DE-ACTIVATING OPTFSPOINT
            /*if(useHSer){
                if(expNode.optFSPoint==null)
                    expNode.optFSPoint = optDOFs;
                else
                    expNode.optFSPoint.assign(optDOFs);
            }*/
            
            /*if(LSBE<-0.001){
                System.err.println("NEGATIVE VALUE ENCOUNTERED FOR LSBE: "+LSBE);
                System.err.println("Outputting LSBObjFunction to LSBOF.dat");
                KSParser.outputObject(of, "LSBOF.dat");
                System.err.println("optDOFs: "+optDOFs);
                System.exit(0);
            }*/
            
            numFS++;

            return LSBE;
        }

        
        
        CETObjFunction getNodeObjFunc(QueueNode expNode) {
            //not for splitBySlack
            
            int dLevel = expNode.level;
            int[] conf = expNode.confSoFar;
            
            int AANums[] = new int[dLevel+1];
            int rots[] = new int[dLevel+1];
            for(int level=0; level<=dLevel; level++){

                if(conf[level]>=0){
                    int redRot = pairwiseMinEnergyMatrix.resOffsets[level]+conf[level];
                    //long-form reduced index for rotamer at this level

                    AANums[level] = pairwiseMinEnergyMatrix.indicesEMatrixAA[redRot];
                    rots[level] = pairwiseMinEnergyMatrix.indicesEMatrixRot[redRot];
                }
                else{//special for not-fully-assigned levels
                    AANums[level] = -1;
                    rots[level] = conf[level];
                }
            
            }
            
            //DEBUG!!!!  To use when checking minimization for a particular conformation
            /*
            AANums = new int[] {0, 1, 1, 20, 4, 1, 4};
            rots = new int[] {0, 1, 2, 0, 3, 0, 3};
            */
            
            
            ContSCObjFunction cof = null;//to set DOFs for SVE
                        
            if(es.useSVE){//will need cof to set DOFs for sve
                //set AA types and rotamers up to now
                //as in RotamerSearch
                int curAANum[] = new int[m.numberOfResidues];//current AA number for each residue in the molecule, for cof
                applyRotamers(dLevel, AANums, rots, curAANum);
                cof = new ContSCObjFunction(m,m.numberOfStrands,null,strandRot,curAANum,false,getTransRotStrands(dLevel));
            }
            
            
            
            return CETM.getObjFunc(AANums,rots,false,false,cof);//This can include h bounds if they're set up
        }
        
        
        boolean[] getTransRotStrands(int dLevel){
            //when minimizing a partial conformation,
            //figure out which strands may be allowed to rotate and translate
            boolean transRotStrands[] = new boolean[m.numberOfStrands];
            
            //if a strand can rotate and translate and has either assigned or template residues,
            //let it rotate and translate
            
            for(int str=0; str<m.numberOfStrands; str++){
                if(m.strand[str].rotTrans){
                    
                    if( m.strand[str].numberOfResidues > strandMut[str].length )//strand has template residues
                        transRotStrands[str] = true;
                    
                    for(int i=0; i<=dLevel; i++){
                        if(str==mutRes2Strand[i])//strand has an assigned residue
                            transRotStrands[str] = true;
                    }
                    
                }
            }
            
            return transRotStrands;
        }
        
        

        
        
        void applyRotamers(int dLevel, int[] AANums, int[] rots, int[] curAANum){
            //apply AA types and rotamers up to the current level
            
            for (int i=0; i<=dLevel; i++){
                int str = mutRes2Strand[i];
                int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
                
                //AA types
                if(m.strand[str].isProtein){
                    String curAAName = m.strand[str].residue[strResNum].name;
                    String newAAName = strandRot[str].rl.getAAName(AANums[i]);
                    if(!curAAName.equalsIgnoreCase(newAAName))
                        strandRot[str].changeResidueType(m,strResNum,newAAName,true);
                }
                
                //fill in curAANum while we're at it
                curAANum[m.strand[str].residue[strResNum].moleculeResidueNumber] = AANums[i];
              
                //rotamers
                
                if(doPerturbations){
                    //For DEEPer curRot is the the current RC
                    boolean validRC = ((StrandRCs)strandRot[str]).applyRC(m, strResNum, rots[i]);
                    if(!validRC)
                        throw new RuntimeException("Error: invalid RC " + rots[i] + " at residue " + strResNum +
                                " of strand " + str );
                }
                else if (strandRot[str].rl.getNumRotForAAtype(AANums[i])!=0)//not GLY or ALA
                    strandRot[str].applyRotamer(m, strResNum, rots[i]);
                
                //for gly or ala don't need to do anything
                
                //also make assigned residues flexible
                m.strand[str].residue[strResNum].flexible = true;
            }
          
            //make the other residues flexible
            for(int i=dLevel+1; i<numTreeLevels; i++){
                int str = mutRes2Strand[i];
                int strResNum = strandMut[str][mutRes2StrandMutIndex[i]];
                m.strand[str].residue[strResNum].flexible = false;
            }
            
        }
        
        

	
	//Called by gCompute(.)
	private double gSumMinPVE(int topLevel, int startLevel, int index1, int conf[]){
		
		int index2;
		double sum = 0.0f;
		
		for (int level=startLevel; level<=topLevel; level++){
			
			index2 = nodeIndexOffset[level] + conf[level];	//s at j
			sum += pairwiseMinEnergyMatrix.getPairwiseE( index1, index2 );
		}
		
		return sum;
	}






        private boolean isPruned( int newIndex, int conf[], int topLevel ){
            //Checks if the RC denoted by newIndex if incompatible with the conformation in conf, up to topLevel
            //(as determined by pruned pairs or triples)
            //newIndex should be at a level above topLevel

            if ( splitFlags != null ){

                //See if (newIndex, index in conf) is a pruned pair for each index in conf
                for(int a=0; a<=topLevel; a++){
                    int index2 = nodeIndexOffset[a] + conf[a];
                    if(splitFlags[newIndex][index2])
                        return true;

                    if ( tripleFlags != null){//If split flags is null, the triple flags are also expected to be null, so they will not be checked either
                        for(int b=0; b<a; b++){
                            int index3 = nodeIndexOffset[b] + conf[b];
                            if(tripleFlags[newIndex][index2][index3])
                                return true;
                        }
                    }
                }
            }

            return false;

        }


        private boolean isPruned ( int newIndex, int oldIndex, int conf[], int topLevel ){
            //Checks if the RC denoted by newIndex if incompatible with the current conformation up to topLevel
            //plus oldIndex (as determined by pruned pairs or triples)
            //newIndex should be a higher level than oldIndex, which should be above topLevel


            if ( splitFlags != null ){

                //First check newIndex against conf
                for(int a=0; a<=topLevel; a++){
                    int index2 = nodeIndexOffset[a] + conf[a];
                    if(splitFlags[newIndex][index2])
                        return true;

                    if ( tripleFlags != null){//If split flags is null, the triple flags are also expected to be null, so they will not be checked either
                        for(int b=0; b<a; b++){
                            int index3 = nodeIndexOffset[b] + conf[b];
                            if(tripleFlags[newIndex][index2][index3])
                                return true;
                        }
                    }
                }

                //Check newIndex against oldIndex
                if( splitFlags[newIndex][oldIndex] )
                    return true;

                //Check newIndex against triples involving oldIndex and conf
                if ( tripleFlags != null){
                    for(int b=0; b<=topLevel; b++){
                        int index3 = nodeIndexOffset[b] + conf[b];
                        if(tripleFlags[newIndex][oldIndex][index3])
                            return true;
                    }
                }

            }

            return false;
        }

        public int getNumNodes(){
            return curExpansion.getNumNodes();
        }
         
        
        
}