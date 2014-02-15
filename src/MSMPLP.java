import java.util.ArrayList;
import java.util.Collections;


public class MSMPLP {
	//number of residues under consideration
	private int numResidues;

	//number of rotamers possible for each residue (given by the assigned AA type)
	private int numRotForRes[] = null;

	//the offset in the array index for each level
	private int nodeIndexOffset[] = null;
	
	//the reduced min pairwise energy matrix
	private double [][][][] unifiedMinEnergyMatrix = null;
	
	MSMPLP(int aNumResidues, int aNumRotForRes[], ReducedEnergyMatrix aPairwiseMinEnergyMatrix){
		numResidues = aNumResidues;

		numRotForRes = new int [numResidues];
		nodeIndexOffset = new int [numResidues];
		int numTotalNodes = 0;


		for (int i=0; i<numResidues; i++){
			nodeIndexOffset[i] = numTotalNodes;
			numRotForRes[i] = aNumRotForRes[i];
			numTotalNodes += aNumRotForRes[i];
		}
		unifiedMinEnergyMatrix = MSMPLP.mergeIntraAndPairMats(numResidues, numRotForRes, nodeIndexOffset, aPairwiseMinEnergyMatrix);
	}
	
	// Computes the low-energy bound using the EMPLP algorithm
	// availableRots is a "numRes*rots" matrix that tells us which 
	// rotamers at each position are available.  Only those rotamers
	// will be considered for the calculation.
	public double optimizeEMPLP(int partialConf[], int iterations){

		
		// availableRots is just a list of lists with the list of rotamers available for the calculation at each rotamer; it is more convenient
		// 		than partialConf.  This code might be unnecessary but it makes the algorithm more elegant, IMO
		int availableRots[][] = new int[numResidues][];
		for(int res = 0; res < numResidues; res++){			
			if(partialConf[res]== -1){ 
				availableRots[res] = new int[numRotForRes[res]];
				for(int rot =0; rot < numRotForRes[res]; rot++){
					availableRots[res][rot] = rot;
				}
			}
			else{ // residue has a rotamer already assigned in the partialConf.
				availableRots[res] = new int[1];
				availableRots[res][0] = partialConf[res];
			}
		}
		
		// Here comes the actual algorithm
		
		// lambda is the message matrix for EMPLP; we set all messages to zero
		double lambda[][][] = CreateMatrix.create3DMsgMat(numResidues, numRotForRes, 0.0f);
		
		// The belief on each rotamer.
		double belief[][] = CreateMatrix.create2DRotMatrix(numResidues, numRotForRes, 0.0f);
			
		
		
	    // EMPLP algorithm.  
	    // The algorithm should loop until convergence.  For now we do it I times. 
	    // Complexity of EMPLP: O(I*numRes*numRes*rotsPerRes*rotsPerRes*numRes)
		for(int i = 0; i < iterations; i++){
			for (int resI = 0; resI < numResidues; resI++){
				for(int resJ = resI+1; resJ < numResidues; resJ++){
					// We first update lambda[resJ][resI][rotIR] and immediately after we update lambda[resI][resJ][rotJS]
					//for(int rotIR = 0; rotIR < numRotForRes[resI]; rotIR++){
					for(int rotIR: availableRots[resI]){
						belief[resI][rotIR] -= lambda[resJ][resI][rotIR];
						
						ArrayList <Double> msgsFromRotsAtJ_to_rotIR = new ArrayList<Double>();
						//for (int rotJS = 0; rotJS < numRotForRes[resJ]; rotJS++){
						for(int rotJS: availableRots[resJ]){
							msgsFromRotsAtJ_to_rotIR.add(belief[resJ][rotJS]-lambda[resI][resJ][rotJS] + unifiedMinEnergyMatrix[resI][rotIR][resJ][rotJS]);
						}
						
						lambda[resJ][resI][rotIR] = -0.5f*belief[resI][rotIR] + 0.5f*Collections.min(msgsFromRotsAtJ_to_rotIR);
						belief[resI][rotIR] += lambda[resJ][resI][rotIR];						
					}
					// Now we update lambda[resI][resJ][rotJS]
					//for(int rotJS = 0; rotJS < numRotForRes[resJ]; rotJS++){
					for(int rotJS: availableRots[resJ]){
						belief[resJ][rotJS] -= lambda[resI][resJ][rotJS];
						
						ArrayList <Double> msgsFromRotsAtI_to_rotJS = new ArrayList<Double>();
						//for (int rotIR = 0; rotIR < numRotForRes[resI]; rotIR++){
						for(int rotIR: availableRots[resI]){
							msgsFromRotsAtI_to_rotJS.add(belief[resI][rotIR]-lambda[resJ][resI][rotIR] + unifiedMinEnergyMatrix[resI][rotIR][resJ][rotJS]);
						}
						
						lambda[resI][resJ][rotJS] = -0.5f*belief[resJ][rotJS] + 0.5f*Collections.min(msgsFromRotsAtI_to_rotJS);
						belief[resJ][rotJS] += lambda[resI][resJ][rotJS];
					}
				}				
			}
			
		}
		double Ebound = 0.0f;
		for (int resI = 0; resI < numResidues; resI++){
			Ebound += Util.min(belief[resI]);
		}
		//System.out.println("MPLP energy = "+ Ebound);
		return Ebound;
		
		
		
		
	}
	
	// The original dual derivation of MPLP did not consider intra-energies.  Thus, 
	//	I found that the easiest way to solve this would be to unify the intra and pair energies by 
	//	dividing the intra interaction between n-1 and adding it to the edge..
	// Note that the unified emat is 4D whereas the standard matrices used by A* in OSPREY are 2D.
	static double [][][][] mergeIntraAndPairMats(int numRes, int rotsPerRes [], int nodeIndexOffset[], ReducedEnergyMatrix pairEmat ){
		double unifiedEmat [][][][] = CreateMatrix.create4DRotMatrix(numRes, rotsPerRes, 0.0f);
		for(int resI = 0; resI < numRes; resI++){
			for(int rotIR = 0; rotIR < rotsPerRes[resI]; rotIR++){
				for (int resJ = 0; resJ < numRes; resJ++){
					for (int rotJS = 0; rotJS < rotsPerRes [resJ]; rotJS++){
						unifiedEmat [resI][rotIR][resJ][rotJS]= pairEmat.getPairwiseE(nodeIndexOffset[resI]+rotIR,nodeIndexOffset[resJ]+rotJS);
						unifiedEmat [resI][rotIR][resJ][rotJS]+= pairEmat.getIntraAndShellE(nodeIndexOffset[resI]+rotIR)/(numRes-1);
						unifiedEmat [resI][rotIR][resJ][rotJS]+= pairEmat.getIntraAndShellE(nodeIndexOffset[resJ]+rotJS)/(numRes-1);
						unifiedEmat [resJ][rotJS][resI][rotIR]=  unifiedEmat[resI][rotIR][resJ][rotJS];
					}
				}
			}
		}
		return unifiedEmat;
	}
	
	

}
