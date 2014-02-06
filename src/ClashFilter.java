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
//	ClashFilter.java
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


public class ClashFilter {
    
    CETMatrix cetm;
    int[] seq;
    double cutoff;
    
    boolean allClashing = false;//clash for all sequences

    
    public ClashFilter(CETMatrix lsbm, int[] seq, double cutoff) {
        
        this.cetm = lsbm;
        this.seq = seq;
        this.cutoff = cutoff;
        
        
        //constructor should change seq at any single-AA-type res;
        //this will sub in the alternate ligand
        for(int res=0; res<lsbm.numRes; res++){
            int soleAA = -1;
            boolean options = false;//keep seq as it is: multiple options
            
            for(int AA=0; AA<lsbm.intraAndShellBounds[res].length; AA++){
                if(lsbm.intraAndShellBounds[res][AA]!=null){
                    if(soleAA==-1)
                        soleAA = AA;
                    else{
                        options = true;
                        break;//keep seq as it is
                    }
                }
            }
            
            if(!options){
                //at most 1 AA type!
                if(soleAA==-1)//no AA types available...count as all clashing
                    allClashing = true;
                else
                    seq[res] = soleAA;
            }
        }
    }
    
    
    
    
    boolean hasNonClashing(){
        //is there a non-clashing conf with the given sequence?
        
        if(allClashing)
            return true;
        
        //going to start with EPolys
        //but ellipses can probably speed this up substantially
        
        //we'll do DFS since the first non-clashing conf we find counts
        return nonClashingHelper(-1,new int[cetm.numRes]);
    }
    
    
    boolean nonClashingHelper(int level, int[] conf){
        //given a conformation at the given sequence, defined up to the given level,
        //is there a non-clashing conformation within the definition?
        //conf is RC indices for the given sequence
        
        if(hasClashes(level,conf))
            return false;
        
        if(level==cetm.numRes-1)//fully defined
            return true;
        else {
            for(int rot=0; rot<cetm.intraAndShellBounds[level+1][seq[level+1]].length; rot++){
                if(cetm.intraAndShellBounds[level+1][seq[level+1]][rot] != null){
                    conf[level+1] = rot;
                    if(nonClashingHelper(level+1,conf))
                        return true;
                }
            }
            return false;
        }
    }
    
    
    boolean hasClashes(int level, int[] conf){
        //build LSBObjFunction up to level (like g-score) and minimize
        //later probs add h-score (envelope!)
        //based on MSAStar.getNodeObjFunc

        /*(if(splitBySlack)
            return LSBM.getObjFunc(null,conf,useHSer);//rotamer sets handled by AStarAxe
            //dLevel = numTreeLevels-1;//we list RCs or RC sets at all residues...negative numbers
        //for not fully assigned
        */
        
        //first check if there are any unavoidable clashes involving unassigned residues
        //ok i guess we can also do this later
        /*for(int res=level+1; res<lsbm.numRes; res++){
            for(int res2=0; res2<res; res2++){
                consider;
            }
        }*/
        
        if(level==-1)
            return false;//nothing defined...can't have clashes yet
        
        
        //have any explicitly assigned pairwise clashes appeared with the new assignment?
        for(int res=0; res<level; res++){
            if(cetm.pairwiseBounds[level][seq[level]][conf[level]][res][seq[res]][conf[res]]==null)
                return true;
        }
        
        //ContSCObjFunction of = ((ContSCObjFunction)ccdMin.objFcn);
        int AANums[] = new int[level+1];
        int rots[] = new int[level+1];
        for(int res=0; res<=level; res++){
            AANums[res] = seq[res];
            rots[res] = conf[res];
        }
            
        CETObjFunction lof = cetm.getObjFunc(AANums,rots,false,false,null);
        CCDMinimizer lmin = new CCDMinimizer(lof,false);
        DoubleMatrix1D optDOFs = lmin.minimize();
        
        
        if(optDOFs==null)//no valid conf
            return true;
        
        double bestVal = lof.getValue(optDOFs);

        if(bestVal>cutoff)//too high best energy
            return true;
        else
            return false;
    }
            
            
    
}