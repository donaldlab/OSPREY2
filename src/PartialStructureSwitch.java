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
//	PartialStructureSwitch.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import java.util.*;
import java.io.*;


//This perturbation puts part of a different experimental structure into the structure being used,
//e.g. a loop in a different conformation or a ligand in a different pose
//It is like a full structure switch except only the affected residues are switched
//Consequently, the alternate structure need not have the same atoms as the original structure in other regions
//A tripeptide closure on the three residues with lowest B factor is used to fix chain discontinuities
public class PartialStructureSwitch extends Perturbation {

    //lowest-RMSD loop closure

    
    ArrayList< double[]> params = new ArrayList<double[]>();//Perturbation parameters for predecessors and current perturbation
    //Indices: 1. Which set of parameters this is ("parameter state" index) 2. Which perturbation (predecessors in order, then current)


    int numProtRes;//How many protein residues are affected
    boolean nonProtLig;//Is there a non-protein ligand?

    int numStructs;//Number of structures

    String structInfo[];//Information on structures

    ChainSection sections[][];//Information on contiguous protein chain sections of alternate structures
    //for use in placement.  Indices: 1. Structure #  2. Section #


    String PDBList[];

    RotMatrix rm = new RotMatrix();



    //Stored coordinates for residues in specific parameter states (already placed):

    ArrayList< ArrayList < HashMap < String, double[]>>> flexBB = new ArrayList<ArrayList<HashMap<String,double[]>>>();
    //Alternate backbones for affected protein residues, (first list: Parameter state indices; second list: residue (# among affected residues))

    ArrayList< double[]> ligCoord = new ArrayList< double[]>();//Ligand coordinates.  First index: parameter state.  Second index: Coordinates

    
    int ligResNum = -1;//Molecule residue number for the ligand; -1 if the ligand is not affected by this perturbation or if there is none
    //We only need this is the ligand is not protein

    boolean useTC = false;//Use a tripeptide closure for two-anchor sections

    double anchorDiffTol = 0.2f;//Maximum distance allowed between transformed coordinates based on 2 different anchors
    //in a two-anchor section switched without a tripeptide closure


    PartialStructureSwitch(Molecule m, int resList[], BufferedReader br){



        resDirectlyAffected = resList;
        type = "PARTIAL STRUCTURE SWITCH";
        this.m = m;

        if( m.strand[ m.residue[ resDirectlyAffected [resDirectlyAffected.length - 1] ].strandNumber ].isProtein ){
            numProtRes = resDirectlyAffected.length;
            nonProtLig = false;
        }
        else{
            numProtRes = resDirectlyAffected.length - 1;
            nonProtLig = true;
            ligResNum = resDirectlyAffected[resDirectlyAffected.length-1];
        }


        //Break the affected protein residues up into continuous sections of protein chain
        ArrayList<Integer> sectionStarts = new ArrayList<Integer>();//Where in affectedRes the sections start
        ArrayList<Integer> sectionSizes = new ArrayList<Integer>();//How many residues they contain

        sectionStarts.add(0);
        int numSections = 1;
        int sizeSoFar = 1;
        Residue lastRes = m.residue[resDirectlyAffected[0]];

        for( int a=1; a<resDirectlyAffected.length; a++ ){

            Residue thisRes = m.residue[resDirectlyAffected[a]];

            if( ( thisRes.strandNumber == lastRes.strandNumber )
                    && ( thisRes.strandResidueNumber == lastRes.strandResidueNumber + 1)
                    && ( m.checkNBonded(resDirectlyAffected[a]) ) ){//If thisRes and lastRes are connected they go in the same section

                sizeSoFar++;
            }
            else{
                sectionSizes.add( sizeSoFar );
                sectionStarts.add( a );
                sizeSoFar = 1;
                numSections++;
            }

            lastRes = thisRes;
        }

        sectionSizes.add(sizeSoFar);
        

        //Now read in the alternate structures
        try{

            /*The extra information for this perturbation is as follows:
             - Number of structures (1 line; includes original): [number] structures
             - For each structure, the file name, and then the PDB residue numbers affected in that structure
             (including for the ligand, since it might be different: i.e. it shouldn't just say "LIG")


*
             Example entry for this perturbation:

            PARTIAL STRUCTURE SWITCH
            5 6 7
            3 structures
            1ABC.pdb 6 7 8
            2DEF.pdb 8 9 10
            3 states
            0 0
            1 1
            2 2

            This indicates that the loop being replaced is residues 5-7 in the original structure
            but residues 6-8 in 1ABC and 8-10 in 2DEF.

            */

            String line = br.readLine();

            StringTokenizer st = new StringTokenizer(line," ");
            if( st.countTokens() != 2 ){
                System.err.println("Error reading partial structure switch number of structures: " + line);
                System.exit(1);
            }

            numStructs = Integer.valueOf(st.nextToken());
            structInfo = new String[numStructs-1];


            sections = new ChainSection[numStructs][numSections];
            PDBList = new String[numStructs];
            PDBList[0] = "default structure";

            //Store the affected residues & anchor atoms from the different structures (store all coordinates for the non-protein residues; backbones for flexible residues)

            //We reserve structure # 0 for the unperturbed state
            for( int structIndex=1; structIndex<numStructs; structIndex++){

                structInfo[structIndex-1] = br.readLine();
                
                st = new StringTokenizer(structInfo[structIndex-1]," ");
                if( st.countTokens() != resDirectlyAffected.length+1 ){
                    System.err.println("Error reading partial structure switch structure information: " + structInfo[structIndex-1]);
                    System.exit(1);
                }

                String PDBFile = st.nextToken();


                Molecule m2 = new Molecule();

                FileInputStream is = new FileInputStream( PDBFile );
                new PDBChemModel(m2,is);

                PDBList[structIndex] = PDBFile;

                int resAffected2[] = new int[resDirectlyAffected.length];//Residues affected in the alternate structure

                for(int a=0; a<resDirectlyAffected.length; a++)
                    resAffected2[a] = m2.mapPDBresNumToMolResNum( Integer.valueOf(st.nextToken()) );

                for( int sec=0; sec<numSections; sec++ ){

                    ChainSection cs = new ChainSection();

                    cs.startRes = sectionStarts.get(sec);
                    cs.numRes = sectionSizes.get(sec);

                    Residue beginRes = m2.residue[ resAffected2[cs.startRes] ];
                    Residue endRes = m2.residue[ resAffected2[cs.startRes + cs.numRes - 1] ];
                    boolean nonCovalentAnchor = false;//The section is not covalently bound to other residues so find a non-covalently-bound anchor
                    cs.anchor = new Atom[3];
                    cs.twoAnchors = false;

                    if( m.strand[beginRes.strandNumber].isProtein ){//isProtein flags are not set for m2 but the same strands should be protein as in m

                        cs.isProtein = true;

                        for(int a=cs.startRes; a<cs.startRes+cs.numRes ;a++)//Store backbone information
                            cs.backbones.add( storeResBB( m2.residue[resAffected2[a]], m2 ) );

                        boolean beginResNTerm = ! m2.checkNBonded(resAffected2[cs.startRes]);//beginRes is at the N terminus
                        boolean endResCTerm = ! m2.checkCBonded(resAffected2[cs.startRes + cs.numRes - 1]);//endRes is at the C terminus

                        //Now get anchors
                        //Note that the Atom.coord arrays and m2.actualCoordinate arrays are the same for m2
                        //So we can use atoms directly from m2 as anchor atoms
                        if( ! beginResNTerm ){
                            cs.anchor[0] = beginRes.getAtomByName("N");
                            Residue prevRes = m2.residue[ resAffected2[cs.startRes] - 1 ];
                            cs.anchor[1] = prevRes.getAtomByName("C");
                            cs.anchor[2] = prevRes.getAtomByName("CA");
                        }

                        if( ! endResCTerm ){
                            cs.anchor2 = new Atom[3];
                            cs.anchor2[0] = endRes.getAtomByName("C");
                            Residue nextRes = m2.residue[ resAffected2[cs.startRes + cs.numRes - 1] + 1 ];
                            cs.anchor2[1] = nextRes.getAtomByName("N");
                            cs.anchor2[2] = nextRes.getAtomByName("CA");
                        }


                        if( beginResNTerm && endResCTerm ){
                            nonCovalentAnchor = true;
                        }
                        else if ( beginResNTerm && ( ! endResCTerm )  ){
                            cs.anchor = cs.anchor2;
                        }
                        else if ( ! ( beginResNTerm || endResCTerm ) ){//Need both anchors

                            cs.twoAnchors = true;

                            if( cs.numRes < 3 && useTC ){//We can't do a tripeptide closure with less than 3 residues
                                System.err.println( "Tripeptide closure not used for partial structure switch: too few residues in section "
                                        + sec + " of " + PDBList[structIndex] );
                                useTC = false;
                            }

                            if (useTC) {
                                
                                //We will need a tripeptide closure
                                cs.tcRes = new int[3];//These will be the 3 residues with the lowest CA B-factors in the section
                                double B[] = new double[3];


                                for(int a=0;a<3;a++){
                                    cs.tcRes[a] = cs.startRes+a;
                                    B[a] = m2.residue[resAffected2[cs.startRes+a]].getAtomByName("CA").BFactor;
                                }


                                for(int a=cs.startRes+3; a<cs.startRes + cs.numRes; a++){

                                    int worst = 0;//Which of the three residues in tcRes (0, 1, or 2) has the worst (biggest) B-factor
                                    if( ( B[2] >= B[1] ) && ( B[2] >= B[0] ) )
                                        worst = 2;
                                    else if( ( B[1] >= B[2] ) && ( B[1] >= B[0] ) )
                                        worst = 1;

                                    double newB = m2.residue[resAffected2[a]].getAtomByName("CA").BFactor;

                                    if( newB < B[worst] ){
                                        B[worst] = newB;
                                        cs.tcRes[worst] = a;
                                    }
                                }


                                Arrays.sort(cs.tcRes);//tcRes needs to be in ascending order

                                int tcMolResNum[] = new int[3];
                                for(int a=0;a<3;a++)
                                    tcMolResNum[a] = resDirectlyAffected[cs.tcRes[a]];//Molecule residue number in the original molecule

                                cs.tc = new TripeptideClosure( m2, tcMolResNum );
                            }
                        }

                    }
                    else{
                        cs.isProtein = false;
                        cs.twoAnchors = false;//This will be the ligand and won't require two anchors

                        int numAtoms=0;

                        for(int a=cs.startRes; a<cs.startRes+cs.numRes ;a++)
                            numAtoms += m2.residue[resAffected2[a]].numberOfAtoms;

                        cs.coord = new double[3*numAtoms];
                        int place=0;

                        for(int a=cs.startRes; a<cs.startRes+cs.numRes ;a++){//Store coordinates
                            int firstAtNum = m2.residue[resAffected2[a]].atom[0].moleculeAtomNumber;
                            int coordCount = 3*m2.residue[resAffected2[a]].numberOfAtoms;
                            System.arraycopy(m2.actualCoordinates, 3*firstAtNum, cs.coord, place, coordCount);
                            place += coordCount;
                        }


                        nonCovalentAnchor = true;
                    }



                    if( nonCovalentAnchor ){

                        double centroid[] = new double[3];
                        int numAtoms = 0;
                        for( int a=cs.startRes; a<cs.startRes+cs.numRes; a++ ){
                            Residue res = m2.residue[resAffected2[a]];

                            for( Atom at : res.atom ){
                                centroid = rm.add(centroid, at.coord);
                                numAtoms++;
                            }
                        }

                        centroid = rm.scale(centroid, 1/(double)numAtoms);

                        //Create an anchor consisting of the three closest CAs to this section's centroid unaffected by this perturbation
                        double dist[] = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};

                        for( Strand str : m2.strand ){

                            if( m.strand[str.number].isProtein ){//Again we assume the same residues are protein in m2 and m

                                for( Residue res : str.residue ){

                                    int worst = 0;//Which of the three distances in dist (0, 1, or 2) is largest (i.e. worst)
                                    if( ( dist[2] >= dist[1] ) && ( dist[2] >= dist[0] ) )
                                        worst = 2;
                                    else if( ( dist[1] >= dist[2] ) && ( dist[1] >= dist[0] ) )
                                        worst = 1;


                                    Atom CA = res.getAtomByName("CA");
                                    double checkDist = rm.norm( rm.subtract( CA.coord, centroid ) );

                                    if( checkDist < dist[worst] ){
                                        cs.anchor[worst] = CA;
                                        dist[worst] = checkDist;
                                    }

                                }


                            }

                        }

                    }

                    sections[structIndex][sec] = cs;
                }

            }
        }
        catch(Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
    }



    public boolean doPerturbationMotion(double param){

        if( param < 0 || param >= sections.length ){
            System.err.println("Bad partial structure switch parameter: " + param + ".  Ignoring perturbation state.");
            return false;
        }

        
        int state = -1;//If the predecessors & this perturbation have been in this state before this will indicate which set of known parameters this is

        for( int a=0; ( a<params.size() ) && ( state == -1 ); a++ ){

            state = a;

            for( int b=0; b<predecessors.length; b++ ){
                if( params.get(state)[b] != m.perts[predecessors[b]].curParam ){//Predecessor parameters
                    state = -1;
                    break;
                }
            }

            if(state != -1){
                if( params.get(state)[predecessors.length] != param )
                    state = -1;
            }

        }

        if(state == -1){//New state
            boolean outcome = placeRes(param);
            if(!outcome)//placement failed
                return false;
            
            state = params.size() - 1;
        }


        for(int a=0; a<numProtRes; a++)
            restoreResBB( m.residue[resDirectlyAffected[a]], flexBB.get(state).get(a) );

        if( nonProtLig )
            System.arraycopy(ligCoord.get(state), 0, m.actualCoordinates, 3*m.residue[ligResNum].atom[0].moleculeAtomNumber, ligCoord.get(state).length );


        return true;
    }



    public boolean placeRes(double param){
    //Places the affected residues from the alternate structure indicated by param
        //in the molecule m and stores their coordinates in flexBB, ligCoords
        //If a residue has other stuff in the strand it is placed based on that
        //Any others are aligned using their closest three atoms
        //Either way these are the "anchor" atoms

        int structIndex = (int)param;

        ArrayList < HashMap < String, double[]>> stateBB = new ArrayList < HashMap < String, double[]>>();
        //The array list of hashmaps representing the backbone coordinates for all residues in this parameter state

        for(int sec=0; sec<sections[structIndex].length; sec++){
            ChainSection cs = sections[structIndex][sec];

            double[][] anchorStarts = new double[3][3];//Anchor coordinates from alternate structure
            double[][] anchorTargets = new double[3][3];//Anchor coordinates from cuurent structure
            
            for(int a=0;a<3;a++){
                anchorStarts[a] = cs.anchor[a].coord.clone();
                anchorTargets[a] = m.getActualCoord( m.residue[cs.anchor[a].moleculeResidueNumber].getAtomByName(cs.anchor[a].name).moleculeAtomNumber );
            }

            double matrix1[][] = rm.getSuperposingRotMatrix( rm.subtract( anchorStarts[1], anchorStarts[0] ),
                    rm.subtract( anchorTargets[1], anchorTargets[0] ),
                    rm.subtract( anchorStarts[2], anchorStarts[0] ),
                    rm.subtract( anchorTargets[2], anchorTargets[0] ) );
            //Matrix based on first anchor
            //The transformation consists of subtracting anchorStarts[0], applying matrix1, and then adding anchorTargets[0]


            if(cs.twoAnchors){//Need to handle the second anchor and the tripeptide closure
                //As currently implemented, this will only happen in protein strands, so we only handle the protein case


                //Make the rotation and translation for the second anchor

                double[][] anchor2Starts = new double[3][3];//Anchor coordinates from alternate structure
                double[][] anchor2Targets = new double[3][3];//Anchor coordinates from cuurent structure

                for(int a=0;a<3;a++){
                    anchor2Starts[a] = cs.anchor2[a].coord.clone();
                    anchor2Targets[a] = m.getActualCoord( m.residue[cs.anchor2[a].moleculeResidueNumber].getAtomByName(cs.anchor2[a].name).moleculeAtomNumber );
                }

                double matrix2[][] = rm.getSuperposingRotMatrix( rm.subtract( anchor2Starts[1], anchor2Starts[0] ),
                    rm.subtract( anchor2Targets[1], anchor2Targets[0] ),
                    rm.subtract( anchor2Starts[2], anchor2Starts[0] ),
                    rm.subtract( anchor2Targets[2], anchor2Targets[0] ) );




                for(int a=0; a<cs.numRes; a++){
                    
                    //Some of these atoms will be moved again by the tripeptide closure if useTC == true

                    HashMap<String, double[]> resBB = new HashMap<String, double[]>();

                    double newCoord[];


                    for( String s : cs.backbones.get(a).keySet() ){

                        if(useTC){

                            if( a + cs.startRes < cs.tcRes[1] )//Translate and rotate residues before the middle of the closure using the first anchor
                                newCoord = transformEntry( cs.backbones.get(a), s, matrix1, anchorStarts[0], anchorTargets[0] );

                            else if ( a + cs.startRes > cs.tcRes[1] )//Use the second anchor after the middle of the closure
                                newCoord = transformEntry( cs.backbones.get(a), s, matrix2, anchor2Starts[0], anchor2Targets[0] );

                            else if ( s.equalsIgnoreCase("C") || s.equalsIgnoreCase("O") )//Split the middle residue according to which CA its pieces will need to rotate around during tripeptide closure
                                newCoord = transformEntry( cs.backbones.get(a), s, matrix2, anchor2Starts[0], anchor2Targets[0] );
                            else
                                newCoord = transformEntry( cs.backbones.get(a), s, matrix1, anchorStarts[0], anchorTargets[0] );
                        }
                        else{
                            //Find the rotated and translated coordinates for each residue using both anchors
                            //Then average them weighted by the distance to each anchor
                            //If the two sets of transformed coordinates for any atom are farther apart than anchorDiffTol then the perturbation fails

                            double[] newCoord1 = transformEntry( cs.backbones.get(a), s, matrix1, anchorStarts[0], anchorTargets[0] );
                            double[] newCoord2 = transformEntry( cs.backbones.get(a), s, matrix2, anchor2Starts[0], anchor2Targets[0] );

                            if ( rm.norm( rm.subtract(newCoord1, newCoord2)) > anchorDiffTol ){
                                System.err.println("Error: Partial structure switch to " + PDBList[structIndex]
                                + " failed, section " + sec + ".  Ignoring perturbation state." );
                                return false;
                            }

                            double dist1 = rm.norm( rm.subtract(newCoord1, anchorTargets[0]) );
                            double dist2 = rm.norm( rm.subtract(newCoord2, anchor2Targets[0]) );

                            double wt1 = dist2/(dist1+dist2);//Weights for average (wt1/wt2 = dist2/dist1; wt1+wt2=1)
                            double wt2 = dist1/(dist1+dist2);

                            newCoord = rm.add( rm.scale(newCoord1, wt1) , rm.scale(newCoord2, wt2) );
                        }

                        resBB.put(s, newCoord );
                    }

                    stateBB.add(resBB);
                }



                if(useTC){

                    //Calculate the tripeptide closure
                    double firstN[] = stateBB.get( cs.tcRes[0] ).get("N");//Coordinates of the first closure residue's N (as transformed in the loop that just finished)
                    double firstCA[] = stateBB.get( cs.tcRes[0] ).get("CA");//Its CA
                    double firstC[] = stateBB.get( cs.tcRes[0] ).get("C");

                    double midN[] = stateBB.get( cs.tcRes[1] ).get("N");//Starting coordinates of the middle N (transformed as above but tripeptide closure not yet performed)
                    double midCA[] = stateBB.get( cs.tcRes[1] ).get("CA");
                    double midC[] = stateBB.get( cs.tcRes[1] ).get("C");

                    double lastN[] = stateBB.get( cs.tcRes[2] ).get("N");
                    double lastCA[] = stateBB.get( cs.tcRes[2] ).get("CA");
                    double lastC[] = stateBB.get( cs.tcRes[2] ).get("C");


                    double r_soln_n[][][] = new double[16][3][3];
                    double r_soln_a[][][] = new double[16][3][3];
                    double r_soln_c[][][] = new double[16][3][3];

                    int numSoln = cs.tc.solve_3pep_poly( firstN, firstCA, lastCA,
                            lastC, r_soln_n, r_soln_a, r_soln_c );

                    if( numSoln == 0 ){//Tripeptide closure failed
                        System.err.println("Error: Tripeptide closure for partial structure switch to " + PDBList[structIndex]
                                + " failed, section " + sec + ".  Ignoring perturbation state.");
                        return false;
                    }
                    else{

                        //Find the tripeptide closure solution that changes the conformation of the alternate structure the least
                        //This will be done by a least-squares comparison to the coordinates before closure

                        int bestSoln = -1;
                        double bestSum = Double.POSITIVE_INFINITY;

                        for(int soln=0; soln<numSoln; soln++){

                            double sum = rm.normsq( rm.subtract( firstC, r_soln_c[soln][0] ) );
                            sum += rm.normsq( rm.subtract( midN, r_soln_n[soln][1] ) );
                            sum += rm.normsq( rm.subtract( midCA, r_soln_a[soln][1] ) );
                            sum += rm.normsq( rm.subtract( midC, r_soln_c[soln][1] ) );
                            sum += rm.normsq( rm.subtract( lastN, r_soln_n[soln][2] ) );

                            if(sum < bestSum)
                                bestSoln = soln;
                        }


                        //Now rotate and translate everything else into place using the best solution


                        double newFirstCA[] = r_soln_a[bestSoln][0];//New position of the first CA of the tripeptide closure; used with matrix3 to transform the first part of the closure
                        double newLastCA[] = r_soln_a[bestSoln][2];//New position of the third CA (used with matrix4 to transform the second part)


                        double matrix3[][] = rm.getSuperposingRotMatrix( rm.subtract( midCA, firstCA ), rm.subtract( r_soln_a[bestSoln][1], firstCA),
                            rm.subtract( midN, firstCA ), rm.subtract( r_soln_n[bestSoln][1], firstCA) );

                        double matrix4[][] = rm.getSuperposingRotMatrix( rm.subtract( midC, lastCA ), rm.subtract( r_soln_c[bestSoln][1], lastCA),
                            rm.subtract( lastN, lastCA ), rm.subtract( r_soln_n[bestSoln][2], lastCA) );



                        for(int b=cs.tcRes[0]; b<=cs.tcRes[2]; b++){
                            //Loop through each residue within the tripeptide closure region
                            //b is an index in resDirectlyAffected
                            HashMap<String, double[]> resBB = stateBB.get(b);

                            double newCoord[];

                            for( String s : resBB.keySet() ){


                                int whichTransformation = 0;
                                //0 indicates no transformation; 3 indicates using matrix3 and firstCA --> newFirstCA;
                                //4 indicates using matrix4 and lastCA --> newLastCA

                                if( b > cs.tcRes[0] && b < cs.tcRes[1] )
                                    whichTransformation = 3;
                                else if ( b > cs.tcRes[1] && b < cs.tcRes[2] )
                                    whichTransformation = 4;
                                else if ( b == cs.tcRes[1] ){
                                    if( s.equalsIgnoreCase("C") || s.equalsIgnoreCase("O") )
                                        whichTransformation = 4;
                                    else
                                        whichTransformation = 3;
                                }
                                else if(b == cs.tcRes[0] && (s.equalsIgnoreCase("C") || s.equalsIgnoreCase("O") ) )
                                     whichTransformation = 3;
                                else if(b == cs.tcRes[2] && (s.equalsIgnoreCase("N") || s.equalsIgnoreCase("H") ) )
                                    whichTransformation = 4;
                                //Otherwise there is no transformation


                                if(whichTransformation == 3){
                                    newCoord = transformEntry( resBB, s, matrix3, firstCA, newFirstCA );
                                    resBB.put(s, newCoord);
                                }
                                else if(whichTransformation == 4){
                                    newCoord = transformEntry( resBB, s, matrix4, lastCA, newLastCA );
                                    resBB.put(s, newCoord);
                                }
                            }
                        }
                    }
                }

            }
            else{//Just put everything in place using the first anchor

                if(cs.isProtein){
                    for(int a=0;a<cs.backbones.size();a++){//Translate and rotate each residue
                        HashMap<String, double[]> resBB = new HashMap<String, double[]>();
                        for( String s : cs.backbones.get(a).keySet() ){
                            double newCoord[] = transformEntry( cs.backbones.get(a), s, matrix1, anchorStarts[0], anchorTargets[0] );
                            resBB.put(s, newCoord );
                        }
                        stateBB.add(resBB);
                    }
                }
                else{
                    double stateLigCoord[] = new double[cs.coord.length];
                    double x[] = new double[3];
                    for(int a=0; a<cs.coord.length; a+=3){//Translate and rotate each atom
                        System.arraycopy(cs.coord, a, x, 0, 3);
                        x = transform ( x, matrix1, anchorStarts[0], anchorTargets[0] );
                        System.arraycopy(x, 0, stateLigCoord, a, 3);
                    }
                    ligCoord.add(stateLigCoord);
                }

            }
        }


        //Store the coordinates and parameters

        int state = params.size();//The state number this information will be stored as

        //Store parameters
        params.add(new double[predecessors.length+1]);

        for(int b=0;b<predecessors.length;b++)
            params.get(state)[b] = m.perts[predecessors[b]].curParam;

        params.get(state)[predecessors.length] = param;

        //Store backbones
        flexBB.add(stateBB);

        //Non-protein ligand coordinates have already been stored if relevant

        return true;

    }



    private double[] transform( double[] x, double[][] mtx, double[] anchorStart, double[] anchorTarget ){
            //This transformation is used heavily by placeRes
            //We subtract anchorStart from x, rotate it using mtx, and then add anchorTarget
        return rm.add( rm.applyRotMatrix( mtx, rm.subtract(x, anchorStart ) ), anchorTarget );
    }


    private double[] transformEntry( HashMap<String, double[]> hm, String s, double[][] mtx, double[] anchorStart, double[] anchorTarget ){
        //Transforms an entry from a backbone coordinate hashmap
        //Works like transform for a set of coordinates, but just rotates a vector (i.e. the "CACB" entry)

        if( s.equalsIgnoreCase("CACB") )
            return rm.applyRotMatrix( mtx, hm.get(s) );
        else
            return rm.add( rm.applyRotMatrix( mtx, rm.subtract( hm.get(s), anchorStart ) ), anchorTarget );
    }


    public void setDefaultParams(){

        minParams = new double[numStructs];
        maxParams = new double[numStructs];

        for(int a=0;a<numStructs;a++)
            minParams[a] = maxParams[a] = a;
    }



    private class ChainSection {//A section of chain from an alternate structure, to be placed into the main structure m as a unit

        boolean isProtein;
        boolean twoAnchors;

        int numRes;//Number of residues in this section
        int startRes;//Which residue (given as an index in affectedRes) this sections starts at

        ArrayList< HashMap < String, double[]>> backbones = new ArrayList< HashMap <String, double[]>>();//Backbone coordinates, taken directly from the alternate structure

        double coord[] = null;//If this is not protein, all the coordinates are stored here.  Currently only the ligand can be handled this way.


        Atom anchor[];//We can place the section by aligning these three atoms with their equivalents in their target structures

        //If the section is covalently bound to other residues on both sides, it needs this second anchor (also three atoms)
        Atom anchor2[] = null;

        //The residues in tcRes (given as indices in affectedRes) are involved in a tripeptide closure, tc, to close the gap
        int tcRes[] = null;
        TripeptideClosure tc = null;

        //Atoms before the closed sections are placed using the first anchor; those after, using the second
    }


    @Override
    public void writeExtraInfo( BufferedWriter bw ){

        try{
            bw.append(numStructs + " structures");
            bw.newLine();
            for(int a=0;a<numStructs-1;a++){
                bw.append(structInfo[a]);
                bw.newLine();
            }
        }
        catch(IOException e){
            System.err.println(e.getMessage());
        }
    }


    
    @Override
    public boolean isParamAngle(){
        return false;
    }

}
