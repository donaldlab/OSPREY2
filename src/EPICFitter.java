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
//	EPICFitter.java
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
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;


//this class generates EPIC polynomial fits for a given RC (intra+shell energy) or pair (pairwise energy)
//it can try various orders until it gets a good one
public class EPICFitter {

    EPICSettings es;
    
    //template for the generated ContETerms
    int[] DOFNums;
    int numDOFs;
    DoubleMatrix1D DOFmax;
    DoubleMatrix1D DOFmin;
    DoubleMatrix1D center;
    double minE;
    
    CCDMinimizer ccdMin;
    
    EPoly PCTemplate = null;//template to use for EPolyPCs
    //needs to be set (by quadratic fitting) before trying to make any of these
    
    
    static int sampPerParam = 10;
    
    public EPICFitter( CCDMinimizer ccdMin, EPICSettings eset,
            DoubleMatrix1D cen, double me ){
        //given the CCDMinimizer used to minimize for a rotamer pair (or intra+shell)
        //construct an EpicFitter for it
        
        ContSCObjFunction of = (ContSCObjFunction)ccdMin.objFcn;
        this.ccdMin = ccdMin;
        
        DOFNums = of.DOFNums.clone();
        numDOFs = of.numDOFs;
        DOFmax = ccdMin.DOFmax;
        DOFmin = ccdMin.DOFmin;

        es = eset;
        center = cen;
        minE = me;
    }

    
    
    
    
    EPoly doFit(FitParams fp){

        int numParams = SeriesFitter.getNumParams(numDOFs, false, fp.order);
        int numSamples = sampPerParam * numParams;
        
        EPolyPC PCFit = null;


        if(fp.PCOrder>fp.order){//need to do PC fit

            PCFit = new EPolyPC(PCTemplate,fp.order,fp.PCOrder,fp.PCFac);

            int numPCs = SeriesFitter.countTrue(PCFit.isPC);

            fp.numPCParams = 0;
            
            for(int n=fp.order+1; n<=fp.PCOrder; n++)
                fp.numPCParams += SeriesFitter.getNumParamsForOrder(numPCs,n);
            
            numParams += fp.numPCParams;

            numSamples = sampPerParam*numParams;
        }


        
        //PARAMETER CAP: save time
        if(numParams>2000){
            System.out.println("ABORTING EPICFITTER.DOFIT BECAUSE THERE ARE TOO MANY PARAMETERS: "+numParams);
            return null;
        }

        

        ContSCObjFunction of = (ContSCObjFunction)ccdMin.objFcn;
        DoubleMatrix1D[] sampRel = new DoubleMatrix1D[numSamples];
        DoubleMatrix1D[] sampAbs = new DoubleMatrix1D[numSamples];
        double trueVal[] = new double[numSamples];

        generateSamples(numSamples,sampRel,sampAbs,trueVal);

        //no point doing high1S if all samples are below bCutoff
        boolean allBelowCutoff = true;
        for(int s=0; s<numSamples; s++){
            if(trueVal[s]>=es.EPICThresh1){
                allBelowCutoff = false;
                break;
            }
        }
        
        SparseVDWEnergy sve = null;
        double bCutoffs[] = new double[numSamples];
        double bCutoffs2[] = new double[numSamples];
        double baseShift = 0;//SVE at center
        
        if(fp.explicitVDWCutoff==0){
            Arrays.fill(bCutoffs, es.EPICThresh1);
            Arrays.fill(bCutoffs2, es.EPICThresh2);
        }
        //if using double bCutoff for provability also set bCutoffs array
        //according to that here!
        else{//>0            
            //SparseVDWEnergy sve = new SparseVDWEnergy(of.m,of.efunc,explicitVDWCutoff,of.curDOFVals);
            sve = new SparseVDWEnergy(of,fp.explicitVDWCutoff,sampAbs);
            //((MultiTermEnergyFunction)of.efunc).addTermWithCoeff( sve, -1);
            
            of.setDOFs(center);
            baseShift = sve.getEnergy();//VDW contribution at startVec
            
            
            for(int s=0; s<numSamples; s++){
                of.setDOFs(sampAbs[s]);
                double shift = sve.getEnergy();
                trueVal[s] -= shift - baseShift;
                bCutoffs[s] = es.EPICThresh1 - shift + baseShift;//this keeps bCutoff at the same place
                //accounting for the sparse VDW energy
                
                bCutoffs2[s] = es.EPICThresh2 - shift + baseShift;
            }
        }
        

        double weights[] = null;
            //trying out weighted least squares
        weights = new double[numSamples];
        for(int s=0; s<numSamples; s++){
            if(trueVal[s]>1)
                weights[s] = 1/trueVal[s];
            else
                weights[s] = 1;
        }


        double lambda = 0;

        double[] seriesCoeffs = null;//used for order>2
        EPoly ans = null;

        if(fp.PCOrder>fp.order){//use principal components
            //fit in PC basis
            DoubleMatrix1D ySamp[] = new DoubleMatrix1D[numSamples];
            for(int s=0; s<numSamples; s++)
                ySamp[s] = PCFit.toPCBasis(sampRel[s]);

            if(allBelowCutoff){
                System.out.println("Analytical:");
                PCFit.coeffs = SeriesFitter.fitSeries(ySamp,trueVal,weights,lambda,
                    false,PCFit.fullOrder,PCFit.PCOrder,PCFit.isPC,false,null,null);
            }
            else {
                //ITERATIVE
                PCFit.coeffs = SeriesFitter.fitSeriesIterative(ySamp,trueVal,weights,lambda,
                    false,PCFit.fullOrder,bCutoffs,bCutoffs2,PCFit.PCOrder,PCFit.isPC);
            }
            
            ans = PCFit;
        }
        else{
            
            if(allBelowCutoff){
                System.out.println("Analytical:");
                seriesCoeffs = SeriesFitter.fitSeries(sampRel,trueVal,weights,lambda,
                        false,fp.order);
            }
            else {
                seriesCoeffs = SeriesFitter.fitSeriesIterative(sampRel,trueVal,weights,
                        lambda,false,fp.order,bCutoffs,bCutoffs2,fp.order,null);
            }
            
            ans = new EPoly(DOFNums,numDOFs,DOFmax,DOFmin,center,minE,seriesCoeffs,fp.order);
        }
        
        
        //compress SVE
        if(fp.explicitVDWCutoff>0){
            //store of and m in their current state!
            try{ans.sveOF = (ContSCObjFunction)KSParser.deepCopy(of);}
            catch(Exception e){
                System.err.println("woah sveOF deepCopy failed");
                System.err.println(e.getMessage());
            }
            ans.sveOF.efunc = null;
            
            //won't need these big objects that aren't needed to setDOFs...
            ans.sveOF.de = null;
            ans.sveOF.strandRot = null;//won't need this
            
            sve.m = ans.sveOF.m;
            
            //now clear out the atoms we won't need for more compact storage...
            //we need atoms in residues that are either (1) flexible, or (2) 
            //involved in perturbations that affect flexible residues
            HashSet<Residue> resNeeded = new HashSet<>();
            for(Residue res : sve.m.residue){
                if(res.flexible){
                    resNeeded.add(res);
                    for(int p : res.perts){
                        for(int molResNum : sve.m.perts[p].resDirectlyAffected){
                            //the resDirectlyAffected need to be complete
                            //because they are used to apply the perturbation
                            resNeeded.add(sve.m.residue[molResNum]);
                        }
                    }
                }
            }
            
            for(Residue res : sve.m.residue){
                
                if(!resNeeded.contains(res)){
                    
                    for(int atNum=0; atNum<res.numberOfAtoms; atNum++){
                        int molAtNum = res.atom[atNum].moleculeAtomNumber;
                        res.atom[atNum] = null;
                        sve.m.atom[molAtNum] = null;
                    }
                }
            }
            
            sve.DOFVals = ans.sveOF.curDOFVals;
            ans.sve = sve;

            ans.baseSVE = baseShift;
        }
        
        
        
        //DEBUG!!!!!
            //checking training set mean residual!
        //similar to crossValidateSeries
        /*if(fp.PCOrder<=fp.order){
            double meanResidual = 0;
            double weightSum = 0;


            for(int s=0; s<numSamples; s++){

                DoubleMatrix1D x = sampAbs[s];

                double realVal = ccdMin.objFcn.getValue(x) - ans.minE;
                //actual energy, relative to voxel minimum

                double sampBCutoff = es.EPICThresh1;
                double sampBCutoff2 = es.EPICThresh2;
                double serVal = ans.evaluate(x,false);

                if(ans.sve!=null){
                    ans.sveOF.setDOFs(x);
                    double shift = ans.sve.getEnergy();
                    realVal -= shift - baseShift;
                    sampBCutoff -= shift - baseShift;
                    sampBCutoff2 -= shift - baseShift;
                    serVal -= shift - baseShift;
                    //we subtract off the SVE contribution (shift-baseShift) from everything
                    //(-baseShift because we set the SVE zero point at center)
                }

                double weight = 1;
                if(realVal>1)
                    weight = 1/realVal;//1/Math.min(realVal,sampBCutoff);
                weightSum += weight;


                if(realVal>=sampBCutoff){
                    if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, false)){
                        meanResidual += weight * (serVal-sampBCutoff)*(serVal-sampBCutoff);
                    }
                    if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, true)){
                        meanResidual += weight * (realVal-serVal)*(realVal-serVal);
                    }
                }
                else
                    meanResidual += weight*(realVal-serVal)*(realVal-serVal);
            }

            meanResidual /= weightSum;//numSamples;
            System.out.println("CHECK TRAINING SET MEAN RESIDUAL:"+meanResidual);
        }
        */    //DEBUG!!!!
        
        
        return ans;
    }




    double crossValidateSeries(EPoly fit, FitParams fp){
        //(note this is much like generateSamples!!)
        //return mean residual
        //used in addFitSeries and augmentSeries
        
        if(fit==null)//parameter cap will do this
            return Double.POSITIVE_INFINITY;//obviously an absent fit is no good
        
        double meanResidual = 0;
        double weightSum = 0;

        double baseShift = 0;
        if(fit.sve!=null){
            fit.sveOF.setDOFs(center);
            baseShift = fit.sve.getEnergy();
        }
        
        
        /*double relMax[] = new double[numDOFs];//maximum shifts of degrees of freedom relative to minimum point (startVec)
        double relMin[] = new double[numDOFs];
        for(int dof=0; dof<numDOFs; dof++){
            relMax[dof] = DOFmax.get(dof) - center.get(dof);
            relMin[dof] = DOFmin.get(dof) - center.get(dof);
        }*/
        
        int numSamples = sampPerParam * fp.numParams();
        
        double trueVal[] = new double[numSamples];
        DoubleMatrix1D sampRel[] = new DoubleMatrix1D[numSamples];
        DoubleMatrix1D sampAbs[] = new DoubleMatrix1D[numSamples];
        generateSamples(numSamples,sampRel,sampAbs,trueVal);
        
        for(int s=0; s<numSamples; s++){

            /*DoubleMatrix1D dx = DoubleFactory1D.dense.make(numDOFs);
            DoubleMatrix1D x = DoubleFactory1D.dense.make(numDOFs);
            for(int dof=0; dof<numDOFs; dof++){
                double top = relMax[dof];
                double bottom = relMin[dof];
                dx.set(dof, bottom + Math.random()*(top-bottom));
                x.set(dof, center.get(dof)+dx.get(dof));
            }

            double realVal = ccdMin.objFcn.getValue(x) - fit.minE;
            //actual energy, relative to voxel minimum
            */
            DoubleMatrix1D x = sampAbs[s];
            double realVal = trueVal[s];

            double sampBCutoff = es.EPICThresh1;
            double sampBCutoff2 = es.EPICThresh2;
            double serVal = fit.evaluate(x,false);
            
            if(fit.sve!=null){
                fit.sveOF.setDOFs(x);
                double shift = fit.sve.getEnergy();
                realVal -= shift - baseShift;
                sampBCutoff -= shift - baseShift;
                sampBCutoff2 -= shift - baseShift;
                serVal -= shift - baseShift;
                //we subtract off the SVE contribution (shift-baseShift) from everything
                //(-baseShift because we set the SVE zero point at center)
            }
            
            double weight = 1;
            if(realVal>1)
                weight = 1/realVal;//1/Math.min(realVal,sampBCutoff);
            weightSum += weight;

            
            /*if(serVal<sampBCutoff||realVal<sampBCutoff)
                //meanResidual += (realVal-bv)*(realVal-bv);
                meanResidual += weight * (realVal-serVal)*(realVal-serVal);
            else if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, true))//bCutoff2 upper restraint
                meanResidual += weight*(realVal-serVal)*(realVal-serVal);*/
            if(realVal>=sampBCutoff){
                if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, false)){
                    meanResidual += weight * (serVal-sampBCutoff)*(serVal-sampBCutoff);
                }
                if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, true)){
                    meanResidual += weight * (realVal-serVal)*(realVal-serVal);
                }
            }
            else
                meanResidual += weight*(realVal-serVal)*(realVal-serVal);
        }

        meanResidual /= weightSum;//numSamples;
        System.out.println("CV MEAN RESIDUAL:"+meanResidual);
        
        //Let's return the mean residual
        return meanResidual;
        //END CV
    }



    static void analyzeLSBRecord(ArrayList<double[]> LSBRecord){
        //Given the triples (old lower bound, LSB lower bound, minimized energy)
        //for a bunch of conformations
        //provide statistics on them


        double binMaxs[] = {0.1,0.5,0.75,0.9,0.99,1,1.01,1.1,1.5,Double.POSITIVE_INFINITY};
        int numBins = binMaxs.length;
        double binMins[] = new double[numBins];
        binMins[0] = Double.NEGATIVE_INFINITY;
        System.arraycopy(binMaxs,0,binMins,1,numBins-1);


        double avgSR = 0;//average of slack recovery fractions
        int binCounts[] = new int[numBins];//how many conformations are in each bin
        int numOverHundredth = 0;//Number of LSB bounds > 0.01 kcal/mol over true energy
        int numOverTenth = 0;//>0.1 kcal/mol over
        int numOverHalf = 0;//>0.5 kcal/mol over
        int numEnum=0;//how many conformations we'd need to enumerate
        //using the LSB lower bound



        
        double minMinE = Double.POSITIVE_INFINITY;//Actual GMEC energy
        for(double[] rec: LSBRecord)
            minMinE = Math.min(minMinE,rec[2]);



        for(double[] rec : LSBRecord){

            double slackRecovered = (rec[1]-rec[0])/(rec[2]-rec[0]);
            //fraction of slack recovered by LSB

            avgSR += slackRecovered;

            if(rec[1]<=minMinE)
                numEnum++;
            if(rec[1]>rec[2]+0.01)
                numOverHundredth++;
            if(rec[1]>rec[2]+0.1)
                numOverTenth++;
            if(rec[1]>rec[2]+0.5)
                numOverHalf++;

            for(int bin=0; bin<numBins; bin++){
                if(slackRecovered<=binMaxs[bin]&&slackRecovered>binMins[bin]){
                    binCounts[bin]++;
                    break;
                }
            }
            
            
        }

        avgSR /= LSBRecord.size();

        System.out.println("ANALYSIS OF LSB:");
        System.out.println("Total conformation count: "+LSBRecord.size());
        System.out.println("Average slack recovery fraction: "+avgSR);
        System.out.println(numOverHundredth+" LSBs > 0.01 over true E; "+numOverTenth+" >0.1 over, "
                +numOverHalf+" >0.5 over");
        System.out.println(numEnum+" need to be enumerated based on LSBs");
        System.out.println("Bin_max Bin_count");
        for(int bin=0; bin<binMaxs.length; bin++){
            System.out.println(binMaxs[bin]+" "+binCounts[bin]);
        }
    }



    void generateSamples(int numSamples,
            DoubleMatrix1D[] sampRel, DoubleMatrix1D[] sampAbs, double[] trueVal){

        //Generate samples relative to startVec (sampRel) and absolute (sampAbs)
        //and give their energies (trueVal), relative to baseE
        //we'll try to get half the samples above es.EPICThresh1
        //to prevent overfitting

        //by default, uniform voxel sampling is used
        //shift to Metropolis if we end up drawing numSamples/2 samples over bCutoff
        //and by that time less than maxUniFrac*numSamples/2 have been drawn under bCutoff
        //double maxUniFrac = 0.05;
        
        //Let's instead say if uniform draw fails enough times (maxFailCount) in a row we bail out and do Metropolis
        int maxFailCount = 10000;//big because Metropolis is expensive
        
        //DEBUG!!  To test Metropolis
        //generateSamplesMetropolis(numSamples,sampRel,sampAbs,trueVal,0);
        //return;
        //DEBUG!!!
        
        
        
        int countOverCutoff = 0;//how many of our samples are over the cutoff
        
        ContSCObjFunction of = (ContSCObjFunction)ccdMin.objFcn;
                
        double relMax[] = new double[numDOFs];//maximum shifts of degrees of freedom relative to minimum point (startVec)
        double relMin[] = new double[numDOFs];
        for(int dof=0; dof<numDOFs; dof++){
            relMax[dof] = ccdMin.DOFmax.get(dof) - center.get(dof);
            relMin[dof] = ccdMin.DOFmin.get(dof) - center.get(dof);
        }


        
        for(int s=0; s<numSamples; s++){

            if(countOverCutoff<numSamples/2){//normal draw
                uniformVoxelSample(s,sampRel,sampAbs,trueVal,of,relMin,relMax);
                if(trueVal[s]>es.EPICThresh1)
                    countOverCutoff++;
            }
            else {//force sub-threshold draw
                /*if(s-countOverCutoff<maxUniFrac*countOverCutoff){//not getting enough sub-threshold samples
                    //switch to Metropolis
                    generateSamplesMetropolis(numSamples,sampRel,sampAbs,trueVal,s);
                    break;
                }
                else {//getting enough some-threshold samples to keep up the uniform method
                    */
                    int failCount = 0;//how many uniform samples ended up above the cutoff in a rows
                    
                    do {
                        uniformVoxelSample(s,sampRel,sampAbs,trueVal,of,relMin,relMax);
                        
                        if(failCount>maxFailCount){//bail and do Metropolis
                            generateSamplesMetropolis(numSamples,sampRel,sampAbs,trueVal,s);
                            System.out.println("Drew "+numSamples+" samples of which "+countOverCutoff+" are over bCutoff");
                            return;
                        }
                        failCount++;
                            
                    } while(trueVal[s]>es.EPICThresh1);
                //}
            }
        }
        
        System.out.println("Drew "+numSamples+" samples of which "+countOverCutoff+" are over bCutoff");
    }
            

    void generateSamplesMetropolis(int numSamples,
            DoubleMatrix1D[] sampRel, DoubleMatrix1D[] sampAbs, double[] trueVal, int startSampNum){
        //Use Metropolis (in SubThreshSampler) to get samples
        //way less efficient than uniform sampling if much of the voxel is below EPICThresh1
        //but can hopefully deal with very small hypervolumes being below EPICThresh1
        //this function can start drawing at 
        //so we can have uniformly drawn samples before it (e.g. the above-the-cutoff ones)
        
        SubThreshSampler sts = new SubThreshSampler(es.EPICThresh1+minE,ccdMin.objFcn,DOFmin,DOFmax);
        
        sts.burnIn(center);//start at center because it's known to be below the threshold
        
        for(int s=startSampNum; s<numSamples; s++){
            sampAbs[s] = sts.nextSample();
            trueVal[s] = ccdMin.objFcn.getValue(sampAbs[s]) - minE;
            sampRel[s] = sampAbs[s].copy();
            sampRel[s].assign(center,Functions.minus);
        }
    }
    
    void uniformVoxelSample(int s, DoubleMatrix1D[] sampRel, DoubleMatrix1D[] sampAbs, double[] trueVal,
            ObjectiveFunction of, double[] relMin, double[] relMax){
        //Draw sample # s for generateSamples, filling in trueVal, sampRel, and sampAbs
        //Generate vector relative to minimum
        DoubleMatrix1D dx = DoubleFactory1D.dense.make(numDOFs);
        //and absolute
        DoubleMatrix1D x = DoubleFactory1D.dense.make(numDOFs);

        for(int dof=0; dof<numDOFs; dof++){
            double top = relMax[dof];
            double bottom = relMin[dof];

            dx.set(dof, bottom + Math.random()*(top-bottom));
            x.set(dof, center.get(dof)+dx.get(dof));
        }

        trueVal[s] = of.getValue(x) - minE;

        sampRel[s] = dx;
        sampAbs[s] = x;
    }
    
    
    
    
    //Default sample generation is not limited to some scale
    /*void generateSamples(int numSamples,
            DoubleMatrix1D[] sampRel, DoubleMatrix1D[] sampAbs, double[] trueVal){

        generateSamples(numSamples,sampRel,sampAbs,trueVal,0,null,0);
    }


    void generateSamples( int numSamples, 
            DoubleMatrix1D[] sampRel, DoubleMatrix1D[] sampAbs, double[] trueVal, 
            double scale, ContETerm baseFit, double errorThresh){
        //Generate samples relative to startVec (sampRel) and absolute (sampAbs)
        //and give their energies (trueVal), relative to baseE
        //if scale!=0 we stay within distance scale of the center in each dimension
        //if baseFit==null, then we'll try to get half the samples above bCutoff
        //if baseFit!=null, then we'll try to get half the samples above bCutoff with
        //baseBound being off by at least errorThresh
        //(these are to prevent overfitting by avoiding having too few samples below bCutoff
        // or, for fit augmentation, too few samples really near the base fit)


        boolean bThresh = true;//make sure half the samples are below bCutoff
        int maxTriesForThresh = 10000;//how many times we're willing to redraw each sample to meet thresholds

        int countOverCutoff = 0;//how many of our samples are over the cutoff
        int countBelowErrorThresh = 0;//how many are below errorThresh
        int countBeyondThresh = 0;//how many are either above cutoff or below error thresh
        //(we'll aim for countBeyondThresh being half the samples or less)

        ContSCObjFunction of = (ContSCObjFunction)ccdMin.objFcn;
                
        double relMax[] = new double[numDOFs];//maximum shifts of degrees of freedom relative to minimum point (startVec)
        double relMin[] = new double[numDOFs];
        for(int dof=0; dof<numDOFs; dof++){
            relMax[dof] = ccdMin.DOFmax.get(dof) - center.get(dof);
            relMin[dof] = ccdMin.DOFmin.get(dof) - center.get(dof);
        }


        for(int s=0; s<numSamples; s++){

            //COULD TAKE THINGS OUT OF LOOP...

            int triesForThresh = 0;//how many times we try to get under bCutoff and under errorThresh
            int triesAboveCutoff = 0;//how many of the tries were above bCutoff (used to decide if we should reduce scale)

            boolean needRedraw = false;//sample needs to be redrawn
            //so we can get enough samples

            boolean belowErrThresh = false;//sample is below error threshold (if present)

            do{
                //Generate vector relative to minimum
                DoubleMatrix1D dx = DoubleFactory1D.dense.make(numDOFs);
                //and absolute
                DoubleMatrix1D x = DoubleFactory1D.dense.make(numDOFs);

                for(int dof=0; dof<numDOFs; dof++){
                    double top = relMax[dof];
                    double bottom = relMin[dof];

                    if(scale!=0){
                        top = Math.min(top, scale);
                        bottom = Math.max(bottom, -scale);
                    }

                    dx.set(dof, bottom + Math.random()*(top-bottom));
                    x.set(dof, center.get(dof)+dx.get(dof));
                }

                trueVal[s] = of.getValue(x) - minE;

                sampRel[s] = dx;
                sampAbs[s] = x;

                needRedraw = false;
                belowErrThresh = false;
                    
                if(trueVal[s]>=es.EPICThresh1){
                    if(bThresh){
                        needRedraw = true;
                        triesAboveCutoff++;
                    }
                }
                else if(baseFit != null){
                    //error thresh applies only if trueVal[s] is below bCutoff
                    double baseFitVal = baseFit.evaluate(x,false);
                    if( (trueVal[s]-baseFitVal)*(trueVal[s]-baseFitVal) < errorThresh ){
                        needRedraw = true;
                        belowErrThresh = true;
                    }
                }


                triesForThresh++;
                
            } while( needRedraw && (triesForThresh<maxTriesForThresh)
                    && (countBeyondThresh>=numSamples/2) );


            if(trueVal[s] >= es.EPICThresh1)
                countOverCutoff++;
            if(belowErrThresh)
                countBelowErrorThresh++;
            if(trueVal[s]>=es.EPICThresh1 || belowErrThresh)
                countBeyondThresh++;
        }


        System.out.println("Drew "+numSamples+" samples of which "+countOverCutoff+" are over bCutoff, "
                +countBelowErrorThresh+" under errorThresh (if applicable), "+countBeyondThresh+" one or the other");
        
        if(countBeyondThresh>numSamples/2)
            System.out.println("Warning: Maxed out trying to get enough sub-threshold samples!");
    }*/

    
    boolean useSVE = true;
    boolean usePC = true;

    
    FitParams raiseFitOrder(FitParams fp){
        //given a FitParams object, return a higher-order one
        //to follow the standard EPIC progression of fit orders
        //(modified by useSVE, usePC as needed)
        //all these fits are meant to be performed iteratively (with thresholds) though
        //analytical evaluation will kick in if no samples exceed the first threshold
        //this function is intended to be used repetitively until a good enough fit is found, starting
        //with a FitParams.quadratic()
        
        //since these are standard EPIC fits we don't includeConst
        
        
        if(useSVE){
            //we try SVE fits after non-PC, polynomial-only fits of the same order
            if(fp.explicitVDWCutoff==0 && fp.order>=fp.PCOrder){
                double cutoff = 3;
                if(fp.order>2)//4 or 6
                    cutoff = 4;
                return new FitParams(numDOFs,fp.order,0,fp.order,false,cutoff);
            }
        }
        
        
        if(usePC && fp.order<6){
            //after quadratic or quartic fits without principal components (SVE or not)
            //we try first PCFac=0.1 and then 0.01
            if(fp.PCOrder==fp.order)//start PC with 0.1
                return new FitParams(numDOFs,fp.order,0.1,fp.order+2,false,0);
            else if(fp.PCFac==0.1)//move on to 0.01
                return new FitParams(numDOFs,fp.order,0.01,fp.order+2,false,0);
        }
        
        
        //if we get here we need to raise the main degree of the polynomial
        //we go 2, 4, 6, and then give up
        
        if(fp.order>=6)//give up
            return null;
        else//raise order by 2
            return new FitParams(numDOFs,fp.order+2,0,fp.order+2,false,0);
    }
    
    
    
    EPoly blank(){
        //return a EPoly on no continuous degrees of freedom
        return new EPoly(DOFNums,numDOFs,DOFmax,DOFmin,center,minE,null,2);
        //arbitrarily calling it quadratic (doesn't matter since no variables in polynomial)
    }
    

}
