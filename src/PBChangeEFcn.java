
import java.io.BufferedReader;
import java.io.InputStreamReader;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author mhall44
 */
public class PBChangeEFcn extends EnergyFunction {
    //This energy function is meant to capture continuous variation in Poisson-Boltzmann energies due to dihedral angle changes
    //(so the molecule involved is not expected to have mutations)
    //we start with a structure (ideal rotamers)
    //then for intra+shell energies, we return Poisson-Boltzmann energy of the molecule w/ changed dihedrals (expected to be at one res) - orig E
    //for pairwise, we do total PB E - (E with only res1 changed) - (E with only res2 changed) + orig E
    //since this is just a first investigation into representing Poisson-Boltzmann energies,
    //it is not particularly speed-optimized.  PB calc will likely be the bottleneck anyway.  
    
    Molecule m;//The molecule whose conformational changes we will be calculating energies for
    Molecule mCopy;//A copy of the molecule used for outputting PDB files, which Delphi will read
    double origE;
    
    boolean pairwise = false;//is this a pairwise energy?
    //if so here are the molecule residue numbers for the pair:
    int res1;
    int res2;
    
    static String delphiFolder = "OSPREY_delphi";
    
    public PBChangeEFcn(Molecule molec, int str1, int strResNum1, int str2, int strResNum2){
        //if pairwise, then we must have valid strand and strand res. num. arguments
        //non-pairwise flagged by str2==-1 (e.g. this is found for shell runs in RotamerSearch)
        
        m = molec;
        
        try{
            mCopy = (Molecule)KSParser.deepCopy(molec);
        }
        catch(Exception e){
            throw new RuntimeException("ERROR deepCopying molec in new PBChangeEFcn()");
        }
        
        mCopy.backupCoordinates = new double[mCopy.actualCoordinates.length];
        fullCopy(mCopy.actualCoordinates,mCopy.backupCoordinates);
        
        origE = mCopyE();
        
        if(str2!=-1){//pairwise
            pairwise = true;
            res1 = m.strand[str1].residue[strResNum1].moleculeResidueNumber;
            res2 = m.strand[str2].residue[strResNum2].moleculeResidueNumber;
        }
    }

    
    @Override
    public double getEnergy() {
        if(pairwise){
            fullCopy(m.actualCoordinates,mCopy.actualCoordinates);
            double ans = mCopyE();
            ans -= singleResE(res1);
            ans -= singleResE(res2);
            ans += origE;
            return ans;
        }
        else {
            //copy over coordinate and evaluate energy
            fullCopy(m.actualCoordinates,mCopy.actualCoordinates);
            return mCopyE()-origE;
        }
    }
    
    double singleResE(int molResNum){
        //create a hybrid molecule that has coordinates of m at molResNum, original (backed-up) coordinates elsewhere
        fullCopy(mCopy.backupCoordinates,mCopy.actualCoordinates);
        Residue res = m.residue[molResNum];
        for(Atom at : res.atom){
            int offset = 3*at.moleculeAtomNumber;
            System.arraycopy(m.actualCoordinates, offset, mCopy.actualCoordinates, offset, 3);
        }
        return mCopyE();   
    }
    
    
    void fullCopy(double[] a, double[] b){
        //full copy of double array
        System.arraycopy(a,0,b,0,a.length);
    }
    
    
    double mCopyE(){
        //get the current PB energy of mCopy
        //mCopy doesn't use Atom.coord except when saving the molecule, so we can just output actualCoordinates using saveMolec
        //without messing anything up
        
        //save the molecule to a PDB file and run Delphi on it
        new KSParser().saveMolecule(mCopy, delphiFolder+"/struct.pdb", 0);
        double E = 0;
        
        try{
            Process p = Runtime.getRuntime().exec(delphiFolder+"/getE");
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            E = Double.valueOf( br.readLine().trim() );
            E *= RotamerSearch.constRT;//convert from thermal energies (Delphi unit) to kcal/mol (OSPREY unit)
            br.close();
        }
        catch(Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        return E;
    }
    
    
    

    @Override
    public void setupPartialComputation(int[][] residues) {
        //not doing partial comp...would require PB perturbation theory probably
    }

    @Override
    public double getEnergy(int a) {
        return getEnergy();
    }
    
    
    
}
