/*
 * JMCSHandler.java
 *
 * Created on January 28, 2007, 11:37 AM
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package chemlib.ebi.molecule.mcs.helper;

//~--- JDK imports ------------------------------------------------------------
import chemlib.ebi.core.tools.EBIAtomContainerManipulator;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;
import java.util.Vector;
//~--- non-JDK imports --------------------------------------------------------

import chemlib.ebi.core.tools.descriptors.MoleculeSanityCheck;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.BondTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.MDLReader;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK  * @contact asad@ebi.ac.uk
 */
public class MolHandler {

    private IAtomContainer Mol = null;
    private IAtomContainerSet FragmentMolSet = null;
    private int BondNumber = 0;
    private int AtomNumber = 0;
    private int startHatom_num = 0;
    private Vector<String> AtomString = null;
    private Vector<String> char_table = null;
    //With H
    private Vector<Integer> int_table = null;
    //with H
    private Vector<IAtom> AtomsWithH = null;
    //with H
    private Vector<String> char_tab = null;
    //without H
    private Vector<Integer> int_tab = null;
    //without H
    private Vector<IAtom> AtomsWithoutH = null;
    //without H
    private boolean FragmentFlag;
    private boolean removeHydrogen = true;

    /** Creates a new instance of JMCSHandler
     * @param MolFile Mol file name
     * @param cleanMolecule
     * @param removeHydrogen 
     *  
     */
    public MolHandler(String MolFile, boolean cleanMolecule, boolean removeHydrogen) {

        MDLReader MolRead;

        //Mol2Reader Mol2Read;
        this.removeHydrogen = removeHydrogen;
        char_tab = new Vector<String>();
        int_tab = new Vector<Integer>();

        //if (checkForH(Mol) > 0) {



        try {


            FileInputStream ReadMolecule;

            ReadMolecule = new FileInputStream(MolFile);
            MolRead = new MDLReader(new InputStreamReader(ReadMolecule));
            Mol = (IMolecule) MolRead.read(new Molecule());


            if (cleanMolecule) {
                MoleculeSanityCheck.fixAromaticity((IMolecule) Mol);
            }


            BondTools.makeUpDownBonds(Mol);


//            /*Remove Hydrogen by Asad*/
            if (removeHydrogen) {
//                removeHydrogenAtoms(Mol);
                Mol = EBIAtomContainerManipulator.removeHydrogens(Mol);

            }

            if (Mol.getAtomCount() > 0) {
                FragmentFlag = ConnectivityChecker.isConnected(Mol);
            }
            FragmentMolSet = DefaultChemObjectBuilder.getInstance().newMoleculeSet();


            //System.out.println("isConnected : " + FragmentFlag);

            if (!FragmentFlag) {

                /*System.err.println("The molecule is not connected, " +
                "Fragment Matcher Will Handle this _molecule");*/

                FragmentMolSet.add(ConnectivityChecker.partitionIntoMolecules(Mol));
                FragmentMolSet.setID(Mol.getID());

            } else {

                FragmentMolSet.addAtomContainer(Mol);
                FragmentMolSet.setID(Mol.getID());
            }

            boolean hydrogenFlag = true;



            BondNumber = Mol.getBondCount();
            AtomNumber = Mol.getAtomCount();

            setAtomString();
            setBondConnectionTable();
            setSymbolConnectionTable();


            enumerate_startHatom_num();




            if ((AtomNumber == 2) && (BondNumber == 1)) {


                if ((char_table.get(0).compareToIgnoreCase("H") == 0) && (char_table.get(1).compareToIgnoreCase("H") == 0)) {
                    {
                        System.out.println("H2 Molekuel");
                        hydrogenFlag = false;
                    }
                }
            }

            if (removeHydrogen) {
                discard_H_bonds(hydrogenFlag);


                //if (checkForH(Mol) > 0) {
            } else {
                use_H_bonds();

            }


        } catch (IOException ex) {
            Logger.getLogger(MolHandler.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException e) {
            System.err.println(e);
        }
    }

    /**
     * 
     * @param _molecule Molecule AtomContainer
     * @param cleanMolecule
     * @param removeHydrogen
     */
    public MolHandler(IAtomContainer _molecule, boolean cleanMolecule, boolean removeHydrogen) {

        String ID = _molecule.getID();
        this.removeHydrogen = removeHydrogen;

        char_tab = new Vector<String>();
        int_tab = new Vector<Integer>();


        this.Mol = _molecule;
        boolean no_H2 = true;


        if (cleanMolecule) {
            MoleculeSanityCheck.fixAromaticity((IMolecule) Mol);
        }


        //
//         /*Hydrogen are always removed for this molecule before mapping*/

        if (removeHydrogen) {
            try {
                Mol = (IMolecule) EBIAtomContainerManipulator.removeHydrogensAndPreserveAtomID(Mol);
                Mol.setID(ID);
//            removeHydrogenAtoms(Mol);
//            System.out.println("Bonds: " + Mol.getBondCount());
            } catch (Exception ex) {
                Logger.getLogger(MolHandler.class.getName()).log(Level.SEVERE, null, ex);
            }
//            removeHydrogenAtoms(Mol);
//            System.out.println("Bonds: " + Mol.getBondCount());

        }

//        System.out.println("AtomCount : " + Mol.getAtomCount());
        if (Mol.getAtomCount() > 0) {
            FragmentFlag = ConnectivityChecker.isConnected(Mol);
        }
        FragmentMolSet = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
//        System.out.println("isConnected : " + FragmentFlag);
        if (!FragmentFlag) {

            /*System.err.println("The molecule is not connected, " +
            "Fragment Matcher Will Handle this _molecule");*/


            FragmentMolSet.add(ConnectivityChecker.partitionIntoMolecules(Mol));
            FragmentMolSet.setID(Mol.getID());
        } else {

            FragmentMolSet.addAtomContainer(Mol);
            FragmentMolSet.setID(Mol.getID());
        }

        //System.out.println(FragmentMolSet.getMoleculeCount());

        BondNumber = Mol.getBondCount();
        AtomNumber = Mol.getAtomCount();

        setAtomString();
        setBondConnectionTable();
        setSymbolConnectionTable();
        enumerate_startHatom_num();


        if ((AtomNumber == 2) && (BondNumber == 1)) {


            if ((char_table.get(0).compareToIgnoreCase("H") == 0) && (char_table.get(1).compareToIgnoreCase("H") == 0)) {
                {
                    System.out.println("H2 Molekuel");
                    no_H2 = false;
                }
            }
        }

        if (removeHydrogen) {
            discard_H_bonds(no_H2);
        } else {
            use_H_bonds();
        }


    }

    /**
     * 
     * @return Bond count
     */
    public int getBondCount() {
        return BondNumber;
    }

    /**
     * 
     * @return Atom count
     */
    public int getAtomCount() {
        return AtomNumber;
    }

    /**
     * 
     * @return atoms as List
     */
    public List<String> getAtomString() {
        //System.out.println("Atom String returned from _molecule: " + char_tab  + " Size:" + char_tab.size());
        return AtomString;
    }

    /**
     * 
     * @return Hydrogen count
     */
    public int getHydrogenNumberStart() {
        return startHatom_num;
    }

    /**
     * 
     * @return bond connection as vector
     */
    public Vector<Integer> getBondConnectionTable() {
        return int_table;
    }

    /**
     * 
     * @return atoms connection table
     */
    public Vector<String> getSymbolConnectionTable() {
        return char_table;
    }

    /**
     * 
     * @return Atoms including Hydrogens
     */
    public Vector<IAtom> getAtomStereoParityConnectionTable() {
        return AtomsWithH;
    }

    /**
     * 
     * @return Bonds between atoms including Hydrogens
     */
    public Vector<Integer> getBondConnectionTableWithoutHydrogen() {
        return int_tab;
    }

    /**
     * 
     * @return Atoms excluding Hydrogens
     */
    public Vector<String> getSymbolConnectionTableWithoutHydrogen() {

        return char_tab;
    }

    /**
     * 
     * @return Bonds between atoms excluding Hydrogens
     */
    public Vector<IAtom> getAtomStereoParityConnectionTableWithoutHydrogen() {
        return AtomsWithoutH;
    }

    private void enumerate_startHatom_num() {
        /*startHatom_num = AtomNumber;
        //System.out.println("Atom Number " + AtomNumber);
        for (int atom = 0; atom < AtomNumber; atom++) {
        if ((Mol.getAtom(atom).getSymbol()).equals("H")) {
        startHatom_num = atom;
        }
        }*/

        startHatom_num = AtomNumber;
        int a = AtomNumber - 1;

        while ((a >= 0) && (AtomString.get(a).equals("H"))) {
            startHatom_num--;
            a--;
        }
    }

    private void setAtomString() {

        AtomString = new Vector<String>();

        for (int atom = 0; atom < AtomNumber; atom++) {
            String symbol = Mol.getAtom(atom).getSymbol();

            // System.out.println(symbol);

            AtomString.add(symbol);
        }

        //System.err.println("In Mol: getString " +char_tab.size()+ " "+ char_tab);
        //this.char_tab = getString(temp);
        //char_tab= temp;//We have changed temp as Array list instead of String ASAD

    }

    private void setBondConnectionTable() {

        int_table = new Vector<Integer>();

        IAtomContainer ac = Mol;

        //System.out.println();

        for (int i = 0; i < ac.getBondCount(); i++) {
            IBond bond = ac.getBond(i);

            /* This will fetch the connected ATOM as integer and its Bond order ex: 2 as double, 1 as single */

            // System.out.println(ac.getAtomCount(bond.getAtom(0))+" "+ac.getAtomCount(bond.getAtom(1))+" "+(int)bond.getOrder());
            int_table.add(ac.getAtomNumber(bond.getAtom(0)));    // Plus one because Java Indexing is one less

            int_table.add(ac.getAtomNumber(bond.getAtom(1)));    // Plus one because Java indexing is one less
            //Uncomment for new version*/

            Order BondOrder = bond.getOrder();
            int_table.add(BondOrder.ordinal() + 1);
            //System.out.println(BondOrder + " = " + BondOrder.ordinal());

            //Uncomment for old cdk 1.01 version
            //Double d = bond.getOrder();
            //int_table.add(d.intValue());

        }




        /* This will fetch the Connected ATOM Symbol */

        // System.out.println(bond.getAtom(0).getSymbol()+" "+bond.getAtom(1).getSymbol());

    }

    private int checkForH(IAtomContainer mol) {
        int hCount = 0;


        for (int i = 0; i <
                mol.getAtomCount(); i++) {


            if (mol.getAtom(i).getSymbol().equals("H")) {
                hCount++;
            }

        }


        return hCount;
    }

    private void setSymbolConnectionTable() {
        // char_table = new Vector<Character>();
        char_table = new Vector<String>();
        AtomsWithH =
                new Vector<IAtom>();
        AtomsWithoutH =
                new Vector<IAtom>();

        //System.out.println();

        IAtomContainer ac = Mol;

        for (int i = 0; i <
                ac.getBondCount(); i++) {
            IBond bond = ac.getBond(i);

            /* This will fetch the Connected ATOM Symbol */
            IAtom Atom1 = bond.getAtom(0);
            IAtom Atom2 = bond.getAtom(1);


            char_table.add(Atom1.getSymbol());

            char_table.add(Atom2.getSymbol());




            // System.out.println(bond.getAtom(0).getSymbol()+" "+bond.getAtom(1).getSymbol());
            // throw new AssertionError("Error");
        }

        for (int i = 0; i <
                ac.getAtomCount(); i++) {

            //System.out.println(" Parity" + ac.getAtomParity(ac.getAtom(i)));

            //ac.getAtom(i).setStereoParity();
            AtomsWithH.add(ac.getAtom(i));


            if (!ac.getAtom(i).getSymbol().equalsIgnoreCase("H")) {
                AtomsWithoutH.add(ac.getAtom(i));
            }

        }



    }

    private int use_H_bonds() {

        int count_bonds = 0;

        for (int x = 0; x <
                BondNumber; x++) {

            char_tab.add(char_table.get(x * 2 + 0));
            char_tab.add(char_table.get(x * 2 + 1));

            //AtomsWithoutH.add(AtomsWithH.get(x * 2 + 0));
            //AtomsWithoutH.add(AtomsWithH.get(x * 2 + 1));

            int_tab.add(int_table.get(x * 3 + 0));
            int_tab.add(int_table.get(x * 3 + 1));
            int_tab.add(int_table.get(x * 3 + 2));


            count_bonds++;
        }


        BondNumber = count_bonds;

        return 0;
    }

    private int discard_H_bonds(boolean is_no_H2) {

        int count_bonds = 0;

        if (is_no_H2) {
            for (int x = 0; x <
                    BondNumber; x++) {
                if ((char_table.elementAt(x * 2 + 0).equalsIgnoreCase("H")) || (char_table.elementAt(x * 2 + 1).equalsIgnoreCase("H"))) {
                    AtomNumber--;
                }

                if ((!char_table.elementAt(x * 2 + 0).equalsIgnoreCase("H")) && (!char_table.elementAt(x * 2 + 1).equalsIgnoreCase("H"))) {
                    char_tab.add(char_table.get(x * 2 + 0));
                    char_tab.add(char_table.get(x * 2 + 1));

                    //AtomsWithoutH.add(AtomsWithH.get(x * 2 + 0));
                    //AtomsWithoutH.add(AtomsWithH.get(x * 2 + 1));

                    int_tab.add(int_table.get(x * 3 + 0));
                    int_tab.add(int_table.get(x * 3 + 1));
                    int_tab.add(int_table.get(x * 3 + 2));

                    count_bonds++;

                }


            }
        } else {

            for (int x = 0; x <
                    BondNumber; x++) {

                char_tab.add(char_table.get(x * 2 + 0));
                char_tab.add(char_table.get(x * 2 + 1));

                //AtomsWithoutH.add(AtomsWithH.get(x * 2 + 0));
                //AtomsWithoutH.add(AtomsWithH.get(x * 2 + 1));

                int_tab.add(int_table.get(x * 3 + 0));
                int_tab.add(int_table.get(x * 3 + 1));
                int_tab.add(int_table.get(x * 3 + 2));


                count_bonds++;

            }


        }

        BondNumber = count_bonds;

        return 0;
    }

    public IAtomContainer getMolecule() {

//        System.out.println("MCSMolHandle Mol:" + Mol.getID() + "  " + Mol.getAtomCount());
        return Mol;
    }

    public boolean getFragmentFlag() {

        return FragmentFlag;
    }

    /**
     * 
     * @return
     */
    public IAtomContainerSet getFragmentedMolecule() {

        return FragmentMolSet;
    }
}
//~ Formatted by Jindent --- http://www.jindent.com

