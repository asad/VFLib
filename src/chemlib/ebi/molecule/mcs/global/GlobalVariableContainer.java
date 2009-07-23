/*
 * GlobalVariableContainer.java
 *
 * Created on February 8, 2007, 9:56 AM
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package chemlib.ebi.molecule.mcs.global;

import java.io.ObjectStreamException;
import java.util.List;
import java.util.Vector;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

public class GlobalVariableContainer {

    private int atom_number1;
    private int atom_number2;
    private int atom_num_H_1;
    private int atom_num_H_2;
    private int bond_number1;
    private int bond_number2;
    private static List<String> atomstr1 = null;
    private static List<String> atomstr2 = null;
    private static Vector<Integer> i_tab1 = null;
    private static Vector<Integer> i_tab2 = null;
    private static Vector<String> c_tab1 = null;
    private static Vector<String> c_tab2 = null;
    private static Vector<IAtom> ReactantStereoParityConnectionTable = null;
//    private static Vector<IAtom> ReactantStereoParityConnectionTableWithoutHydrogen = null;
    private static Vector<IAtom> ProductStereoParityConnectionTable = null;
//    private static Vector<IAtom> ProductStereoParityConnectionTableWithoutHydrogen = null;
    private Integer best_MAPPING_size = 0;
    private Integer clique_count = 0;
    private String osName = null;
    private static IAtomContainer ReactantMol = null;
    private static IAtomContainer ProductMol = null;
    /** Creates a new instance of GlobalVariableContainer */
    private static GlobalVariableContainer INSTANCE = null;

    /**
     * 
     * @return
     */
    public static synchronized GlobalVariableContainer getInstance() {
        if (INSTANCE == null) {

            // it's ok, we can call this constructor
            INSTANCE = new GlobalVariableContainer();
        }

        return INSTANCE;
    }

    protected GlobalVariableContainer() {

        i_tab1 = new Vector<Integer>();
        i_tab2 = new Vector<Integer>();
        c_tab1 = new Vector<String>();
        c_tab2 = new Vector<String>();


        ReactantStereoParityConnectionTable = new Vector<IAtom>();
//        ReactantStereoParityConnectionTableWithoutHydrogen = new Vector<IAtom>();
        ProductStereoParityConnectionTable = new Vector<IAtom>();
//        ProductStereoParityConnectionTableWithoutHydrogen = new Vector<IAtom>();

        this.best_MAPPING_size = 0;
        this.clique_count = 0;

        this.atom_number1 = 0;
        this.atom_number2 = 0;
        this.atom_num_H_1 = 0;
        this.atom_num_H_2 = 0;
        this.bond_number1 = 0;
        this.bond_number2 = 0;

        ReactantMol = null;
        ProductMol = null;

        best_MAPPING_size = 0;
        clique_count = 0;


        this.osName = System.getProperty("os.name");
        //System.out.println("OS Name: " + osName);


    }

    public void Clear() {
        i_tab1 = new Vector<Integer>();
        i_tab2 = new Vector<Integer>();
        c_tab1 = new Vector<String>();
        c_tab2 = new Vector<String>();

        ReactantStereoParityConnectionTable = new Vector<IAtom>();
//        ReactantStereoParityConnectionTableWithoutHydrogen = new Vector<IAtom>();
        ProductStereoParityConnectionTable = new Vector<IAtom>();
//        ProductStereoParityConnectionTableWithoutHydrogen = new Vector<IAtom>();

        this.best_MAPPING_size = 0;
        this.clique_count = 0;

        this.atom_number1 = 0;
        this.atom_number2 = 0;
        this.atom_num_H_1 = 0;
        this.atom_num_H_2 = 0;
        this.bond_number1 = 0;
        this.bond_number2 = 0;

        ReactantMol = null;
        ProductMol = null;
        best_MAPPING_size = 0;
        clique_count = 0;

        this.osName = System.getProperty("os.name");

    }

    /**
     * 
     * @param atomStereoParityConnectionTable
     */
    public void setProductAtomStereoParity(Vector<IAtom> atomStereoParityConnectionTable) {
        ReactantStereoParityConnectionTable = atomStereoParityConnectionTable;
    }

//    /**
//     *
//     * @param atomStereoParityConnectionTableWithoutHydrogen
//     */
//    public void setProductAtomStereoParityWithoutHydrogen(Vector<IAtom> atomStereoParityConnectionTableWithoutHydrogen) {
//        ReactantStereoParityConnectionTableWithoutHydrogen = atomStereoParityConnectionTableWithoutHydrogen;
//    }
    public void setProductMolecule(IAtomContainer molecule) {
        ProductMol = molecule;
    }

    /**
     * 
     * @param atomStereoParityConnectionTable
     */
    public void setReactantAtomStereoParity(Vector<IAtom> atomStereoParityConnectionTable) {
        ProductStereoParityConnectionTable = atomStereoParityConnectionTable;
    }

//    /**
//     *
//     * @param atomStereoParityConnectionTableWithoutHydrogen
//     */
//    public void setReactantAtomStereoParityWithoutHydrogen(Vector<IAtom> atomStereoParityConnectionTableWithoutHydrogen) {
//        ProductStereoParityConnectionTableWithoutHydrogen = atomStereoParityConnectionTableWithoutHydrogen;
//    }
    public void setReactantMolecule(IAtomContainer molecule) {
        ReactantMol = molecule;
    }

    private Object readResolve() throws ObjectStreamException {
        return INSTANCE;
    }

    /**
     * 
     * @param atom_number
     */
    public void setReactantAtomSize(int atom_number) {

        this.atom_number1 = atom_number;

    }

    /**
     * 
     * @param atom_number
     */
    public void setProductAtomSize(int atom_number) {

        this.atom_number2 = atom_number;
    }

    /**
     * 
     * @param atom_num_H
     */
    public void setReactantHydrogenCount(int atom_num_H) {

        this.atom_num_H_1 = atom_num_H;
    }

    /**
     * 
     * @param atom_num_H
     */
    public void setProductHydrogenCount(int atom_num_H) {

        this.atom_num_H_2 = atom_num_H;
    }

    /**
     * 
     * @param bond_number
     */
    public void setReactantBondNumber(int bond_number) {
        this.bond_number1 = bond_number;
    }

    /**
     * 
     * @param bond_number
     */
    public void setProductBondNumber(int bond_number) {
        this.bond_number2 = bond_number;
    }

    /**
     * 
     * @param atomstr
     */
    public void setReactantElement(List<String> atomstr) {
        atomstr1 = java.util.Collections.synchronizedList(new Vector<String>(atomstr));
    }

    /**
     * 
     * @param atomstr
     */
    public void setProductElement(List<String> atomstr) {
        atomstr2 = java.util.Collections.synchronizedList(new Vector<String>(atomstr));
    }

    /**
     * 
     * @param i_tab
     */
    public void setReactantBondConnectionTable(Vector<Integer> i_tab) {
        i_tab1 = i_tab;
    }

    /**
     * 
     * @param i_tab
     */
    public void setProductBondConnectionTable(Vector<Integer> i_tab) {
        i_tab2 = i_tab;
    }

    /**
     * 
     * @param c_tab
     */
    public void setReactantAtomConnectionTable(Vector<String> c_tab) {
        //System.out.println("Atom String Educt set in GlobalVarialeContainer: " + c_tab + " Size:" + c_tab.size());
        c_tab1 = c_tab;
    }

    /**
     * 
     * @param c_tab
     */
    public void setProductAtomConnectionTable(Vector<String> c_tab) {
        // System.out.println("Atom String Product set in GlobalVarialeContainer: " + c_tab  + " Size:" + c_tab.size());
        c_tab2 = c_tab;
    }

    /**
     * 
     * @param MAPPING_size 
     */
    public void setCliqueSize(int MAPPING_size) {
        this.best_MAPPING_size = 0;
        this.best_MAPPING_size = MAPPING_size;
    }

    /**
     * 
     * @param clique_size
     */
    public void setCliqueCount(int clique_size) {
        //System.out.println("Best Size was Set to:" + clique_size);
        this.clique_count = 0;
        this.clique_count = clique_size;
    }

    /*******************
     * Get Methods*
     ******************
     * @return 
     */
    public int getMappingSize() {
        return this.best_MAPPING_size;
    }

    /**
     * 
     * @return
     */
    public int getCliqueCount() {
        return this.clique_count;
    }

    /**
     * 
     * @return
     */
    public int getReactantAtomSize() {

        return atom_number1;

    }

    /**
     * 
     * @return
     */
    public int getProductAtomSize() {

        return atom_number2;
    }

    /**
     * 
     * @return
     */
    public int getReactantHydrogenNumber() {

        return atom_num_H_1;
    }

    /**
     * 
     * @return
     */
    public int getProductHydrogenNumber() {

        return atom_num_H_2;
    }

    /**
     * 
     * @return
     */
    public int getReactantBondNumber() {
        return bond_number1;
    }

    /**
     * 
     * @return
     */
    public int getProductBondNumber() {
        return bond_number2;
    }

    public IAtomContainer getReactantMolecule() {
        return ReactantMol;
    }

    public IAtomContainer getProductMolecule() {
        return ProductMol;
    }

    /**
     * 
     * @return
     */
    public List<String> getReactantElement() {


        return atomstr1;
    }

    /**
     * 
     * @return
     */
    public List<String> getProductElement() {

        return atomstr2;
    }

    /**
     * 
     * @return
     */
    public Vector<Integer> getReactantBondConnectionTable() {
        return i_tab1;
    }

    /**
     * 
     * @return
     */
    public Vector<Integer> getProductBondConnectionTable() {
        return i_tab2;
    }

    /**
     * 
     * @return
     */
    public Vector<String> getReactantAtomConnectionTable() {
        return c_tab1;
    }

    /**
     * 
     * @return
     */
    public Vector<String> getProductAtomConnectionTable() {
        return c_tab2;
    }

    /**
     * 
     * @return 
     */
    public Vector<IAtom> getProductAtom() {
        return ReactantStereoParityConnectionTable;
    }

    /**
     * 
     * @return 
     */
    public Vector<IAtom> getReactantAtom() {
        return ProductStereoParityConnectionTable;
    }

    /**
     * 
     * @return
     */
    public String getOSName() {

        return osName;
    }

}
