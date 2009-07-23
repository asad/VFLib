/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.mcs.algorithm.vflib.mcs;

/*
 * McGregorVFLib.java
 *
 * Created on 01 February 2007, 10:51
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
import chemlib.ebi.molecule.mcs.global.BondType;
import chemlib.ebi.molecule.mcs.helper.BinaryTree;
import java.io.IOException;
import java.util.Map;
import java.util.Stack;
import java.util.Vector;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK  * @contact asad@ebi.ac.uk
 */
public class McGregorVFLib {

    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;
    private BinaryTree last;
    private BinaryTree first;
    private Stack<Vector<Integer>> BESTARCS;
    private Vector<Integer> MARCS;
    private Vector<Integer> i_globalA;
    private Vector<Integer> i_globalB;
    private Vector<String> c_globalA;
    private Vector<String> c_globalB;
    Vector<String> c_tab1_copy;
    Vector<String> c_tab2_copy;
    private int nNum_globalA;
    private int nNum_globalB;
    private int bestarcsleft;
    private int globalMCSSize = 0;
    private Vector<Vector<Integer>> _mappings = null;
    private int neighbor_bondnum_A = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
    private int set_bondnum_A = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
    private int neighbor_bondnum_B = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1
    private int set_bondnum_B = 0; //number of remaining molecule B bonds after the clique search, which aren't neighbors
    /*This should be more or equal to all the atom types*/
    String[] SignROW = {"$1", "$2", "$3", "$4", "$5", "$6", "$7", "$8", "$9", "$10", "$11", "$12",
        "$13", "$15", "$16", "$17", "$18", "$19", "$20", "$21", "$22", "$23", "$24",
        "$25", "$26", "$27", "$28", "$29", "$30", "$31", "$32", "$33", "$34", "$35", "$36",
        "$37", "$38", "$39", "$40", "$41", "$42", "$43", "$44", "$45", "$46",
        "$47", "$48", "$49", "$50", "$51", "$52", "$53", "$54", "$55"
    };
    protected boolean new_matrix;
    protected boolean bondTypeFlag = BondType.getInstance().getBondSensitiveFlag();

    /**
     * Creates z new instance of McGregorVFLib
     * @param ac1
     * @param ac2 
     * @param _mappings_org
     */
    public McGregorVFLib(IAtomContainer ac1, IAtomContainer ac2, Vector<Vector<Integer>> _mappings_org) {

        //System.out.println("In McGregorBondTypeInSensitive Class");
//        gvc = GlobalVariableContainer.getInstance();
        this.ac1 = ac1;
        this.ac2 = ac2;
        _mappings = new Vector<Vector<Integer>>(_mappings_org);
        nNum_globalA = 0;
        nNum_globalB = 0;
        bestarcsleft = 0;

        this.globalMCSSize = 0;

        MARCS = new Vector<Integer>();

        BESTARCS = new Stack<Vector<Integer>>();

        //Initialization of global vectors
        i_globalA = new Vector<Integer>();
        i_globalB = new Vector<Integer>();
        c_globalA = new Vector<String>();
        c_globalB = new Vector<String>();

        c_tab1_copy = new Vector<String>();
        c_tab2_copy = new Vector<String>();

        new_matrix = false;



    }

    public void McGregor_IterationStart(Map<Integer, Integer> present_Mapping) throws IOException {

        //protected void McGregor_IterationStart(Vector<Integer> clique_vector) throws IOException {

        this.globalMCSSize = present_Mapping.size();

        c_tab1_copy.clear();
        generate_c_tab1_copy();


        c_tab2_copy.clear();
        generate_c_tab2_copy();


        //find mapped atoms of both molecules and store these in mapped_atoms
        Vector<Integer> mapped_atoms = new Vector<Integer>();
//        System.out.println("\nMapped Atoms");
        for (Map.Entry<Integer, Integer> map : present_Mapping.entrySet()) {
//            System.out.println("i:" + map.getKey() + " j:" + map.getValue());
            mapped_atoms.add(map.getKey());
            mapped_atoms.add(map.getValue());
        }
        int clique_siz = present_Mapping.size();


        Vector<Integer> i_bond_neighborsA = new Vector<Integer>();
        Vector<Integer> i_bond_setA = new Vector<Integer>();
        Vector<String> c_bond_neighborsA = new Vector<String>();
        Vector<String> c_bond_setA = new Vector<String>();

        Vector<Integer> i_bond_neighborsB = new Vector<Integer>();
        Vector<Integer> i_bond_setB = new Vector<Integer>();
        Vector<String> c_bond_neighborsB = new Vector<String>();
        Vector<String> c_bond_setB = new Vector<String>();

        //find unmapped atoms of molecule A
        Vector<Integer> unmapped_atoms_molA = new Vector<Integer>();

        int unmapped_numA = 0;
        boolean atomA_is_unmapped = true;
//
//        System.out.println("\n---------------\n" +
//                "Mapped Atoms: " + mapped_atoms.size());

        for (int a = 0; a < ac1.getAtomCount(); a++) {
            //Atomic list are only numbers from 1 to atom_number1

            for (Integer key : present_Mapping.keySet()) {
                if (key == a) {
                    atomA_is_unmapped = false;
                }
            }


            if (atomA_is_unmapped) {
//                System.out.println("UnMapped Atoms: " + a +
//                        " (" + ac1.getAtom(a).getSymbol() + ")");
                unmapped_atoms_molA.add(unmapped_numA, a);
                unmapped_numA++;
            }
            atomA_is_unmapped = true;
        }


//        System.out.println("neighbor_bondnum_A Before:" + neighbor_bondnum_A);

//        System.out.println("clique Size " + mappingSize);

//        System.out.println("unmapped_atoms_molA: " + unmapped_atoms_molA.size());

        int SR_count = 0;
        processQueryMolecule(
                unmapped_atoms_molA,
                clique_siz,
                i_bond_neighborsA,
                i_bond_setA,
                c_bond_neighborsA,
                c_bond_setA,
                mapped_atoms,
                SR_count);


        //find unmapped atoms of molecule B
        Vector<Integer> unmapped_atoms_molB = new Vector<Integer>();
        int unmapped_numB = 0;
        boolean atomB_is_unmapped = true;

//        System.out.println("neighbor_bondnum_A After:" + neighbor_bondnum_A);
//
//        System.out.println("\n---------------\n");
        for (int a = 0; a < ac2.getAtomCount(); a++) {
            for (Integer value : present_Mapping.values()) {

                if (a == value) {
                    atomB_is_unmapped = false;
                }
            }
            if (atomB_is_unmapped) {
//                System.out.println("UnMapped Atoms: " + a +
//                        " (" + ac2.getAtom(a).getSymbol() + ")");
                unmapped_atoms_molB.add(unmapped_numB, a);
                unmapped_numB++;
            }
            atomB_is_unmapped = true;
        }

//        System.out.println("unmapped_atoms_molB: " + unmapped_atoms_molB.size());

        //Extract bonds which are related with unmapped atoms of molecule B.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A

        processTargetMolecule(
                c_bond_neighborsA,
                unmapped_atoms_molB,
                clique_siz,
                i_bond_neighborsA,
                i_bond_neighborsB,
                i_bond_setB,
                c_bond_neighborsB,
                c_bond_setB,
                mapped_atoms,
                SR_count);


        boolean dummy = false;
//
//        System.out.println("neighbor_bondnum_A: " + neighbor_bondnum_A);
//        System.out.println("neighbor_bondnum_B: " + neighbor_bondnum_B);
//        System.out.println("set_bondnum_A: " + set_bondnum_A);
//        System.out.println("set_bondnum_B: " + set_bondnum_B);


        Iterator(dummy, present_Mapping.size(), mapped_atoms, neighbor_bondnum_A, neighbor_bondnum_B, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, set_bondnum_A, set_bondnum_B, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);

        //System.exit(1); //uncomment to debug


    }

    private void processQueryMolecule(
            Vector<Integer> unmapped_atoms_molA,
            int mappingSize,
            Vector<Integer> i_bond_neighborsA,
            Vector<Integer> i_bond_setA,
            Vector<String> c_bond_neighborsA,
            Vector<String> c_bond_setA,
            Vector<Integer> mapped_atoms,
            int SR_count) {


//        int mappingSize = clique_vector.size();

        boolean bond_considered = false;
        boolean normal_bond = true;

//        System.out.println("\n" + c_tab1_copy + "\n");

        IAtomContainer reactant = ac1;

        for (int a = 0; a < reactant.getBondCount(); a++) {


            Integer indexI = reactant.getAtomNumber(reactant.getBond(a).getAtom(0));
            Integer indexJ = reactant.getAtomNumber(reactant.getBond(a).getAtom(1));
            String AtomI = reactant.getBond(a).getAtom(0).getSymbol();
            String AtomJ = reactant.getBond(a).getAtom(1).getSymbol();
            Integer order = reactant.getBond(a).getOrder().ordinal() + 1;

//            System.out.println(AtomI + "= , =" + AtomJ );
            for (Integer unMappedAtomIndex : unmapped_atoms_molA) {



                if (unMappedAtomIndex.equals(indexI)) {
//                    System.out.println("Unmapped Atom is equal to Reaction Bond table");
//                    System.out.println("unMappedAtomIndex " + unMappedAtomIndex + " " + "indexI: " + indexI);
//                    System.out.println("unMappedAtomIndex=IndexI " + reactant.getAtom(unMappedAtomIndex).getSymbol());
                    for (int c = 0; c < mappingSize; c++) {

//                        System.out.println("\n*****\nmapped_atoms.elementAt(c * 2): " + mapped_atoms.elementAt(c * 2));
//                        System.out.println("indexJ: " + indexJ);

                        if (mapped_atoms.elementAt(c * 2).equals(indexJ)) {

//                            System.out.println("i_bond_neighbor_atoms_A: ");
//                            System.out.println(indexI +
//                                    " " + indexJ +
//                                    " " + order);
//

                            i_bond_neighborsA.add(indexI);
                            i_bond_neighborsA.add(indexJ);
                            i_bond_neighborsA.add(order);

                            if (c_tab1_copy.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0) {


                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(SignROW[SR_count]);
                                c_bond_neighborsA.add("X");
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));

                                change_char_bonds(indexJ, SignROW[SR_count],
                                        reactant.getBondCount(), ac1, c_tab1_copy);

                                int cor_atom = search_corresponding_atom(mappingSize, indexJ, 1, mapped_atoms);
                                change_char_bonds(cor_atom, SignROW[SR_count], ac2.getBondCount(), ac2, c_tab2_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 2));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighbor_bondnum_A++;
                        }
                    }
                    if (normal_bond) {
                        i_bond_setA.add(indexI);
                        i_bond_setA.add(indexJ);
                        i_bond_setA.add(order);
                        c_bond_setA.add(c_tab1_copy.get(a * 4 + 0));
                        c_bond_setA.add(c_tab1_copy.get(a * 4 + 1));
                        c_bond_setA.add("X");
                        c_bond_setA.add("X");
                        set_bondnum_A++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                //Does a ungemaptes atom at second position in the connection occur?
                if (unMappedAtomIndex.equals(indexJ)) {
                    for (int c = 0; c < mappingSize; c++) {


                        if (mapped_atoms.elementAt(c * 2 + 0).equals(indexI)) {



                            i_bond_neighborsA.add(indexI);
                            i_bond_neighborsA.add(indexJ);
                            i_bond_neighborsA.add(order);

                            if (c_tab1_copy.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0) {


                                c_bond_neighborsA.add(SignROW[SR_count]);
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add("X");
                                change_char_bonds(indexI, SignROW[SR_count],
                                        ac1.getBondCount(), ac1, c_tab1_copy);

                                int cor_atom = search_corresponding_atom(mappingSize, indexI, 1, mapped_atoms);
                                change_char_bonds(cor_atom, SignROW[SR_count], ac2.getBondCount(), ac2, c_tab2_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 2));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighbor_bondnum_A++;
                            //System.out.println("Neighbor");
                            //System.out.println(neighbor_bondnum_A);
                        }
                    }
                    if (normal_bond) {

//                        System.out.println("\n\n normal_bond:" + normal_bond);


                        i_bond_setA.add(indexI);
                        i_bond_setA.add(indexJ);
                        i_bond_setA.add(order);
                        c_bond_setA.add(AtomI);
                        c_bond_setA.add(AtomJ);
                        c_bond_setA.add("X");
                        c_bond_setA.add("X");
                        set_bondnum_A++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (bond_considered) {
                    break;
                }
            }
            bond_considered = false;
        }

        /*******************************************************************************///        System.out.println("Neighbor A & B");
    }

    private void processTargetMolecule(
            Vector<String> c_bond_neighborsA,
            Vector<Integer> unmapped_atoms_molB,
            int mappingSize,
            Vector<Integer> i_bond_neighborsA,
            Vector<Integer> i_bond_neighborsB,
            Vector<Integer> i_bond_setB,
            Vector<String> c_bond_neighborsB,
            Vector<String> c_bond_setB,
            Vector<Integer> mapped_atoms,
            int SR_count) {

        int unmapped_numB = unmapped_atoms_molB.size();
        boolean bond_considered = false;
        boolean normal_bond = true;

        IAtomContainer product = ac2;

        for (int a = 0; a < product.getBondCount(); a++) {


            Integer indexI = product.getAtomNumber(product.getBond(a).getAtom(0));
            Integer indexJ = product.getAtomNumber(product.getBond(a).getAtom(1));
            String AtomI = product.getBond(a).getAtom(0).getSymbol();
            String AtomJ = product.getBond(a).getAtom(1).getSymbol();
            Integer order = product.getBond(a).getOrder().ordinal() + 1;

            for (int b = 0; b < unmapped_numB; b++) {
                if (unmapped_atoms_molB.elementAt(b).equals(indexI)) {
                    for (int c = 0; c < mappingSize; c++) {
                        if (mapped_atoms.elementAt(c * 2 + 1).equals(indexJ)) {
                            i_bond_neighborsB.add(indexI);
                            i_bond_neighborsB.add(indexJ);
                            i_bond_neighborsB.add(order);
                            if (c_tab2_copy.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(SignROW[SR_count]);
                                c_bond_neighborsB.add("X");
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                change_char_bonds(indexJ, SignROW[SR_count], ac2.getBondCount(), ac2, c_tab2_copy);
                                int cor_atom = search_corresponding_atom(mappingSize, indexJ, 2, mapped_atoms);
                                //Commented by Asad
//                                change_char_bonds(cor_atom, SignROW[SR_count], neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA);
                                change_char_bonds(cor_atom, SignROW[SR_count], ac1.getBondCount(), ac1, c_tab1_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add("X");
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighbor_bondnum_B++;
                        }
                    }
                    if (normal_bond) {
                        i_bond_setB.add(indexI);
                        i_bond_setB.add(indexJ);
                        i_bond_setB.add(order);
                        c_bond_setB.add(c_tab2_copy.get(a * 4 + 0));
                        c_bond_setB.add(c_tab2_copy.get(a * 4 + 1));
                        c_bond_setB.add("X");
                        c_bond_setB.add("X");
                        set_bondnum_B++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (unmapped_atoms_molB.elementAt(b) == indexJ) {
                    for (int c = 0; c < mappingSize; c++) {
                        if (mapped_atoms.elementAt(c * 2 + 1).equals(indexI)) {
                            i_bond_neighborsB.add(indexI);
                            i_bond_neighborsB.add(indexJ);
                            i_bond_neighborsB.add(order);
                            if (c_tab2_copy.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                                c_bond_neighborsB.add(SignROW[SR_count]);
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add("X");
                                change_char_bonds(indexI, SignROW[SR_count], ac2.getBondCount(), ac2, c_tab2_copy);
                                int cor_atom = search_corresponding_atom(mappingSize, indexI, 2, mapped_atoms);
//                                change_char_bonds(cor_atom, SignROW[SR_count], neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA);
                                change_char_bonds(cor_atom, SignROW[SR_count], ac1.getBondCount(), ac1, c_tab1_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 2));
                                c_bond_neighborsB.add("X");
                            }
                            normal_bond = false;
                            neighbor_bondnum_B++;
                        }
                    }
                    if (normal_bond) {

//                        System.out.println("\n\n normal_bond:" + normal_bond);

                        i_bond_setB.add(indexI);
                        i_bond_setB.add(indexJ);
                        i_bond_setB.add(order);
                        c_bond_setB.add(AtomI);
                        c_bond_setB.add(AtomJ);
                        c_bond_setB.add("X");
                        c_bond_setB.add("X");
                        set_bondnum_B++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (bond_considered) {
                    break;
                }
            }
            bond_considered = false;
        }

    }

    private int Iterator(boolean MAPPING_check,
            int mapped_atom_num,
            Vector<Integer> mapped_atoms_org,
            int neighbor_bondnum_A,
            int neighbor_bondnum_B,
            Vector<Integer> i_bond_neighbor_atoms_A,
            Vector<Integer> i_bond_neighbor_atoms_B,
            Vector<String> c_bond_neighborsA,
            Vector<String> c_bond_neighborsB,
            int set_num_A, int set_num_B,
            Vector<Integer> i_bond_setA,
            Vector<Integer> i_bond_setB,
            Vector<String> c_bond_setA,
            Vector<String> c_bond_setB) throws IOException {


        Vector<Integer> mapped_atoms = new Vector<Integer>(mapped_atoms_org);

        //check possible mappings:
        boolean no_further_mapping_possible = true;

        for (int row = 0; row < neighbor_bondnum_A; row++) {

            for (int column = 0; column < neighbor_bondnum_B; column++) {
                String G1A = c_bond_neighborsA.get(row * 4 + 0);
                String G2A = c_bond_neighborsA.get(row * 4 + 1);
                String G1B = c_bond_neighborsB.get(column * 4 + 0);
                String G2B = c_bond_neighborsB.get(column * 4 + 1);

                //Get the atom index from the i_bond neighbor vactor

                int Index_I = i_bond_neighbor_atoms_A.get(row * 3 + 0);
                int Index_IPlus1 = i_bond_neighbor_atoms_A.get(row * 3 + 1);

                //Get the atoms
                IAtom R1_A = ac1.getAtom(Index_I);
                IAtom R2_A = ac1.getAtom(Index_IPlus1);
                IBond ReactantBond = ac1.getBond(R1_A, R2_A);

                //Get the atom index from the i_bond neighbor vactor
                int Index_J = i_bond_neighbor_atoms_B.get(column * 3 + 0);
                int Index_JPlus1 = i_bond_neighbor_atoms_B.get(column * 3 + 1);

                //Get the atoms
                IAtom P1_B = ac2.getAtom(Index_J);
                IAtom P2_B = ac2.getAtom(Index_JPlus1);
                IBond ProductBond = ac2.getBond(P1_B, P2_B);


                if (bondTypeFlag && bondMatch(ReactantBond, ProductBond)) {
                    if ((G1A.compareToIgnoreCase(G1B) == 0 && G2A.compareToIgnoreCase(G2B) == 0) || (G1A.compareToIgnoreCase(G2B) == 0 && G2A.compareToIgnoreCase(G1B) == 0)) {

                        no_further_mapping_possible = false;
                    }
//                      
                } else if (!bondTypeFlag) {
//                  System.out.println("Eins bei: " + G1A + " " + G2A + " " + G1B + " " + G2B);

                    if (((G1A.compareToIgnoreCase(G1B) == 0) && (G2A.compareToIgnoreCase(G2B) == 0)) || ((G1A.compareToIgnoreCase(G2B) == 0) && (G2A.compareToIgnoreCase(G1B) == 0))) {

//                    System.out.println("Eins bei: " + G1A + " " + G2A + " " + G1B + " " + G2B);
                        no_further_mapping_possible = false;
                    }


                }


            }
        }

        if (neighbor_bondnum_A == 0 || neighbor_bondnum_B == 0 || MAPPING_check || no_further_mapping_possible) {
//            System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
            try {
                if (mapped_atom_num >= globalMCSSize) {
//                    System.out.println("Hello-1");
                    if (mapped_atom_num > globalMCSSize) {
//                        System.out.println("Hello-2");
                        this.globalMCSSize = mapped_atom_num;
//                        System.out.println("best_MAPPING_size: " + gvc.getMappingSize());
                        _mappings.clear();
                        //_mappings=new Vector<Vector<Integer>>();
                    }
                    _mappings.addElement(mapped_atoms);
                }




            } catch (Exception ex) {
                ex.printStackTrace();
            }
            //System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
            //NoMoreMappingFlag = false;
            return 0;
        }

        i_globalA.clear();
        i_globalB.clear();
        c_globalA.clear();
        c_globalB.clear();

        //redefining of global vectors and variables
        nNum_globalA = neighbor_bondnum_A; //N global variable defined

        nNum_globalB = neighbor_bondnum_B; //N global variable defined

        i_globalA.addAll(i_bond_neighbor_atoms_A);
        i_globalB.addAll(i_bond_neighbor_atoms_B);
        c_globalA.addAll(c_bond_neighborsA);
        c_globalB.addAll(c_bond_neighborsB);
        MARCS.clear();

        MARCS.setSize(neighbor_bondnum_A * neighbor_bondnum_B);
        for (int i = 0; i < neighbor_bondnum_A * neighbor_bondnum_B; i++) {

            MARCS.setElementAt(0, i);
        }

        for (int row = 0; row < neighbor_bondnum_A; row++) {
            for (int column = 0; column < neighbor_bondnum_B; column++) {

                String G1A = c_bond_neighborsA.get(row * 4 + 0);
                String G2A = c_bond_neighborsA.get(row * 4 + 1);
                String G1B = c_bond_neighborsB.get(column * 4 + 0);
                String G2B = c_bond_neighborsB.get(column * 4 + 1);

                int Index_I = i_bond_neighbor_atoms_A.get(row * 3 + 0);
                int Index_IPlus1 = i_bond_neighbor_atoms_A.get(row * 3 + 1);

                IAtom R1_A = ac1.getAtom(Index_I);
                IAtom R2_A = ac1.getAtom(Index_IPlus1);
                IBond ReactantBond = ac1.getBond(R1_A, R2_A);

                int Index_J = i_bond_neighbor_atoms_B.get(column * 3 + 0);
                int Index_JPlus1 = i_bond_neighbor_atoms_B.get(column * 3 + 1);

                IAtom P1_B = ac2.getAtom(Index_J);
                IAtom P2_B = ac2.getAtom(Index_JPlus1);
                IBond ProductBond = ac2.getBond(P1_B, P2_B);

                if ((G1A.compareToIgnoreCase(G1B) == 0 && G2A.compareToIgnoreCase(G2B) == 0) || (G1A.compareToIgnoreCase(G2B) == 0 && G2A.compareToIgnoreCase(G1B) == 0)) {


                    if (bondTypeFlag && bondMatch(ReactantBond, ProductBond)) {
                        MARCS.setElementAt(1, row * neighbor_bondnum_B + column);

                    } else if (!bondTypeFlag) {
                        MARCS.setElementAt(1, row * neighbor_bondnum_B + column);
                    }
                }


            }
        }
        first = last = new BinaryTree();
        last.setValue(-1);
        last.equal = null;
        last.not_equal = null;

        bestarcsleft = 0;

        startsearch();
        Stack<Vector<Integer>> BESTARCS_copy = new Stack<Vector<Integer>>();


        BESTARCS_copy.addAll(BESTARCS);
        while (!BESTARCS.empty()) {
            BESTARCS.pop();
        }

        while (!BESTARCS_copy.empty()) {

            Vector<Integer> MARCS_vector = new Vector<Integer>(BESTARCS_copy.peek());
            Vector<Integer> new_MAPPING = find_mcgregor_MAPPING(MARCS_vector, mapped_atom_num, mapped_atoms, neighbor_bondnum_A, i_bond_neighbor_atoms_A, neighbor_bondnum_B, i_bond_neighbor_atoms_B);

            int new_MAPPING_size = (new_MAPPING.size() / 2);
            boolean no_further_MAPPINGS = false;
            if (mapped_atom_num == new_MAPPING_size) {
                no_further_MAPPINGS = true;
            }


            int new_neighbor_numA = 0; //instead of neighbor_bondnum_A

            int new_neighbor_numB = 0; //instead of neighbor_bondnum_B

            Vector<Integer> new_i_neighborsA = new Vector<Integer>(); //instead of i_bond_neighbor_atoms_A

            Vector<Integer> new_i_neighborsB = new Vector<Integer>(); //instead of i_bond_neighbor_atoms_B

            Vector<String> new_c_neighborsA = new Vector<String>(); //instead of c_bond_neighborsA

            Vector<String> new_c_neighborsB = new Vector<String>(); //instead of c_bond_neighborsB

            new_i_neighborsA.clear();
            new_i_neighborsB.clear();
            new_c_neighborsA.clear();
            new_c_neighborsB.clear();


            //new values for set_num_A + set_num_B
            //new arrays for i_bond_setA + i_bond_setB + c_bond_setB + c_bond_setB

            set_bondnum_A = 0; //instead of set_num_A

            set_bondnum_B = 0; //instead of set_num_B

            Vector<Integer> new_i_bond_setA = new Vector<Integer>(); //instead of i_bond_setA

            Vector<Integer> new_i_bond_setB = new Vector<Integer>(); //instead of i_bond_setB

            Vector<String> new_c_bond_setA = new Vector<String>(); //instead of c_bond_setA

            Vector<String> new_c_bond_setB = new Vector<String>(); //instead of c_bond_setB
            Vector<String> c_setB_copy = new Vector<String>();
            Vector<String> c_setA_copy = new Vector<String>();
//            c_setB_copy.clear();
            generate_c_setA_copy(set_num_A, c_bond_setA, c_setA_copy);
            generate_c_setB_copy(set_num_B, c_bond_setB, c_setB_copy);
//            Vector<String> c_setA_copy = new Vector<String>(c_bond_setA);

            //find unmapped atoms of molecule A
            Vector<Integer> unmapped_atoms_molA = new Vector<Integer>();
            int unmapped_numA = 0;
            boolean atomA_is_unmapped = true;

            for (int a = 0; a < ac1.getAtomCount(); a++) {
                for (int b = 0; b < new_MAPPING_size; b++) {
                    if (a == new_MAPPING.elementAt(b * 2 + 0)) {
                        atomA_is_unmapped = false;
                    }

                }
                if (atomA_is_unmapped) {
                    unmapped_atoms_molA.add(a);
                    unmapped_numA++;

                }


                atomA_is_unmapped = true;
            }


            //The special signs must be transfered to the corresponding atoms of molecule B

            int SR_count = 0;
            boolean bond_considered = false;
            boolean normal_bond = true;
            for (int a = 0; a < set_num_A; a++) {

                int _elementAt_a = i_bond_setA.get(a * 3 + 0);
                for (int b = 0; b < unmapped_numA; b++) {
                    Integer unMappedAtomIndex = unmapped_atoms_molA.elementAt(b);
                    if (unMappedAtomIndex == _elementAt_a) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.elementAt(c * 2 + 0) == (i_bond_setA.elementAt(a * 3 + 1))) {

                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
                                new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                if (c_setA_copy.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0) {

                                    new_c_neighborsA.add(SignROW[SR_count]);
                                    new_c_neighborsA.add("X");
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    change_char_bonds(i_bond_setA.get(a * 3 + 1), SignROW[SR_count], set_num_A, i_bond_setA, c_setA_copy);
                                    int cor_atom = search_corresponding_atom(new_MAPPING_size, i_bond_setA.get(a * 3 + 1), 1, new_MAPPING);
                                    change_char_bonds(cor_atom, SignROW[SR_count], set_num_B, i_bond_setB, c_setB_copy);
                                    SR_count++;

                                } else {

                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add("X");
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 3));

                                }

                                normal_bond = false;
                                new_neighbor_numA++;

                            }
                        }

                        if (normal_bond) {

                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 0));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 1));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 2));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 0));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 1));
                            new_c_bond_setA.add("X");
                            new_c_bond_setA.add("X");
                            set_bondnum_A++;

                        }
                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (unMappedAtomIndex == i_bond_setA.elementAt(a * 3 + 1)) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.elementAt(c * 2 + 0) == i_bond_setA.elementAt(a * 3 + 0)) {

                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
                                if (c_setA_copy.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0) {

                                    new_c_neighborsA.add(SignROW[SR_count]);
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                    new_c_neighborsA.add("X");
                                    change_char_bonds(i_bond_setA.get(a * 3 + 0), SignROW[SR_count], set_num_A, i_bond_setA, c_setA_copy);
                                    int cor_atom = search_corresponding_atom(new_MAPPING_size, i_bond_setA.get(a * 3 + 0), 1, new_MAPPING);
                                    change_char_bonds(cor_atom, SignROW[SR_count], set_num_B, i_bond_setB, c_setB_copy);
                                    SR_count++;

                                } else {

                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 2));
                                    new_c_neighborsA.add("X");

                                }

                                normal_bond = false;
                                new_neighbor_numA++;

                            }


                        }
                        if (normal_bond) {

                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 0));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 1));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 2));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 0));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 1));
                            new_c_bond_setA.add("X");
                            new_c_bond_setA.add("X");
                            set_bondnum_A++;

                        }

                        normal_bond = true;
                        bond_considered = true;
                    }

                    if (bond_considered) {
                        break;
                    }

                }

                bond_considered = false;
            }

            //find unmapped atoms of molecule B

            Vector<Integer> unmapped_atoms_molB = new Vector<Integer>();
            int unmapped_numB = 0;
            boolean atomB_is_unmapped = true;

            for (int a = 0; a < ac2.getAtomCount(); a++) {
                for (int b = 0; b < new_MAPPING_size; b++) {
                    if (a == new_MAPPING.elementAt(b * 2 + 1)) {
                        atomB_is_unmapped = false;
                    }

                }
                if (atomB_is_unmapped) {
                    unmapped_atoms_molB.add(a);
                    unmapped_numB++;

                }


                atomB_is_unmapped = true;
            }

            //The special signs must be transfered to the corresponding atoms of molecule A

            bond_considered = false;
            normal_bond = true;
            for (int a = 0; a < set_num_B; a++) {
                for (int b = 0; b < unmapped_numB; b++) {
                    if (unmapped_atoms_molB.elementAt(b) == (i_bond_setB.elementAt(a * 3 + 0))) {
                        for (int c = 0; c < new_MAPPING_size; c++) {
                            if (new_MAPPING.elementAt(c * 2 + 1) == i_bond_setB.elementAt(a * 3 + 1)) {
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 0));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 1));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 2));
                                if (c_setB_copy.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(SignROW[SR_count]);
                                    new_c_neighborsB.add("X");
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    change_char_bonds(i_bond_setB.get(a * 3 + 1), SignROW[SR_count], set_num_B, i_bond_setB, c_setB_copy);
                                    int cor_atom = search_corresponding_atom(new_MAPPING_size, i_bond_setB.get(a * 3 + 1), 2, new_MAPPING);
                                    change_char_bonds(cor_atom, SignROW[SR_count], new_neighbor_numA, new_i_neighborsA, new_c_neighborsA);
                                    SR_count++;

                                } else {
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add("X");
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 3));
                                }

                                normal_bond = false;
                                new_neighbor_numB++;

                            }


                        }
                        if (normal_bond) {
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 0));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 1));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 2));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 0));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 1));
                            new_c_bond_setB.add("X");
                            new_c_bond_setB.add("X");
                            set_bondnum_B++;

                        }


                        normal_bond = true;
                        bond_considered =
                                true;
                    }
                    if (unmapped_atoms_molB.elementAt(b) == i_bond_setB.elementAt(a * 3 + 1)) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.elementAt(c * 2 + 1) == i_bond_setB.elementAt(a * 3 + 0)) {

                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 0));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 1));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 2));

                                if (c_setB_copy.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                                    new_c_neighborsB.add(SignROW[SR_count]);
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add("X");
                                    change_char_bonds(i_bond_setB.get(a * 3 + 0), SignROW[SR_count], set_num_B, i_bond_setB, c_setB_copy);
                                    int cor_atom = search_corresponding_atom(new_MAPPING_size, i_bond_setB.get(a * 3 + 0), 2, new_MAPPING);
                                    change_char_bonds(cor_atom, SignROW[SR_count], new_neighbor_numA, new_i_neighborsA, new_c_neighborsA);
                                    SR_count++;

                                } else {

                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 2));
                                    new_c_neighborsB.add("X");
                                }

                                normal_bond = false;
                                new_neighbor_numB++;

                            }


                        }

                        if (normal_bond) {

                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 0));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 1));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 2));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 0));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 1));
                            new_c_bond_setB.add("X");
                            new_c_bond_setB.add("X");
                            set_bondnum_B++;

                        }


                        normal_bond = true;
                        bond_considered =
                                true;
                    }

                    if (bond_considered) {
                        break;
                    }

                }
                bond_considered = false;
            }


            // System.out.println("Mapped Atoms before Iterator2: " + mapped_atoms);
            Iterator(no_further_MAPPINGS, new_MAPPING_size, new_MAPPING, new_neighbor_numA, new_neighbor_numB, new_i_neighborsA, new_i_neighborsB, new_c_neighborsA, new_c_neighborsB,
                    set_bondnum_A, set_bondnum_B, new_i_bond_setA, new_i_bond_setB, new_c_bond_setA, new_c_bond_setB);
            BESTARCS_copy.pop();
            //System.out.println("Schleife beendet in Iterator!!!!");
        }

        //}
        //System.out.println("In the Iterator Termination");
        //System.out.println("============+++++++++==============");

        //System.out.println("Mapped Atoms before Iterator Over: " + mapped_atoms);
        return 0;
    }

    private int generate_c_tab1_copy() throws IOException {
        IAtomContainer reactant = ac1;
        for (int a = 0; a < reactant.getBondCount(); a++) {
            String AtomI = reactant.getBond(a).getAtom(0).getSymbol();
            String AtomJ = reactant.getBond(a).getAtom(1).getSymbol();
            c_tab1_copy.addElement(AtomI);
            c_tab1_copy.addElement(AtomJ);
            c_tab1_copy.addElement("X");
            c_tab1_copy.addElement("X");
        }
        return 0;
    }

    private int generate_c_tab2_copy() throws IOException {
        IAtomContainer product = ac2;
        for (int a = 0; a < product.getBondCount(); a++) {
            String AtomI = product.getBond(a).getAtom(0).getSymbol();
            String AtomJ = product.getBond(a).getAtom(1).getSymbol();
            c_tab2_copy.addElement(AtomI);
            c_tab2_copy.addElement(AtomJ);
            c_tab2_copy.addElement("X");
            c_tab2_copy.addElement("X");
        }

        return 0;

    }

    private int change_char_bonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, Vector<Integer> i_bond_neighbors, Vector<String> c_bond_neighbors) {

        for (int a = 0; a < neighbor_bondnum; a++) {
            if ((i_bond_neighbors.elementAt(a * 3 + 0) == (corresponding_atom)) && (c_bond_neighbors.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 2, c_bond_neighbors.get(a * 4 + 0));
                c_bond_neighbors.set(a * 4 + 0, new_symbol);
            }

            if ((i_bond_neighbors.elementAt(a * 3 + 1) == (corresponding_atom)) && (c_bond_neighbors.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 3, c_bond_neighbors.get(a * 4 + 1));
                c_bond_neighbors.set(a * 4 + 1, new_symbol);
            }

        }

        return 0;
    }

    private int change_char_bonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, IAtomContainer ac, Vector<String> c_bond_neighbors) {
        //private int change_char_bonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, Vector<Integer> i_bond_neighbors, Vector<String> c_bond_neighbors) {

        for (int a = 0; a < neighbor_bondnum; a++) {
            IBond bond = ac.getBond(a);
            if ((ac.getAtomNumber(bond.getAtom(0)) == corresponding_atom) && (c_bond_neighbors.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 2, c_bond_neighbors.get(a * 4 + 0));
                c_bond_neighbors.set(a * 4 + 0, new_symbol);
            }

            if ((ac.getAtomNumber(bond.getAtom(1)) == corresponding_atom) && (c_bond_neighbors.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 3, c_bond_neighbors.get(a * 4 + 1));
                c_bond_neighbors.set(a * 4 + 1, new_symbol);
            }

        }
//        for (int a = 0; a < neighbor_bondnum; a++) {
//            if ((i_bond_neighbors.elementAt(a * 3 + 0) == (corresponding_atom)) && (c_bond_neighbors.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0)) {
//                c_bond_neighbors.set(a * 4 + 2, c_bond_neighbors.get(a * 4 + 0));
//                c_bond_neighbors.set(a * 4 + 0, new_symbol);
//            }
//
//            if ((i_bond_neighbors.elementAt(a * 3 + 1) == (corresponding_atom)) && (c_bond_neighbors.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0)) {
//                c_bond_neighbors.set(a * 4 + 3, c_bond_neighbors.get(a * 4 + 1));
//                c_bond_neighbors.set(a * 4 + 1, new_symbol);
//            }
//
//        }

        return 0;
    }

    public boolean bondMatch(IBond ReactantBond, IBond ProductBond) {
        boolean Flag = false;
        int ReactantBondType = ReactantBond.getOrder().ordinal();
        int ProductBondType = ProductBond.getOrder().ordinal();
        if (bondTypeFlag) {
            if ((ReactantBond.getFlag(CDKConstants.ISAROMATIC) == ProductBond.getFlag(CDKConstants.ISAROMATIC)) && (ReactantBondType == ProductBondType)) {
                Flag = true;
            } else if (ReactantBond.getFlag(CDKConstants.ISAROMATIC) && ProductBond.getFlag(CDKConstants.ISAROMATIC)) {
                Flag = true;
            }
        }
        return Flag;
    }

    private Vector<Integer> find_mcgregor_MAPPING(Vector<Integer> MARCS_vector_org, int mapped_atoms_num, Vector<Integer> current_MAPPING_org, int bondnum_A, Vector<Integer> i_bonds_A_org, int bondnum_B, Vector<Integer> i_bonds_B_org) {

        Vector<Integer> MARCS_vector = new Vector<Integer>(MARCS_vector_org);
        Vector<Integer> current_MAPPING = new Vector<Integer>(current_MAPPING_org);
        Vector<Integer> i_bonds_A = new Vector<Integer>(i_bonds_A_org);
        Vector<Integer> i_bonds_B = new Vector<Integer>(i_bonds_B_org);

        Vector<Integer> additional_mapping = new Vector<Integer>();



        for (int x = 0; x < bondnum_A; x++) {
            for (int y = 0; y < bondnum_B; y++) {


                if (MARCS_vector.elementAt(x * bondnum_B + y) == 1) {


                    int Atom1_moleculeA = i_bonds_A.get(x * 3 + 0);
                    int Atom2_moleculeA = i_bonds_A.get(x * 3 + 1);
                    int Atom1_moleculeB = i_bonds_B.get(y * 3 + 0);
                    int Atom2_moleculeB = i_bonds_B.get(y * 3 + 1);

                    IAtom R1_A = ac1.getAtom(Atom1_moleculeA);
                    IAtom R2_A = ac1.getAtom(Atom2_moleculeA);
                    IBond ReactantBond = ac1.getBond(R1_A, R2_A);

                    IAtom P1_B = ac2.getAtom(Atom1_moleculeB);
                    IAtom P2_B = ac2.getAtom(Atom2_moleculeB);
                    IBond ProductBond = ac2.getBond(P1_B, P2_B);


//                    //Bond Order Check Introduced by Asad

                    boolean bMatch = bondMatch(ReactantBond, ProductBond);

                    for (int z = 0; z < mapped_atoms_num; z++) {

                        int Mapped_Atom_1 = current_MAPPING.elementAt(z * 2 + 0);
                        int Mapped_Atom_2 = current_MAPPING.elementAt(z * 2 + 1);

                        if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
                            if (bondTypeFlag && bMatch) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            } else if (!bondTypeFlag) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            }
                        }

                        if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
                            if (bondTypeFlag && bMatch) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            } else if (!bondTypeFlag) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            }
                        }

                        if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
                            if (bondTypeFlag && bMatch) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            } else if (!bondTypeFlag) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            }
                        }

                        if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
                            if (bondTypeFlag && bMatch) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            } else if (!bondTypeFlag) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            }
                        }

                    }
                }
            }
        }


        int additional_mapping_size = additional_mapping.size();

        //add McGregorBondTypeInSensitive mapping to the Clique mapping
        for (int a = 0; a < additional_mapping_size; a = a + 2) {
            current_MAPPING.add(additional_mapping.get(a + 0));
            current_MAPPING.add(additional_mapping.get(a + 1));
        }
        //remove recurring mappings from current_MAPPING

        Vector<Integer> unique_MAPPING = remove_recurring_mappings(current_MAPPING);

        return unique_MAPPING;
    }

    //Function compaires z structure array with itself. Sometimes z mapping occurs several times within the array.
//The function eliminates these recurring mappings. Function is called in function best_solution.
//The function is called by itself as long as the last list element is processed.
    private Vector<Integer> remove_recurring_mappings(Vector<Integer> atom_mapping) {

        //System.out.println("Mapped Atoms remove_recurring_mappings: " + atom_mapping);

        //Vector<Integer> atom_mapping = new Vector<Integer>(atom_mapping_org);

        boolean exist = true;
        Vector<Integer> temp_map = new Vector<Integer>();
        int temp_counter = 0;
        int atom_mapping_size = atom_mapping.size();
        for (int x = 0; x < atom_mapping_size; x = x + 2) {
            int atom = atom_mapping.get(x);
            for (int y = x + 2; y < atom_mapping_size; y = y + 2) {
                if (atom == atom_mapping.elementAt(y)) {
                    exist = false;
                }

            }
            if (exist == true) {
                temp_map.add(atom_mapping.get(x));
                temp_map.add(atom_mapping.get(x + 1));
                temp_counter = temp_counter + 2;
            }

            exist = true;
        }


        return temp_map;
    }

    private void partsearch(int xstart, int ystart, Vector<Integer> TEMPMARCS_ORG) {


        int x = xstart;
        int y = ystart;

        Vector<Integer> TEMPMARCS = new Vector<Integer>(TEMPMARCS_ORG);

        if (TEMPMARCS.elementAt(xstart * nNum_globalB + ystart) == (1)) {

            remove_redundant_arcs(xstart, ystart, TEMPMARCS);
            int arcsleft = 0;

            for (int a = 0; a < nNum_globalA; a++) {
                for (int b = 0; b < nNum_globalB; b++) {

                    if (TEMPMARCS.elementAt(a * nNum_globalB + b) == (1)) {
                        arcsleft++;
                    }

                }
            }

            //test Bestarcsleft and skip rest if needed
            if (arcsleft >= bestarcsleft) {
                do {
                    y++;
                    if (y == nNum_globalB) {
                        y = 0;
                        x++;

                    }



                } while ((x < nNum_globalA) && (TEMPMARCS.elementAt(x * nNum_globalB + y) != 1)); //Correction by ASAD set value minus 1



                if (x < nNum_globalA) {

                    partsearch(x, y, TEMPMARCS);
                    TEMPMARCS.setElementAt(0, x * nNum_globalB + y);
                    partsearch(x, y, TEMPMARCS);

                } else {
                    if (arcsleft > bestarcsleft) {
                        remove_tree_structure(first);
                        first = last = new BinaryTree();
                        last.setValue(-1);
                        last.equal = null;
                        last.not_equal = null;

                        while (!BESTARCS.empty()) {
                            BESTARCS.pop();
                        }

                    }
                    bestarcsleft = arcsleft;

                    if (check_MARCS(TEMPMARCS)) {
                        BESTARCS.push(TEMPMARCS);
                    }

                }
            }
        } else {
            do {
                y++;
                if (y == nNum_globalB) {
                    y = 0;
                    x++;

                }


            } while ((x < nNum_globalA) && (TEMPMARCS.elementAt(x * nNum_globalB + y) != 1)); //Correction by ASAD set value minus 1

            if (x < nNum_globalA) {

                partsearch(x, y, TEMPMARCS);
                TEMPMARCS.setElementAt(0, x * nNum_globalB + y);
                partsearch(x, y, TEMPMARCS);
            } else {
                int arcsleft = 0;
                for (int a = 0; a < nNum_globalA; a++) {
                    for (int b = 0; b < nNum_globalB; b++) {
                        if (TEMPMARCS.elementAt(a * nNum_globalB + b) == (1)) {
                            arcsleft++;
                        }

                    }
                }
                if (arcsleft >= bestarcsleft) {
                    if (arcsleft > bestarcsleft) {
                        remove_tree_structure(first);
                        first = last = new BinaryTree();
                        last.setValue(-1);
                        last.equal = null;
                        last.not_equal = null;
                        while (!BESTARCS.empty()) {
                            BESTARCS.pop();
                        }

                    }
                    bestarcsleft = arcsleft;

                    if (check_MARCS(TEMPMARCS)) {
                        BESTARCS.push(TEMPMARCS);
                    }

                }
            }
        }
    }

    //The function is called in function partsearch. The function is given z temporary matrix and z position (row/column)
//within this matrix. First the function sets all entries to zero, which can be exlcuded in respect to the current
//atom by atom matching. After this the function replaces all entries in the same row and column of the current
//position by zeros. Only the entry of the current position is set to one.
//Return value "count_arcsleft" counts the number of arcs, which are still in the matrix.
    private void remove_redundant_arcs(int row, int column, Vector<Integer> MARCS) {

        //System.err.print("Betrachte: " + c_globalA.get(row*2+0) + c_globalA.get(row*2+1));
        //System.err.println( " und " + c_globalA.get(column*2+0) + c_globalA.get(column*2+1) );
        int G1_atom = i_globalA.get(row * 3 + 0);
        int G2_atom = i_globalA.get(row * 3 + 1);
        int G3_atom = i_globalB.get(column * 3 + 0);
        int G4_atom = i_globalB.get(column * 3 + 1);

        for (int x = 0; x < nNum_globalA; x++) {
            int row_atom1 = i_globalA.get(x * 3 + 0);
            int row_atom2 = i_globalA.get(x * 3 + 1);

            for (int y = 0; y < nNum_globalB; y++) {
                int column_atom3 = i_globalB.get(y * 3 + 0);
                int column_atom4 = i_globalB.get(y * 3 + 1);

                if (((G1_atom == row_atom1) || (G1_atom == row_atom2)) && (!(((column_atom3 == G3_atom) || (column_atom4 == G3_atom)) || ((column_atom3 == G4_atom) || (column_atom4 == G4_atom))))) {

                    MARCS.setElementAt(0, x * nNum_globalB + y);
                }

                if (((G2_atom == row_atom1) || (G2_atom == row_atom2)) && (!(((column_atom3 == G3_atom) || (column_atom4 == G3_atom)) || ((column_atom3 == G4_atom) || (column_atom4 == G4_atom))))) {
                    MARCS.setElementAt(0, x * nNum_globalB + y);
                }

                if (((G3_atom == column_atom3) || (G3_atom == column_atom4)) && (!(((row_atom1 == G1_atom) || (row_atom2 == G1_atom)) || ((row_atom1 == G2_atom) || (row_atom2 == G2_atom))))) {
                    MARCS.setElementAt(0, x * nNum_globalB + y);
                }

                if (((G4_atom == column_atom3) || (G4_atom == column_atom4)) && (!(((row_atom1 == G1_atom) || (row_atom2 == G1_atom)) || ((row_atom1 == G2_atom) || (row_atom2 == G2_atom))))) {
                    MARCS.setElementAt(0, x * nNum_globalB + y);
                }

            }
        }

        for (int v = 0; v < nNum_globalA; v++) {
            MARCS.set(v * nNum_globalB + column, 0);
        }

        for (int w = 0; w < nNum_globalB; w++) {
            MARCS.set(row * nNum_globalB + w, 0);
        }

        MARCS.set(row * nNum_globalB + column, 1);
        //System.err.println("MARCS: " + MARCS);
    }

    /* Mdified function call by ASAD in Java have to check
     *
     */
    private int remove_tree_structure(BinaryTree cur_struc) {

        BinaryTree equal_struc = cur_struc.equal;
        BinaryTree not_equal_struc = cur_struc.not_equal;
        cur_struc = null;


        if (equal_struc != null) {
            remove_tree_structure(equal_struc);
        }

        if (not_equal_struc != null) {
            remove_tree_structure(not_equal_struc);
        }

        return 0;
    }

    //The function is called in function partsearch. The function is given z temporary matrix.
//The function checks whether the temporary matrix is already found by calling the function
//"verify_nodes". If the matrix already exists the function returns false which means that
//the matrix will not be stored. Otherwise the function returns true which means that the
//matrix will be stored in function partsearch.
    //private boolean check_MARCS(Vector<Integer> MARCS) {
    private boolean check_MARCS(Vector<Integer> MARCS_org) {


        Vector<Integer> MARCS_T = new Vector<Integer>(MARCS_org);

        //Vector<Integer> MARCS_T=MARCS;

        Vector<Integer> posnum_list = new Vector<Integer>();
        posnum_list.setSize(nNum_globalA * nNum_globalA); /*TO DO ASAD Initi by 0;*/

        for (int i = 0; i < posnum_list.size(); i++) {
            posnum_list.set(i, 0);
        }

        int y = 0;
        int count_entries = 0;
        for (int x = 0; x < (nNum_globalA * nNum_globalB); x++) {
            if (MARCS_T.elementAt(x) == (1)) {
                posnum_list.setElementAt(x, y++);
                count_entries++;

            }


        }



        boolean flag = false;

        verify_nodes(posnum_list, first, 0, count_entries);
        if (new_matrix) {
            flag = true;
        }

        return flag;

    }

    /* Modified function call by ASAD in Java have to check
     *
     */
    private boolean verify_nodes(Vector<Integer> matrix_org, BinaryTree cur_struc, int x, int field_length) {
        //private boolean verify_nodes(Vector<Integer> matrix, BinaryTree cur_struc, int x, int field_length) {

        Vector<Integer> matrix = new Vector<Integer>(matrix_org);

        if ((matrix.elementAt(x) == cur_struc.getValue()) && (x < field_length)) {
            if (cur_struc.equal != null) {
                new_matrix = false;
                verify_nodes(matrix, cur_struc.equal, x + 1, field_length);
            }

        }
        if (matrix.elementAt(x) != cur_struc.getValue()) {
            if (cur_struc.not_equal != null) {
                verify_nodes(matrix, cur_struc.not_equal, x, field_length);
            }

            if (cur_struc.not_equal == null) {
                cur_struc.not_equal = new BinaryTree();
                cur_struc.not_equal.setValue(matrix.elementAt(x));
                cur_struc.not_equal.not_equal = null;
                int y = 0;


                BinaryTree last_one = cur_struc.not_equal;

                while ((y + x + 1) < field_length) {
                    last_one.equal = new BinaryTree();
                    last_one =
                            last_one.equal;
                    last_one.setValue(matrix.elementAt(y + x + 1));
                    last_one.not_equal = null;
                    y++;

                }


                last_one.equal = null;
                new_matrix = true;
            }

        }
        return true;
    }

    private int search_corresponding_atom(int mapped_atoms_size, int atom_from_other_molecule, int molecule, Vector<Integer> mapped_atoms_org) {


        Vector<Integer> mapped_atoms = new Vector<Integer>(mapped_atoms_org);

        int corresponding_atom = 0;
        for (int a = 0; a < mapped_atoms_size; a++) {
            if (molecule == 1) {
                if (mapped_atoms.elementAt(a * 2 + 0).intValue() == atom_from_other_molecule) {

                    corresponding_atom = mapped_atoms.get(a * 2 + 1);
                }

            }
            if (molecule == 2) {
                if (mapped_atoms.elementAt(a * 2 + 1).intValue() == atom_from_other_molecule) {
                    corresponding_atom = mapped_atoms.get(a * 2 + 0);
                }

            }
        }
        return corresponding_atom;
    }

    // private int generate_c_setB_copy(int bond_number, Vector<String> c_setB_org, Vector<String> c_setA_copy) {
    private int generate_c_setB_copy(int bond_number, Vector<String> c_setB, Vector<String> c_setB_copy) {

        //Vector<String> c_setA = new Vector<String>(c_setB_org);

        for (int a = 0; a < bond_number; a++) {
            c_setB_copy.add(c_setB.get(a * 4 + 0));
            c_setB_copy.add(c_setB.get(a * 4 + 1));
            c_setB_copy.add("X");
            c_setB_copy.add("X");
        }

        return 0;
    }

    private int generate_c_setA_copy(int bond_number, Vector<String> c_setA, Vector<String> c_setA_copy) {

        //Vector<String> c_setA = new Vector<String>(c_setB_org);

        for (int a = 0; a < bond_number; a++) {
            c_setA_copy.add(c_setA.get(a * 4 + 0));
            c_setA_copy.add(c_setA.get(a * 4 + 1));
            c_setA_copy.add("X");
            c_setA_copy.add("X");
        }

        return 0;
    }

    private void startsearch() {
        Vector<Integer> FIXARCS = new Vector<Integer>(nNum_globalA * nNum_globalB);//  Initialize FIXARCS with 0

        FIXARCS.setSize(nNum_globalA * nNum_globalB);

        for (int i = 0; i < nNum_globalA * nNum_globalB; i++) {

            FIXARCS.set(i, 0);
        }

        int x = 0;
        int y = 0;

        while ((x < nNum_globalA) && (MARCS.elementAt(x * nNum_globalB + y) != 1)) {
            y++;
            if (y == nNum_globalB) {
                y = 0;
                x++;

            }

        }

        if (x == nNum_globalA) {
            //vorhergehende Schleife springt am Ende ber das Ziel hinaus
            y = nNum_globalB - 1;
            x = x - 1;
        }
        if (MARCS.elementAt(x * nNum_globalB + y) == 0) {
            partsearch(x, y, MARCS);
        }
        if (MARCS.elementAt(x * nNum_globalB + y) != 0) {
            partsearch(x, y, MARCS);
            MARCS.set(x * nNum_globalB + y, 0);
            partsearch(x, y, MARCS);
        }

    }

    public Vector<Vector<Integer>> getMappings() {

        return this._mappings;
    }

    public int getMCSSize() {

        return this.globalMCSSize;
    }
}
