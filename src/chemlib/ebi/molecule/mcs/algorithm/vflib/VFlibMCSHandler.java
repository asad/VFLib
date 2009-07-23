/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.mcs.algorithm.vflib;

import chemlib.ebi.core.tools.EBIException;

import chemlib.ebi.interfaces.IMCS;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IMapper;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.INode;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IQuery;
import chemlib.ebi.molecule.mcs.algorithm.vflib.query.TemplateCompiler;
import chemlib.ebi.molecule.mcs.helper.MolHandler;
import chemlib.ebi.molecule.mcs.algorithm.vflib.map.VFMCSMapper;
import chemlib.ebi.molecule.mcs.algorithm.vflib.mcs.McGregorVFLib;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author sar
 */
public class VFlibMCSHandler implements IMCS {

    private static Vector<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static Vector<Map<IAtom, IAtom>> allAtomMCS_copy = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static Vector<TreeMap<Integer, Integer>> allMCS = null;
    private static Vector<TreeMap<Integer, Integer>> allMCS_copy = null;
    private boolean removeHydrogen = false;
    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;
    private boolean RONP = false;
    private IQuery query = null;
    private IMapper mapper = null;
    private List<Map<INode, IAtom>> vfLibSolutions = null;
    private int VFMCSSize = 0;

    public VFlibMCSHandler() {


        allAtomMCS = new Vector<Map<IAtom, IAtom>>();
        allAtomMCS_copy = new Vector<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new Vector<TreeMap<Integer, Integer>>();

        allMCS_copy = new Vector<TreeMap<Integer, Integer>>();


    }

    /**
     *
     * @return true if Query/Reactant is a subgraph of Target/Product
     * else false
     * @throws java.io.IOException
     * @throws chemlib.ebi.core.tools.EBIException
     */
    @Override
    public int search_MCS() throws IOException, EBIException {



//        System.out.println("VF Solution Count: " + vfLibSolutions.size());

        for (Map<INode, IAtom> solution : vfLibSolutions) {

            Map<IAtom, IAtom> atomatomMapping = new HashMap<IAtom, IAtom>();
            TreeMap<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();

            for (Map.Entry<INode, IAtom> mapping : solution.entrySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;
                if (RONP) {
                    qAtom = query.getAtom(mapping.getKey());
                    tAtom = mapping.getValue();

                } else {
                    tAtom = query.getAtom(mapping.getKey());
                    qAtom = mapping.getValue();
                }

                Integer qIndex = new Integer(ac1.getAtomNumber(qAtom));
                Integer tIndex = new Integer(ac2.getAtomNumber(tAtom));

//
//                System.out.println("i:" + qIndex + " j:" + tIndex);
//                System.out.println("i:" + qAtom.getSymbol() + " j:" + tAtom.getSymbol());

                if (qIndex != null && tIndex != null) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {

                    throw new EBIException("Atom index pointing to NULL");
                }

            }
            if (!atomatomMapping.isEmpty()) {
                allAtomMCS_copy.add(atomatomMapping);
                allMCS_copy.add(indexindexMapping);
                this.VFMCSSize = atomatomMapping.size();
            }


        }

//
//        System.out.println("count of the VFLib Mapping Stored: " + allAtomMCS_copy.size());
//        System.out.println("Mapping Size of the VFLib Mapping Stored: " + allAtomMCS_copy.firstElement().size());
//        System.out.println("GVC Mapping Size of the VFLib Mapping Stored: " + gvc.getMappingSize());

        boolean flag = mcgregorFlag();

        if (flag) {
            Vector<Vector<Integer>> _mappings = new Vector<Vector<Integer>>();

            for (TreeMap<Integer, Integer> firstPassMappings : allMCS_copy) {

//                System.out.println("\n\n---------------\n");
//                System.out.println("Calling McGregor");
//                System.out.println("for Start size " + firstPassMappings.size());
//                McGregorVFLib mgit = new McGregorVFLib(_mappings);
                McGregorVFLib mgit = new McGregorVFLib(ac1, ac2, _mappings);
                mgit.McGregor_IterationStart(firstPassMappings); //Start McGregorVFLib search

                _mappings = mgit.getMappings();
                mgit = null;
//                System.out.println("After Calling McGregorVFLib, solution Size" + _mappings.firstElement().size() / 2);
//

            }
//            System.out.println("Solution Count After McGregorVFLib " + _mappings.size());
//            System.out.println("Solution Size After Calling McGregorVFLib" + _mappings.firstElement().size() / 2);

            for (Vector<Integer> mapping : _mappings) {

                Map<IAtom, IAtom> atomatomMapping = new HashMap<IAtom, IAtom>();
                TreeMap<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();

                for (int index = 0; index < mapping.size(); index += 2) {
                    IAtom qAtom = null;
                    IAtom tAtom = null;

                    tAtom = ac1.getAtom(mapping.get(index));
                    qAtom = ac2.getAtom(mapping.get(index + 1));


                    Integer qIndex = mapping.get(index);
                    Integer tIndex = mapping.get(index + 1);


                    if (qIndex != null && tIndex != null) {
                        atomatomMapping.put(qAtom, tAtom);
                        indexindexMapping.put(qIndex, tIndex);
                    } else {

                        throw new EBIException("Atom index pointing to NULL");
                    }
                }

                if (!atomatomMapping.isEmpty()) {
                    if (!hasMap(indexindexMapping, allMCS)) {
                        allAtomMCS.add(atomatomMapping);
                        allMCS.add(indexindexMapping);
                    }
                }

            }


        } else if (!allAtomMCS_copy.isEmpty()) {


            allAtomMCS.addAll(allAtomMCS_copy);
            allMCS.addAll(allMCS_copy);
        }

        if (!allAtomMCS.isEmpty()) {
            atomsMCS.putAll(allAtomMCS.firstElement());
            firstMCS.putAll(allMCS.firstElement());


        }

//        System.out.println("Size of the Atom Mapping Stored: " + allAtomMCS.firstElement().size());
//        System.out.println("Size of the Atom Index Mapping Stored: " + allMCS.firstElement().size());

        return 0;

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

    private boolean mcgregorFlag() {
        int commonAtomCount = checkCommonAtomCount(ac1, ac2);
//        System.out.println("ra: " + ratomCount);
//        System.out.println("pa: " + patomCount);
//        System.out.println("ca: " + commonAtomCount);

        if (commonAtomCount > VFMCSSize && commonAtomCount > VFMCSSize) {
            return true;

        } else {
            return false;
        }
//       
    }

    /**
     * Set the JMCS software
     *
     * @param reactant
     * @param product
     * @param removeHydrogen
     */
    @Override
    public void set(IAtomContainer reactant, IAtomContainer product,
            boolean removeHydrogen) {

        this.removeHydrogen = removeHydrogen;
        IAtomContainer mol1 = reactant;
        IAtomContainer mol2 = product;

        MolHandler Reactant = new MolHandler(mol1, false, removeHydrogen);
        MolHandler Product = new MolHandler(mol2, false, removeHydrogen);

        //System.out.println(atom_number1+" " + atom_num_H_1);

        //best_clique_size = 0;

        this.set(Reactant, Product, removeHydrogen);
        // System.exit(1);

    }

    /**
     * Set the JMCS software
     *
     * @param removeHydrogen
     */
    @Override
    public void set(MolHandler Reactant, MolHandler Product, boolean removeHydrogen) {




        /*Remove Hydrogen by Asad*/
        this.removeHydrogen = removeHydrogen;

        ac1 = Reactant.getMolecule();
        ac2 = Product.getMolecule();



//        System.out.println("R Atom count " + ac1.getAtomCount());
//        System.out.println("P Atom count " + ac2.getAtomCount());

        if (ac1.getAtomCount() <= ac2.getAtomCount()) {

//            System.out.println("Query is Reactant");

            query = TemplateCompiler.compile(ac1);
            mapper = new VFMCSMapper(query);
            vfLibSolutions = mapper.getMaps(ac2);
            RONP = true;

        } else {
//            System.out.println("Query is Prouct");
            query = TemplateCompiler.compile(ac2);
            mapper = new VFMCSMapper(query);
            vfLibSolutions = mapper.getMaps(ac1);
            RONP = false;
        }




    }

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     * @param removeHydrogen
     */
    @Override
    public void set(String ReactantMolFileName, String ProductMolFileName,
            boolean removeHydrogen) {



        this.removeHydrogen = removeHydrogen;

        String mol1 = ReactantMolFileName;
        String mol2 = ProductMolFileName;

        MolHandler Reactant = new MolHandler(mol1, false, removeHydrogen);
        MolHandler Product = new MolHandler(mol2, false, removeHydrogen);

        //System.out.println(atom_number1+" " + atom_num_H_1);

        //best_clique_size = 0;

        this.set(Reactant, Product, removeHydrogen);
        // System.exit(1);



    }

    private boolean hasMap(Map<Integer, Integer> map, Vector<TreeMap<Integer, Integer>> mapGlobal) {

        for (Map<Integer, Integer> test : mapGlobal) {

            if (test.equals(map)) {
                return true;
            }
        }


        return false;
    }

    /**
     *
     * @return
     */
    @Override
    public Vector<Map<IAtom, IAtom>> getAllAtomMapping() {
        return allAtomMCS;
    }

    @Override
    public Vector<TreeMap<Integer, Integer>> getAllMapping() {
        return allMCS;
    }

    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        return atomsMCS;
    }

    @Override
    public TreeMap<Integer, Integer> getFirstMapping() {
        return firstMCS;
    }

    private int checkCommonAtomCount(IAtomContainer reactantMolecule, IAtomContainer productMolecule) {
        ArrayList<String> a = new ArrayList<String>();
        for (int i = 0; i < reactantMolecule.getAtomCount(); i++) {
            if (removeHydrogen && !reactantMolecule.getAtom(i).getSymbol().equals("H")) {
                a.add(reactantMolecule.getAtom(i).getSymbol());
            } else {
                a.add(reactantMolecule.getAtom(i).getSymbol());
            }
        }


        int common = 0;
        for (int i = 0; i < productMolecule.getAtomCount(); i++) {

            if (a.contains(productMolecule.getAtom(i).getSymbol())) {
                a.remove(productMolecule.getAtom(i).getSymbol());
                common++;
            }
        }
        return common;
    }
}
