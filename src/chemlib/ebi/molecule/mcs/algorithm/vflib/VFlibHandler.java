/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.mcs.algorithm.vflib;

import chemlib.ebi.core.tools.EBIAtomContainerManipulator;
import chemlib.ebi.core.tools.EBIException;

import chemlib.ebi.interfaces.ISubGraph;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IMapper;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.INode;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IQuery;
import chemlib.ebi.molecule.mcs.algorithm.vflib.query.TemplateCompiler;
import chemlib.ebi.molecule.mcs.helper.MolHandler;
import chemlib.ebi.molecule.mcs.algorithm.vflib.map.VFMapper;
import java.io.IOException;
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
public class VFlibHandler implements ISubGraph {

    private IAtomContainer Reactant;
    private IAtomContainer Product;
    private boolean RonPFlag = false;
    private static Vector<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static Vector<TreeMap<Integer, Integer>> allMCS = null;

    public VFlibHandler() {


        allAtomMCS = new Vector<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new Vector<TreeMap<Integer, Integer>>();

    }

    /**
     *
     * @return true if Query/Reactant is a subgraph of Target/Product
     * else false
     * @throws java.io.IOException
     * @throws chemlib.ebi.core.tools.EBIException
     */
    @Override
    public boolean isSubgraph() throws IOException, EBIException {

        IQuery query = TemplateCompiler.compile(Reactant);

        IMapper mapper = new VFMapper(query);

        List<Map<INode, IAtom>> vfLibSolutions = mapper.getMaps(Product);

//        System.out.println("Size of the Mapping: " + vfLibSolutions.size());

        for (Map<INode, IAtom> solution : vfLibSolutions) {

            Map<IAtom, IAtom> atomatomMapping = new HashMap<IAtom, IAtom>();
            TreeMap<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();

            for (Map.Entry<INode, IAtom> mapping : solution.entrySet()) {

                IAtom qAtom = query.getAtom(mapping.getKey());
                IAtom tAtom = mapping.getValue();

                Integer qIndex = Reactant.getAtomNumber(qAtom);
                Integer tIndex = Product.getAtomNumber(tAtom);

                atomatomMapping.put(qAtom, tAtom);
                indexindexMapping.put(qIndex, tIndex);


            }
            if (atomatomMapping.size() > 0) {
                allAtomMCS.add(atomatomMapping);
                allMCS.add(indexindexMapping);
            }


        }
        if (allAtomMCS.size() > 0) {
            atomsMCS.putAll(allAtomMCS.firstElement());
            firstMCS.putAll(allMCS.firstElement());
        } else {
            return false;
        }

        return true;

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

    /**
     * Set the JMCS software
     *
     * @param reactant
     * @param product
     * @param removeHydrogen
     */
    @Override
    public void set(IAtomContainer reactant, IAtomContainer product, boolean removeHydrogen) {

        this.Reactant = reactant;
        this.Product = product;

        /*Remove Hydrogen by Asad*/
        if (checkForH(Reactant) > 0 && removeHydrogen) {
            Reactant = EBIAtomContainerManipulator.removeHydrogens(reactant);
        }
        if (checkForH(Product) > 0 && removeHydrogen) {
            Product = EBIAtomContainerManipulator.removeHydrogens(product);
        }

    }

    /**
     * Set the JMCS software
     *
     * @param reactant
     * @param product
     * @param removeHydrogen
     */
    @Override
    public void set(MolHandler reactant, MolHandler product, boolean removeHydrogen) {


        this.Reactant = reactant.getMolecule();
        this.Product = product.getMolecule();

        /*Remove Hydrogen by Asad*/
        if (checkForH(Reactant) > 0 && removeHydrogen) {
            Reactant = EBIAtomContainerManipulator.removeHydrogens(reactant.getMolecule());
        }
        if (checkForH(Product) > 0 && removeHydrogen) {
            Product = EBIAtomContainerManipulator.removeHydrogens(product.getMolecule());
        }


    }

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     * @param removeHydrogen
     */
    @Override
    public void set(String ReactantMolFileName, String ProductMolFileName, boolean removeHydrogen) {


        String mol1 = ReactantMolFileName;
        String mol2 = ProductMolFileName;

        this.Reactant = new MolHandler(mol1, false, removeHydrogen).getMolecule();
        this.Product = new MolHandler(mol2, false, removeHydrogen).getMolecule();

        if (checkForH(Reactant) > 0 && removeHydrogen) {
            Reactant = EBIAtomContainerManipulator.removeHydrogens(Reactant);
        }
        if (checkForH(Product) > 0 && removeHydrogen) {
            Product = EBIAtomContainerManipulator.removeHydrogens(Product);
        }

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
}
