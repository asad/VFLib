/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package smsd.algorithm.vflib;

import java.util.Map;
import junit.framework.Assert;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import smsd.AtomAtomMapping;

/**
 *
 * @author Asad
 */
public class VF2SubTest {

    public VF2SubTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test ring match using MCS VF2Plus
     * @throws Exception
     */
    @Test
    public void testVF2Substructure() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Benzene
        IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
        // Napthalene
        IAtomContainer target = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
        //Algorithm is VF2MCS
        //Bond Sensitive is set True
        //Ring Match is set True

        VF2Sub comparison = new VF2Sub(true, true);
        comparison.set(query, target);
        Assert.assertTrue(comparison.isSubgraph());

        Assert.assertEquals(6, comparison.getFirstAtomMapping().getCount());
        Assert.assertEquals(24, comparison.getAllAtomMapping().size());

        /**
         * Print the mapping between molecules
         **/
        System.out.println(" Mappings: ");
        for (AtomAtomMapping atomatomMapping : comparison.getAllAtomMapping()) {
            for (Map.Entry<IAtom, IAtom> mapping : atomatomMapping.getMappings().entrySet()) {
                IAtom sourceAtom = mapping.getKey();
                IAtom targetAtom = mapping.getValue();
                System.out.println(sourceAtom.getSymbol() + " " + targetAtom.getSymbol());
                System.out.println(atomatomMapping.getQueryIndex(sourceAtom) + " " + atomatomMapping.getTargetIndex(targetAtom));
            }
            System.out.println("");
        }
    }
}
