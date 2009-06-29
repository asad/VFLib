/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.mcs.algorithm.vflib;

import chemlib.ebi.molecule.mcs.global.BondType;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IMapper;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.INode;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IQuery;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IState;
import chemlib.ebi.molecule.mcs.algorithm.vflib.map.VFMCSMapper;
import chemlib.ebi.molecule.mcs.algorithm.vflib.map.VFMapper;
import chemlib.ebi.molecule.mcs.algorithm.vflib.map.VFMatch;
import chemlib.ebi.molecule.mcs.algorithm.vflib.map.VFState;
import chemlib.ebi.molecule.mcs.algorithm.vflib.query.TemplateCompiler;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import junit.framework.TestCase;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author sar
 */
public class TestVFlib extends TestCase {

    private IMolecule hexane;
    private IQuery hexaneQuery;
    private IMolecule benzene;
    private IQuery benzeneQuery;
    private IMolecule ATP;
    private IMolecule ADP;

    public TestVFlib() {
    }

    @Override
    protected void setUp() throws Exception {

        ATP = prepareMolecule("/Users/sar/TestData/ATP.mol");
//        ADP = prepareMolecule("/Users/sar/TestData/Mother/C03374_H.mol");
        ADP = prepareMolecule("/Users/sar/TestData/ADP.mol");
//        ATP = prepareMolecule("/Users/sar/TestData/Mother/C05787_H.mol");

        BondType.getInstance().setBondSensitiveFlag(true);
        System.out.println("AtomSize ADP:" + ADP.getAtomCount());
        System.out.println("AtomSize ATP:" + ATP.getAtomCount());

//        IQuery qADP = new VFQuery(ADP);
        IQuery qADP = TemplateCompiler.compile(ADP);
//        IMapper mapper = new VFMapper(qADP);
        IMapper mapper = new VFMCSMapper(qADP);

//        IState s = new VFState(qADP, ATP);
//        VFMatch match = new VFMatch(qADP.getNode(0), ATP.getAtom(0));
//        System.out.println("State count " + s.nextState(match).getMap().size());
//
//        IState s1 = new VFState(qADP, ATP);
//        VFMatch match1 = new VFMatch(qADP.getNode(1), ATP.getAtom(1));
//        System.out.println("State count " + s1.nextState(match1).getMap().size());
//
//        IState s2 = new VFState(qADP, ATP);
//        VFMatch match2 = new VFMatch(qADP.getNode(2), ATP.getAtom(2));
//        System.out.println("State count " + s2.nextState(match2).getMap().size());

        List<Map<INode, IAtom>> maps = mapper.getMaps(ATP);

//        Map<INode, IAtom> maps = mapper.getFirstMap(ATP);

        if (maps.size() > 0) {
            System.out.println("Mapping cou nt " + maps.size() + ", Mapping Size: " + maps.get(0).size());
        } else {
            System.out.println("No Match reported");
        }



//        assertEquals(2, maps.size());

//        hexane = createHexane();
//        System.out.println(" " + hexane.getAtomCount());
//        hexaneQuery = new VFQuery(hexane);
//        System.out.println(" " + hexaneQuery.countNodes());
//        benzene = createBenzene();
//        pyridine = Molecules.createPyridine();
//        toluene4 = create4Toluene();
//        pyridazine = MoleculeKit.readMolfile("[NO NAME]\r\n  CHEMWRIT          2D\r\nCreated with ChemWriter - http://metamolecular.com/chemwriter\r\n  6  6  0  0  0  0  0  0  0  0  0 V2000\r\n    1.8322   -5.0815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    2.6982   -5.5815    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\r\n    3.5643   -5.0815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    3.5643   -4.0815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    2.6982   -3.5815    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\r\n    1.8322   -4.0815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n  1  2  2  0  0  0  0\r\n  2  3  1  0  0  0  0\r\n  3  4  2  0  0  0  0\r\n  4  5  1  0  0  0  0\r\n  5  6  2  0  0  0  0\r\n  6  1  1  0  0  0  0\r\nM  END");
//        chloroisoquinoline4 = MoleculeKit.readMolfile("[NO NAME]\r\n  CHEMWRIT          2D\r\nCreated with ChemWriter - http://metamolecular.com/chemwriter\r\n 11 12  0  0  0  0  0  0  0  0  0 V2000\r\n    0.8800   -1.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    1.7460   -2.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    2.6121   -1.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    2.6121   -0.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    1.7460   -0.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    0.8800   -0.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    3.4781   -2.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    4.3442   -1.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    4.3442   -0.7400    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\r\n    3.4781   -0.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    3.4781   -3.2400    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\r\n  1  2  2  0  0  0  0\r\n  2  3  1  0  0  0  0\r\n  3  4  2  0  0  0  0\r\n  4  5  1  0  0  0  0\r\n  5  6  2  0  0  0  0\r\n  6  1  1  0  0  0  0\r\n  3  7  1  0  0  0  0\r\n  7  8  2  0  0  0  0\r\n  8  9  1  0  0  0  0\r\n  9 10  2  0  0  0  0\r\n 10  4  1  0  0  0  0\r\n  7 11  1  0  0  0  0\r\nM  END");
//        chlorobenzene = MoleculeKit.readMolfile("[NO NAME]\r\n  CHEMWRIT          2D\r\nCreated with ChemWriter - http://metamolecular.com/chemwriter\r\n  7  7  0  0  0  0  0  0  0  0  0 V2000\r\n   -3.3359    0.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   -2.4699    0.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   -1.6038    0.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   -1.6038    1.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   -2.4699    2.2400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   -3.3359    1.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   -0.7378    2.2400    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\r\n  1  2  2  0  0  0  0\r\n  2  3  1  0  0  0  0\r\n  3  4  2  0  0  0  0\r\n  4  5  1  0  0  0  0\r\n  5  6  2  0  0  0  0\r\n  6  1  1  0  0  0  0\r\n  4  7  1  0  0  0  0\r\nM  END");
//        benzeneQuery = new VFQuery(benzene);
//        chlorobenzeneQuery = new VFQuery(chlorobenzene);

//        pyridazineQuery = new VFQuery(pyridazine);
//        naphthalene = Molecules.createNaphthalene();
//        pyridineQuery = new VFQuery(pyridine);
//        toluene = Molecules.createToluene();
//        phenol = Molecules.createPhenol();
//        tolueneQuery = new VFQuery(toluene);
//        acetone = Molecules.createAcetone();
//        acetoneQuery = new VFQuery(acetone);
//        propane = Molecules.createPropane();
//        propaneQuery = new VFQuery(propane);
//        cyclopropane = Molecules.createCyclopropane();
    }

    private org.openscience.cdk.Molecule prepareMolecule(String MOLString) throws Exception {

        // creates CDK Molecule object
        Molecule mol = null;
        try {

            // reads MOL file

            File f = new File(MOLString);
            if (f.exists()) {
                MDLV2000Reader mdlr = new MDLV2000Reader(new FileReader(f));
                mol = (Molecule) mdlr.read(new Molecule());

                CDKHueckelAromaticityDetector.detectAromaticity(mol);
            } else {
                throw new Exception("File not found");
            }

        } catch (Exception e) {
            e.printStackTrace();
            throw new Exception("Rendering of structure(s) failed");
        }

        return mol;
    }

    public void testItShouldFindAllMatchCandidatesInTheRootState() {
        IState state = new VFState(benzeneQuery, benzene);
        int count = 0;

        while (state.hasNextCandidate()) {
            state.nextCandidate();

            count++;

        }



        assertEquals(benzene.getAtomCount() * benzene.getAtomCount(), count);
    }

    public void testItShoudFindAllMatchCandidatesInThePrimaryState() {
        IState state = new VFState(benzeneQuery, benzene);
        VFMatch match = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState newState = state.nextState(match);
        List<VFMatch> candidates = new ArrayList<VFMatch>();

        while (newState.hasNextCandidate()) {
            candidates.add(newState.nextCandidate());
        }

        assertEquals(4, candidates.size());
    }

    public void testItShouldFindAllMatchCandidatesInTheSecondaryState() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);
        List<VFMatch> candidates = new ArrayList<VFMatch>();

        while (state2.hasNextCandidate()) {
            candidates.add(state2.nextCandidate());
        }

        assertEquals(1, candidates.size());
    }

    public void testItShouldMapAllAtomsInTheSecondaryState() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);

        Map<INode, IAtom> map = state2.getMap();

        assertEquals(2, map.size());
        assertEquals(benzene.getAtom(0), map.get(benzeneQuery.getNode(0)));
        assertEquals(benzene.getAtom(1), map.get(benzeneQuery.getNode(1)));
    }

    public void testItShouldFindAllMatchCandidatesFromTheTeriaryState() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);
        VFMatch match2 = new VFMatch(benzeneQuery.getNode(2), benzene.getAtom(2));
        IState state3 = state2.nextState(match2);
        List<VFMatch> candidates = new ArrayList<VFMatch>();

        while (state3.hasNextCandidate()) {
            candidates.add(state3.nextCandidate());
        }

        assertEquals(1, candidates.size());
    }

    public void testItShouldMapAllAtomsInTheTertiaryState() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);
        VFMatch match2 = new VFMatch(benzeneQuery.getNode(2), benzene.getAtom(2));
        IState state3 = state2.nextState(match2);
        Map<INode, IAtom> map = state3.getMap();

        assertEquals(3, map.size());
        assertEquals(benzene.getAtom(0), map.get(benzeneQuery.getNode(0)));
        assertEquals(benzene.getAtom(1), map.get(benzeneQuery.getNode(1)));
        assertEquals(benzene.getAtom(2), map.get(benzeneQuery.getNode(2)));
    }

    public void testItShouldReachGoalWhenAllAtomsAreMapped() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);
        VFMatch match2 = new VFMatch(benzeneQuery.getNode(2), benzene.getAtom(2));
        IState state3 = state2.nextState(match2);
        VFMatch match3 = new VFMatch(benzeneQuery.getNode(3), benzene.getAtom(3));
        IState state4 = state3.nextState(match3);
        VFMatch match4 = new VFMatch(benzeneQuery.getNode(4), benzene.getAtom(4));
        IState state5 = state4.nextState(match4);

        assertFalse(state5.isGoal());

        VFMatch match5 = new VFMatch(benzeneQuery.getNode(5), benzene.getAtom(5));
        IState state6 = state5.nextState(match5);

        assertTrue(state6.isGoal());
    }

    public void testItShouldHaveANextCandidateInTheSecondaryState() {
        IState state = new VFState(benzeneQuery, benzene);
        VFMatch match = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));

        IState nextState = state.nextState(match);

        assertTrue(nextState.hasNextCandidate());
    }

    public void testItShouldMatchHexaneToHexane() {
        IMapper mapper = new VFMapper(hexaneQuery);

        assertTrue(mapper.hasMap(hexane));
    }

    public void testItShouldMatchHexaneToHexaneWhenUsingMolecule() {
        IMapper mapper = new VFMapper(hexane);

        assertTrue(mapper.hasMap(hexane));
    }

    public void testItShouldFindTwoMapsFromHexaneToHexane() {
        IMapper mapper = new VFMapper(hexaneQuery);

        List<Map<INode, IAtom>> maps = mapper.getMaps(hexane);



        assertEquals(2, maps.size());
    }
//
//    public void testItShouldMatchBenzeneToBenzene() {
//        IMapper mapper = new VFMapper(benzeneQuery);
//
//        assertTrue(mapper.hasMap(benzene));
//    }
//
//    public void testItShouldNotMatchHexaneToBenzene() {
//        IMapper mapper = new VFMapper(hexaneQuery);
//
//        assertFalse(mapper.hasMap(benzene));
//    }
//
//    public void testItShouldNotMatchPyridazineToNaphthalene() {
//        IMapper mapper = new VFMapper(pyridazineQuery);
//
//        assertFalse(mapper.hasMap(naphthalene));
//    }
//
//    public void testItShouldNotMatchChlorobenzeneTo4ChloroIsoquinoline() {
//        IMapper mapper = new VFMapper(chlorobenzeneQuery);
//
//        assertFalse(mapper.hasMap(chloroisoquinoline4));
//    }
//
//    public void testItShouldNotMatchBenzeneToPyridine() {
//        IMapper mapper = new VFMapper(benzeneQuery);
//
//        assertFalse(mapper.hasMap(pyridine));
//
//        mapper = new VFMapper(pyridineQuery);
//
//        assertFalse(mapper.hasMap(benzene));
//    }
//
//    public void testItShouldNotMatchTolueneToBenzene() {
//        IMapper mapper = new VFMapper(tolueneQuery);
//
//        assertFalse(mapper.hasMap(benzene));
//    }
//
//    public void testItShouldMatchAcetoneToAcetone() {
//        IMapper mapper = new VFMapper(acetoneQuery);
//
//        assertTrue(mapper.hasMap(acetone));
//    }
//
//    public void testItShouldMatchPropaneToCyclopropane() {
//        IMapper mapper = new VFMapper(propaneQuery);
//
//        assertTrue(mapper.hasMap(cyclopropane));
//    }
//
//
//
//    public void testItShouldNotMatchTolueneToPhenol() {
//        IMapper mapper = new VFMapper(tolueneQuery);
//
//        assertFalse(mapper.hasMap(phenol));
//    }
//
//    public void testItShouldMapSixAtomsOfBenzeneOntoBenzene() {
//        IMapper mapper = new VFMapper(benzeneQuery);
//        Map<INode, Atom> map = mapper.getFirstMap(benzene);
//
//        assertEquals(6, map.size());
//    }
//
//    public void testItShouldCountTwelveMapsForBenzeneOntoBenzene() {
//        IMapper mapper = new VFMapper(benzeneQuery);
//
//        assertEquals(12, mapper.countMaps(benzene));
//    }
//
//    public void testItShouldCountTwoMapsForTolueneOntoToluene() {
//        IMapper mapper = new VFMapper(tolueneQuery);
//
//        assertEquals(2, mapper.countMaps(toluene));
//    }
//
//    public void testItShouldFindTwelveMapsForBenzeneOntoBenzene() {
//        IMapper mapper = new VFMapper(benzeneQuery);
//        List<Map<INode, Atom>> maps = mapper.getMaps(benzene);
//
//        assertEquals(12, maps.size());
//    }
//
//    public void testItShouldFindTwentyFourMapsForBenzeneOntoNaphthalene() {
//        IMapper mapper = new VFMapper(benzeneQuery);
//        List<Map<INode, Atom>> maps = mapper.getMaps(naphthalene);
//
//        assertEquals(24, maps.size());
//    }
//
//    public void testItShouldFindAMapForEquivalentFormsOfToluene() {
//        IMapper mapper = new VFMapper(tolueneQuery);
//        Map<INode, Atom> map = mapper.getFirstMap(toluene4);
//
//        assertEquals(7, map.size());
//    }
//
//    public void testItShouldFindTwoMapsForEquivalentFormsOfToluene() {
//        IMapper mapper = new VFMapper(tolueneQuery);
//        List<Map<INode, Atom>> maps = mapper.getMaps(toluene4);
//
//        assertEquals(2, maps.size());
//    }

    private IMolecule create4Toluene() throws CDKException {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();
        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c1.setID("1");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c2.setID("2");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c3.setID("3");
        IAtom c4 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("4");
        IAtom c5 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("5");
        IAtom c6 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("6");
        IAtom c7 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c7.setID("7");

        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);
        result.addAtom(c4);
        result.addAtom(c5);
        result.addAtom(c6);
        result.addAtom(c7);



        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.DOUBLE);
        IBond bond3 = new Bond(c3, c4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(c4, c5, IBond.Order.DOUBLE);
        IBond bond5 = new Bond(c5, c6, IBond.Order.SINGLE);
        IBond bond6 = new Bond(c6, c1, IBond.Order.DOUBLE);
        IBond bond7 = new Bond(c7, c4, IBond.Order.SINGLE);

        result.addBond(bond1);
        result.addBond(bond2);
        result.addBond(bond3);
        result.addBond(bond4);
        result.addBond(bond5);
        result.addBond(bond6);
        result.addBond(bond7);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(result);

        return result;
    }

    public IMolecule createMethane() {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();
        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        result.addAtom(c1);

        return result;
    }

    public IMolecule createPropane() {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();
        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");


        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);



        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.SINGLE);

        result.addBond(bond1);
        result.addBond(bond2);



        return result;
    }

//    public static Molecule createCyclopropane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c2, c0, 1);
//
//        return result;
//    }
    public IMolecule createHexane() throws CDKException {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();

        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c1.setID("1");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c2.setID("2");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c3.setID("3");
        IAtom c4 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("4");
        IAtom c5 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("5");
        IAtom c6 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("6");

        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);
        result.addAtom(c4);
        result.addAtom(c5);
        result.addAtom(c6);

        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.SINGLE);
        IBond bond3 = new Bond(c3, c4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(c4, c5, IBond.Order.SINGLE);
        IBond bond5 = new Bond(c5, c6, IBond.Order.SINGLE);


        result.addBond(bond1);
        result.addBond(bond2);
        result.addBond(bond3);
        result.addBond(bond4);
        result.addBond(bond5);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(result);

        return result;
    }
//
//    public static Molecule createCyclohexane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 1);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 1);
//        result.connect(c5, c0, 1);
//
//        return result;
//    }
//

    public static IMolecule createBenzene() throws CDKException {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();

        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c1.setID("1");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c2.setID("2");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c3.setID("3");
        IAtom c4 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("4");
        IAtom c5 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("5");
        IAtom c6 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("6");

        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);
        result.addAtom(c4);
        result.addAtom(c5);
        result.addAtom(c6);

        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.DOUBLE);
        IBond bond3 = new Bond(c3, c4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(c4, c5, IBond.Order.DOUBLE);
        IBond bond5 = new Bond(c5, c6, IBond.Order.SINGLE);
        IBond bond6 = new Bond(c6, c1, IBond.Order.DOUBLE);


        result.addBond(bond1);
        result.addBond(bond2);
        result.addBond(bond3);
        result.addBond(bond4);
        result.addBond(bond5);
        result.addBond(bond6);


        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(result);



        return result;
    }
//    public static Molecule createPyridine() {
//        Molecule result = new DefaultMolecule();
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("N");
//
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 2);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 2);
//        result.connect(c5, c6, 1);
//        result.connect(c6, c1, 2);
//
//        return result;
//    }
//
//    public static Molecule createToluene() {
//        Molecule result = new DefaultMolecule();
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("C");
//        Atom c7 = result.addAtom("C");
//
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 2);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 2);
//        result.connect(c5, c6, 1);
//        result.connect(c6, c1, 2);
//        result.connect(c7, c1, 1);
//
//        return result;
//    }
//
//    public static Molecule createPhenol() {
//        Molecule result = new DefaultMolecule();
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("C");
//        Atom c7 = result.addAtom("O");
//
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 2);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 2);
//        result.connect(c5, c6, 1);
//        result.connect(c6, c1, 2);
//        result.connect(c7, c1, 1);
//
//        return result;
//    }
//
//    public static Molecule createNaphthalene() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("C");
//        Atom c7 = result.addAtom("C");
//        Atom c8 = result.addAtom("C");
//        Atom c9 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 2);
//        result.connect(c2, c3, 1);
//        result.connect(c3, c4, 2);
//        result.connect(c4, c5, 1);
//        result.connect(c5, c0, 2);
//        result.connect(c4, c6, 1);
//        result.connect(c6, c7, 2);
//        result.connect(c7, c8, 1);
//        result.connect(c8, c9, 2);
//        result.connect(c9, c5, 1);
//
//        return result;
//    }
//
//    public static Molecule createAcetone() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom o3 = result.addAtom("O");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c1, o3, 2);
//
//        return result;
//    }
//
//    public static Molecule createNeopentane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c0, c2, 1);
//        result.connect(c0, c3, 1);
//        result.connect(c0, c4, 1);
//
//        return result;
//    }
//
//    public static Molecule createCubane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("C");
//        Atom c7 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 1);
//        result.connect(c3, c0, 1);
//
//        result.connect(c4, c5, 1);
//        result.connect(c5, c6, 1);
//        result.connect(c6, c7, 1);
//        result.connect(c7, c4, 1);
//
//        result.connect(c0, c4, 1);
//        result.connect(c1, c5, 1);
//        result.connect(c2, c6, 1);
//        result.connect(c3, c7, 1);
//
//        return result;
//    }
//
//    public static Molecule createBicyclo220hexane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 1);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 1);
//        result.connect(c5, c0, 1);
//        result.connect(c2, c5, 1);
//
//        return result;
//    }
//
//    public static Molecule createEthylbenzeneWithSuperatom() {
//        Molecule result = Molecules.createBenzene();
//        Atom carbon1 = result.addAtom("C");
//        Atom carbon2 = result.addAtom("C");
//        Bond crossingBond = result.connect(result.getAtom(0), carbon1, 1);
//        result.connect(carbon1, carbon2, 1);
//
//        Superatom substructure = result.addSuperatom();
//        substructure.addAtom(carbon1);
//        substructure.addAtom(carbon2);
//        substructure.addCrossingBond(crossingBond);
//        substructure.setCrossingVector(crossingBond, 0.1, 0.1);
//        substructure.setLabel("Ethyl");
//
//        return result;
//    }
}
