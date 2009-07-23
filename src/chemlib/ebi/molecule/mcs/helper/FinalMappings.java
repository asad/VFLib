/*
 * FinalMapping.java
 *
 * Created on 09 February 2007, 11:38
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package chemlib.ebi.molecule.mcs.helper;

import chemlib.ebi.interfaces.IFinalMapping;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Vector;

public class FinalMappings implements IFinalMapping {

    private Vector<TreeMap<Integer, Integer>> _mappings;
    private static FinalMappings INSTANCE = null;

    protected FinalMappings() {
        _mappings = new Vector<TreeMap<Integer, Integer>>();
    }

    synchronized public static FinalMappings getInstance() {
        if (INSTANCE == null) {

            INSTANCE = new FinalMappings();

        }
        return INSTANCE;
    }

    @Override
    synchronized public void addElement(TreeMap<Integer, Integer> v) {
        //System.err.println("Vector added");
        _mappings.addElement(v);
    }

    @Override
    synchronized public void add(TreeMap<Integer, Integer> v) {
        //System.err.println("Vector added");
        _mappings.add(v);
    }

    /**
     * 
     * @param v
     */
    @Override
    synchronized public final void set(Vector<TreeMap<Integer, Integer>> v) {

//        System.err.println("FinalMapping Vector Size before set: " + _mappings.size());
//        System.err.println("Vector Size to be set: " + v.size());

        /*for (Vector<Integer> I : v) {
        _mappings.add(I);
        }*/

        _mappings.clear();

//        System.err.println("Vector Initialized");

        for (TreeMap<Integer, Integer> M : v) {
            _mappings.add(M);
        }
//        System.err.println("FinalMapping Vector Size after set: " + _mappings.size());
    }

    @Override
    synchronized public Iterator<TreeMap<Integer, Integer>> getIterator() {
        Iterator<TreeMap<Integer, Integer>> it = _mappings.iterator();
        return it;
    }

    @Override
    synchronized public void Clear() {

        //System.err.println("Size before Clear: " + final_MAPPINGS.size());

        _mappings.clear();

    //System.err.println("Size after Clear: " + final_MAPPINGS.size());
    }

    @Override
    synchronized public Vector<TreeMap<Integer, Integer>> getFinalMapping() {
        return _mappings;
    }

    @Override
    synchronized public int getSize() {
        return _mappings.size();
    }
}
