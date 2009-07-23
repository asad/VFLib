/*
 * LabelContainer.java
 *
 * Created on 26 Nov 2007, 11:38
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package chemlib.ebi.molecule.mcs.helper;

import java.util.HashMap;
import java.util.Map;

public class LabelContainer {

    private HashMap<String, Integer> _labelMap;
    private int _LabelCounter = 1;
    private static LabelContainer INSTANCE = null;

    protected LabelContainer() {

        // System.err.println("Vector Initialized");
        _labelMap = new HashMap<String, Integer>();

        _labelMap.put("C", _LabelCounter++);
        _labelMap.put("O", _LabelCounter++);
        _labelMap.put("N", _LabelCounter++);
        _labelMap.put("S", _LabelCounter++);

        _labelMap.put("P", _LabelCounter++);
        _labelMap.put("F", _LabelCounter++);
        _labelMap.put("I", _LabelCounter++);
        _labelMap.put("R", _LabelCounter++);

        _labelMap.put("Br", _LabelCounter++);
        _labelMap.put("Cl", _LabelCounter++);
        _labelMap.put("Co", _LabelCounter++);
        _labelMap.put("Fe", _LabelCounter++);
        _labelMap.put("Na", _LabelCounter++);
        _labelMap.put("Ca", _LabelCounter++);

        _labelMap.put("K", _LabelCounter++);
        _labelMap.put("Mg", _LabelCounter++);
        _labelMap.put("Se", _LabelCounter++);
        _labelMap.put("Cu", _LabelCounter++);
        _labelMap.put("Hg", _LabelCounter++);

        _labelMap.put("X", _LabelCounter++);
        _labelMap.put("R", _LabelCounter++);
        _labelMap.put("X1", _LabelCounter++);
        

    }

    synchronized public static LabelContainer getInstance() {
        if (INSTANCE == null) {

            INSTANCE = new LabelContainer();

        }
        return INSTANCE;
    }

    synchronized public void addLabel(String label) {
        //System.err.println("Vector added");
        int lableNo = _labelMap.size() + 1;
        _labelMap.put(label, lableNo);
    }

    synchronized public Integer getLabelID(String label) {

        int id = -1;

        if (!_labelMap.containsKey(label)) {
            int lableNo = _labelMap.size() + 1;
            _labelMap.put(label, lableNo);

        }

        id = _labelMap.get(label);


        return id;
    }

    synchronized public String getLabel(Integer labelID) {

        String id = null;
        boolean flag = false;

        for (Map.Entry map : _labelMap.entrySet()) {

            id = (String) map.getKey();

            if (labelID == map.getValue()) {
                flag = true;
                break;
            }
        }

        if (!flag) {
            id = null;
        }
        return id;
    }

    synchronized public int getSize() {
        return _labelMap.size();
    }
}
