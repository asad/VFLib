/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.mcs.algorithm.vflib.map;

import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IMapper;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.INode;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IQuery;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IState;
import chemlib.ebi.molecule.mcs.algorithm.vflib.query.TemplateCompiler;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/*
 * SMSD (Small Molecule Subgraph Detector)
 *
 * Copyright (c) 2009 EMBL-EBI, Cambridge, UK
 *
 * Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
/**
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */

/*
 * @author Syed Asad Rahman <asad @ ebi.ac.uk> 
 */
public class VFMCSMapper implements IMapper {

    private IQuery query;
    private List<Map<INode, IAtom>> maps;
    private int StoredMCSSize = -1;

    /**
     *
     * @param query
     */
    public VFMCSMapper(IQuery query) {
        this.query = query;
        this.maps = new ArrayList<Map<INode, IAtom>>();
    }

    /**
     *
     * @param molecule
     */
    public VFMCSMapper(IAtomContainer molecule) {
        this.query = TemplateCompiler.compile(molecule);
        this.maps = new ArrayList<Map<INode, IAtom>>();
    }

    /**
     *
     * @param target
     * @return
     */
    @Override
    public boolean hasMap(IAtomContainer target) {
        VFMCSState state = new VFMCSState(query, target);
        maps.clear();

        return mapFirst(state);
    }

    @Override
    public List<Map<INode, IAtom>> getMaps(IAtomContainer target) {
        VFMCSState state = new VFMCSState(query, target);

        maps.clear();

        mapAll(state);

        return new ArrayList<Map<INode, IAtom>>(maps);
    }

    /**
     *
     * @param target
     * @return
     */
    @Override
    public Map<INode, IAtom> getFirstMap(IAtomContainer target) {
        VFMCSState state = new VFMCSState(query, target);

        maps.clear();

        mapFirst(state);

        return maps.isEmpty() ? new HashMap<INode, IAtom>() : maps.get(0);
    }

    /**
     *
     * @param target
     * @return
     */
    @Override
    public int countMaps(IAtomContainer target) {
        VFMCSState state = new VFMCSState(query, target);

        maps.clear();

        mapAll(state);

        return maps.size();
    }

    private void mapAll(IState state) {
        if (state.isDead()) {
            return;
        }

        if (state.isGoal()) {

            Map<INode, IAtom> map = state.getMap();
            if (isMCS(map) && !hasMap(map)) {
                maps.add(state.getMap());
            }
            return;
        }


        while (state.hasNextCandidate()) {
            VFMatch candidate = state.nextCandidate();
            if (state.isMatchFeasible(candidate)) {
                IState nextState = state.nextState(candidate);
                mapAll(nextState);
                nextState.backTrack();

                Map<INode, IAtom> map = state.getMap();
//                Added to Store MCS Solutions

                if (!map.isEmpty() && isMCS(map) && !hasMap(map)) {//!hasSubGraph(map)) {
                    maps.add(state.getMap());
                }


            }
        }
    }

    private boolean mapFirst(IState state) {
        if (state.isDead()) {
            return false;
        }

        if (state.isGoal()) {
            maps.add(state.getMap());
            return true;
        }

        boolean found = false;

        while (!found && state.hasNextCandidate()) {
            VFMatch candidate = state.nextCandidate();
            if (state.isMatchFeasible(candidate)) {
                IState nextState = state.nextState(candidate);
                found = mapFirst(nextState);

                nextState.backTrack();

                Map<INode, IAtom> map = state.getMap();
                if (!map.isEmpty() && isMCS(map) && !hasMap(map)) {//!hasSubGraph(map)
                    maps.add(state.getMap());
                }
            }
        }

        return found;
    }

    //Method added by Asad
    private boolean hasSubGraph(Map<INode, IAtom> map) {

        for (Map<INode, IAtom> storedMappings : maps) {
            boolean MatchFlag = true;
            for (Map.Entry<INode, IAtom> mapping : storedMappings.entrySet()) {
                if (!map.containsKey(mapping.getKey()) && !map.containsValue(mapping.getValue())) {
                    MatchFlag = false;
                    break;
                }
            }

            if (MatchFlag) {

                return true;
            }

        }

        return false;
    }

    //Method added by Asad
    private boolean isMCS(Map<INode, IAtom> map) {

        boolean flag = true;
        int mapSize = map.size();

        if (!maps.isEmpty() && StoredMCSSize > mapSize) {

            flag = false;
        }
        //Comment this if to get all the subgraphs
        if (mapSize > StoredMCSSize) {

            StoredMCSSize = mapSize;
            maps.clear();
        }

        return flag;
    }

    private boolean hasMap(Map<INode, IAtom> map) {
        for (Map<INode, IAtom> test : maps) {
            if (test.equals(map)) {
                return true;
            }
        }

        return false;
    }
}
