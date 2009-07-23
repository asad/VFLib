/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.mcs.algorithm.vflib.map;

import chemlib.ebi.molecule.mcs.algorithm.vflib.validator.VFMatch;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IEdge;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.INode;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IQuery;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IState;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class VFMCSState implements IState {

    private List<VFMatch> candidates;
    private IQuery query;
    private IAtomContainer target;
    private List<INode> queryPath;
    private List<IAtom> targetPath;
    private Map<INode, IAtom> map;

    public VFMCSState(IQuery query, IAtomContainer target) {
        this.map = new HashMap<INode, IAtom>();
        this.queryPath = new ArrayList<INode>();
        this.targetPath = new ArrayList<IAtom>();

        this.query = query;
        this.target = target;
        this.candidates = new ArrayList<VFMatch>();

        loadRootCandidates();
    }

    private VFMCSState(VFMCSState state, VFMatch match) {
        this.candidates = new ArrayList<VFMatch>();
        this.queryPath = new ArrayList<INode>(state.queryPath);
        this.targetPath = new ArrayList<IAtom>(state.targetPath);
        this.map = state.map;
        this.query = state.query;
        this.target = state.target;

        map.put(match.getQueryNode(), match.getTargetAtom());
        queryPath.add(match.getQueryNode());
        targetPath.add(match.getTargetAtom());

        loadCandidates(match);
    }

    @Override
    public void backTrack() {
        if (queryPath.isEmpty() || isGoal()) {
            map.clear();

            return;
        }
        /*Commented by Asad as its seems to map more elements in MCS*/
        if (isHeadMapped()) {
            return;
        }

        map.clear();

        for (int i = 0; i < queryPath.size() - 1; i++) {
            map.put(queryPath.get(i), targetPath.get(i));
        }

    }

    @Override
    public Map<INode, IAtom> getMap() {
        return new HashMap<INode, IAtom>(map);
    }

    @Override
    public boolean hasNextCandidate() {
        return !candidates.isEmpty();
    }

    @Override
    public boolean isDead() {
        return query.countNodes() > target.getAtomCount();
    }

    @Override
    public boolean isGoal() {
        return map.size() == query.countNodes();
    }

    @Override
    public boolean isMatchFeasible(VFMatch match) {
        if (map.containsKey(match.getQueryNode()) || map.containsValue(match.getTargetAtom())) {
            return false;
        }

        if (!matchAtoms(match)) {
            return false;
        }

        if (!matchBonds(match)) {
            return false;
        }

        return true;
    }

    @Override
    public VFMatch nextCandidate() {
        return candidates.remove(candidates.size() - 1);
    }

    @Override
    public IState nextState(VFMatch match) {
        return new VFMCSState(this, match);
    }

    private void loadRootCandidates() {
//        System.out.println("Node count " + query.countNodes());
//        System.out.println("Target count " + target.getAtomCount());
        for (int i = 0; i < query.countNodes(); i++) {
            for (int j = 0; j < target.getAtomCount(); j++) {
                VFMatch match = new VFMatch(query.getNode(i), target.getAtom(j));
                candidates.add(match);
            }
        }

//        System.out.println("candidates count " + candidates.size());
    }

//@TODO Asad Check the Neighbour count
    private void loadCandidates(VFMatch lastMatch) {
        IAtom a = lastMatch.getTargetAtom();
        List<IAtom> targetNeighbors = target.getConnectedAtomsList(a);

        for (INode q : lastMatch.getQueryNode().neighbors()) {
            for (IAtom t : targetNeighbors) {
                VFMatch match = new VFMatch(q, t);

                if (candidateFeasible(match)) {
                    candidates.add(match);
                }
            }
        }


    }

    private boolean candidateFeasible(VFMatch candidate) {
        for (INode queryAtom : map.keySet()) {
            if (queryAtom.equals(candidate.getQueryNode()) ||
                    map.get(queryAtom).equals(candidate.getTargetAtom())) {
                return false;
            }
        }

        return true;
    }

    private boolean matchAtoms(VFMatch match) {
        IAtom a = match.getTargetAtom();
        List<IAtom> targetNeighbors = target.getConnectedAtomsList(a);
        if (match.getQueryNode().countNeighbors() > targetNeighbors.size()) {
            return false;
        }


        return match.getQueryNode().getAtomMatcher().matches(match.getTargetAtom());
    }

    private boolean matchBonds(VFMatch match) {
        if (queryPath.isEmpty()) {
            return true;
        }

        if (!matchBondsToHead(match)) {
            return false;
        }

        for (int i = 0; i < queryPath.size() - 1; i++) {
            IEdge queryBond = query.getEdge(queryPath.get(i), match.getQueryNode());
            IBond targetBond = target.getBond(targetPath.get(i), match.getTargetAtom());

            if (queryBond == null) {
                continue;
            }

            if (targetBond == null) {
                continue;
            }

            if (!matchBond(queryBond, targetBond)) {
                return false;
            }
        }

        return true;
    }

    private boolean matchBond(IEdge edge, IBond targetBond) {


        return edge.getBondMatcher().matches(targetBond);


    }

    private boolean isHeadMapped() {
        INode head = queryPath.get(queryPath.size() - 1);
        for (INode neighbor : head.neighbors()) {
            if (!map.containsKey(neighbor)) {
                return false;
            }
        }

        return true;
    }

    private boolean matchBondsToHead(VFMatch match) {
        INode queryHead = getQueryPathHead();
        IAtom targetHead = getTargetPathHead();

        IEdge queryBond = query.getEdge(queryHead, match.getQueryNode());
        IBond targetBond = target.getBond(targetHead, match.getTargetAtom());

        if (queryBond == null || targetBond == null) {
            return false;
        }

        return matchBond(queryBond, targetBond);
    }

    private INode getQueryPathHead() {
        return queryPath.get(queryPath.size() - 1);
    }

    private IAtom getTargetPathHead() {
        return targetPath.get(targetPath.size() - 1);
    }
}
