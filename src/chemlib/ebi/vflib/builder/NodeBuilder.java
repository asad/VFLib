/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.vflib.builder;

import chemlib.ebi.vflib.interfaces.IAtomMatcher;
import chemlib.ebi.vflib.interfaces.IEdge;
import chemlib.ebi.vflib.interfaces.INode;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class NodeBuilder implements INode {

    private List<INode> neighbors;
    private List<IEdge> edges;
    private IAtomMatcher matcher;

    public NodeBuilder(IAtomMatcher matcher) {
        edges = new ArrayList<IEdge>();
        neighbors = new ArrayList<INode>();
        this.matcher = matcher;
    }

    @Override
    public int countNeighbors() {
        return neighbors.size();
    }

    @Override
    public Iterable<INode> neighbors() {
        return neighbors;
    }

    @Override
    public IAtomMatcher getAtomMatcher() {
        return matcher;
    }

    /**
     * 
     * @return
     */
    @Override
    public List<IEdge> getEdges() {
        return edges;
    }

    @Override
    public void addEdge(EdgeBuilder edge) {
        edges.add(edge);
    }

    @Override
    public void addneighbor(NodeBuilder node) {
        neighbors.add(node);
    }
}

