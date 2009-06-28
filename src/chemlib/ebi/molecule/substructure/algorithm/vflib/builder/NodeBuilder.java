/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.substructure.algorithm.vflib.builder;

import chemlib.ebi.molecule.substructure.algorithm.vflib.interfaces.IAtomMatcher;
import chemlib.ebi.molecule.substructure.algorithm.vflib.interfaces.IEdge;
import chemlib.ebi.molecule.substructure.algorithm.vflib.interfaces.INode;
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

    
    public int countNeighbors() {
        return neighbors.size();
    }

    
    public Iterable<INode> neighbors() {
        return neighbors;
    }

    
    public IAtomMatcher getAtomMatcher() {
        return matcher;
    }

    /**
     * 
     * @return
     */
    
    public List<IEdge> getEdges() {
        return edges;
    }

    
    public void addEdge(EdgeBuilder edge) {
        edges.add(edge);
    }

    
    public void addneighbor(NodeBuilder node) {
        neighbors.add(node);
    }
}

