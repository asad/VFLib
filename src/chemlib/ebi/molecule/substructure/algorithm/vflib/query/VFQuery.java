/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.substructure.algorithm.vflib.query;

import chemlib.ebi.molecule.substructure.algorithm.vflib.interfaces.IAtomMatcher;
import chemlib.ebi.molecule.substructure.algorithm.vflib.interfaces.IBondMatcher;
import chemlib.ebi.molecule.substructure.algorithm.vflib.interfaces.IEdge;
import chemlib.ebi.molecule.substructure.algorithm.vflib.interfaces.INode;
import chemlib.ebi.molecule.substructure.algorithm.vflib.interfaces.IQuery;
import chemlib.ebi.molecule.substructure.algorithm.vflib.builder.EdgeBuilder;
import chemlib.ebi.molecule.substructure.algorithm.vflib.builder.NodeBuilder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;

/*
* MX Cheminformatics Tools for Java
*
* Copyright (c) 2007-2009 Metamolecular, LLC
*
* http://metamolecular.com
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
 * @author Richard L. Apodaca <rapodaca at metamolecular.com>
 * @author Syed Asad Rahman <asad @ ebi.ac.uk> (modified the orignal code)
 */
public class VFQuery implements IQuery {

    private List<INode> nodes;
    private List<IEdge> edges;
    private Map<INode, IAtom> NodesBonds;

    public VFQuery() {
        nodes = new ArrayList<INode>();
        edges = new ArrayList<IEdge>();
        NodesBonds = new HashMap<INode, IAtom>();
    }

    /**
     *
     * @return
     */
    
    public Iterable<IEdge> edges() {
        return edges;
    }

    /**
     *
     * @return
     */
    
    public Iterable<INode> nodes() {
        return nodes;
    }

    /**
     *
     * @param index
     * @return
     */
    
    public INode getNode(int index) {
        return nodes.get(index);
    }

    /**
     *
     * @param index
     * @return
     */
    
    public IEdge getEdge(int index) {
        return edges.get(index);
    }

    /**
     *
     * @param source
     * @param target
     * @return
     */
    
    public IEdge getEdge(INode source, INode target) {
        if (source == target) {
            return null;
        }

        NodeBuilder sourceImpl = (NodeBuilder) source;

        for (IEdge edge : sourceImpl.getEdges()) {
            if (edge.getSource() == target || edge.getTarget() == target) {
                return edge;
            }
        }

        return null;
    }

    /**
     *
     * @param matcher
     * @return
     */
    public INode addNode(IAtomMatcher matcher) {
        NodeBuilder node = new NodeBuilder(matcher);

        nodes.add(node);
        return node;
    }

    /**
     *
     * @param matcher
     * @param atom
     * @return
     */
    public INode addNode(IAtomMatcher matcher, IAtom atom) {
        NodeBuilder node = new NodeBuilder(matcher);
        NodesBonds.put(node, atom);
        nodes.add(node);
        return node;
    }

    /**
     *
     * @param node
     * @return
     */
    
    public IAtom getAtom(INode node) {

        return NodesBonds.get(node);
    }

    
    public int countNodes() {
        return nodes.size();
    }

    /**
     *
     * @return
     */
    
    public int countEdges() {
        return edges.size();
    }

    /**
     * 
     * @param source
     * @param target
     * @param matcher
     * @return
     */
    public IEdge connect(INode source, INode target, IBondMatcher matcher) {
        NodeBuilder sourceImpl = (NodeBuilder) source;
        NodeBuilder targetImpl = (NodeBuilder) target;
        EdgeBuilder edge = new EdgeBuilder(sourceImpl, targetImpl, matcher);

        sourceImpl.addneighbor(targetImpl);
        targetImpl.addneighbor(sourceImpl);

        sourceImpl.addEdge(edge);
        targetImpl.addEdge(edge);

        edges.add(edge);
        return edge;
    }
}
 