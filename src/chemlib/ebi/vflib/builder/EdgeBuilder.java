/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.vflib.builder;

import chemlib.ebi.vflib.interfaces.IBondMatcher;
import chemlib.ebi.vflib.interfaces.IEdge;
import chemlib.ebi.vflib.interfaces.INode;

/**
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class EdgeBuilder implements IEdge {

    private NodeBuilder source;
    private NodeBuilder target;
    private IBondMatcher matcher;

    /**
     * 
     * @param source
     * @param target
     * @param matcher
     */
    public EdgeBuilder(NodeBuilder source, NodeBuilder target, IBondMatcher matcher) {
        this.source = source;
        this.target = target;
        this.matcher = matcher;
    }

    @Override
    public INode getSource() {
        return source;
    }

    @Override
    public INode getTarget() {
        return target;
    }

    @Override
    public IBondMatcher getBondMatcher() {
        return matcher;
    }
}
