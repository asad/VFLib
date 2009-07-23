/*
 * BinaryTree.java
 *
 * Created on January 27, 2007, 9:36 PM
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 *
 * @Revised by Asad, EBI, Java 1.5 Compatable
 */
package chemlib.ebi.molecule.mcs.helper;

//the second part of the program extents the mapping by the McGregor algorithm in case
//that not all atoms of molecule A and molecule B are mapped by the clique approach


public class BinaryTree {

    /** Creates a new instance of BinaryTree */
    public BinaryTree() {
        //  System.err.println("Binary Tree Created");
        equal = null;
        not_equal = null;
    }
    /**
     * equal is initialized as null
     */
    public BinaryTree equal = null;
    /**
     * not equal is initialized as null
     */
    public BinaryTree not_equal = null;
    
    private int Value = -1;

    /**
     * 
     * @param val set the value of the current node
     */
    public void setValue(int val) {
        //    System.err.println("Value Set: " + val);
        this.Value = val;
    }

    /**
     * 
     * @return get the value of the current node
     */
    public int getValue() {
        return this.Value;
    }
}
/**
 * A class representing a tree data structure with exactly 2 child elements or
 * none.
 *
 * @author Silvere Martin-Michiellot
 * @version 1.0
 *///perhaps we should add synchronization support here
//null contents for each node accepted
//we could have some methods for clustering here, see
//http://128.32.125.151/riot/Applications/Clustering/index.html or http://www2.biology.ualberta.ca/jbrzusto/cluster.php
//public class BinaryTree extends Object implements Tree, Cloneable, Serializable {
//    /** DOCUMENT ME! */
//    private Object obj;
//
//    /** DOCUMENT ME! */
//    private BinaryTree parent;
//
//    /** DOCUMENT ME! */
//    private BinaryTree child1;
//
//    /** DOCUMENT ME! */
//    private BinaryTree child2;
//
//    static final long serialVersionUID = -4544769666886838818L;
//
//
//
//    /**
//     * Creates a new BinaryTree object.
//     */
//    public BinaryTree() {
//        obj = null;
//        child1 = null;
//        child2 = null;
//    }
//
//    /**
//     * Creates a new BinaryTree object.
//     *
//     * @param o DOCUMENT ME!
//     */
//    public BinaryTree(Object o) {
//        obj = o;
//        child1 = null;
//        child2 = null;
//    }
//
//    //may return null
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public Object getContents() {
//        return obj;
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param obj DOCUMENT ME!
//     */
//    public void setContents(Object obj) {
//        this.obj = obj;
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public boolean hasChild() {
//        return (child1 != null) || (child2 != null);
//    }
//
//    //if tree is one of the children of this
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param child DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public boolean hasChild(BinaryTree child) {
//        return child1.equals(child) || child2.equals(child);
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public Set<BinaryTree>  getChildren() {
//        HashSet<BinaryTree> result;
//
//        result = new HashSet<BinaryTree>();
//        result.add(child1);
//        result.add(child2);
//
//        return result;
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public BinaryTree getLeftChild() {
//        return child1;
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public BinaryTree getRightChild() {
//        return child2;
//    }
//
//    //all elements should be of class BinaryTree
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param children DOCUMENT ME!
//     *
//     * @throws CircularReferenceException DOCUMENT ME!
//     * @throws IllegalArgumentException DOCUMENT ME!
//     */
//    public void setChildren(Set children) throws CircularReferenceException {
//        Object currentElement;
//        BinaryTree child;
//        Object[] childrenArray;
//        int i;
//
//        if (children != null) {
//            childrenArray = children.toArray();
//
//            if (childrenArray.length < 3) {
//                i = 0;
//
//                while (i < childrenArray.length) {
//                    currentElement = childrenArray[i];
//
//                    if (currentElement instanceof BinaryTree) {
//                        child = (BinaryTree) currentElement;
//
//                        if ((child.getParent() == null) && (child != this)) {
//                            child.setParent(this);
//                        } else {
//                            throw new CircularReferenceException(
//                                "Cannot add BinaryTree child that has already a parent.");
//                        }
//                    } else {
//                        throw new IllegalArgumentException(
//                            "The children Set must contain only BinaryTrees.");
//                    }
//                }
//
//                if (childrenArray.length == 2) {
//                    child1 = (BinaryTree) childrenArray[0];
//                    child2 = (BinaryTree) childrenArray[1];
//                } else {
//                    if (childrenArray.length == 1) {
//                        child1 = (BinaryTree) childrenArray[0];
//                        child2 = null;
//                    } else {
//                        child1 = null;
//                        child2 = null;
//                    }
//                }
//            } else {
//                throw new IllegalArgumentException(
//                    "The children Set must have no more than 2 elements.");
//            }
//        } else {
//            throw new IllegalArgumentException(
//                "The children Set must be not null.");
//        }
//    }
//
//    //we cannot add self or a (distant) parent
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param child1 DOCUMENT ME!
//     * @param child2 DOCUMENT ME!
//     *
//     * @throws CircularReferenceException DOCUMENT ME!
//     */
//    public void setChildren(BinaryTree child1, BinaryTree child2)
//        throws CircularReferenceException {
//        if (child1 != null) {
//            if ((child1.getParent() == null) && (child1 != this)) {
//                child1.setParent(this);
//                this.child1 = child1;
//            } else {
//                throw new CircularReferenceException(
//                    "Cannot add BinaryTree child that has already a parent.");
//            }
//        } else {
//            this.child1 = null;
//        }
//
//        if (child2 != null) {
//            if ((child2.getParent() == null) && (child2 != this)) {
//                child2.setParent(this);
//                this.child2 = child2;
//            } else {
//                throw new CircularReferenceException(
//                    "Cannot add BinaryTree child that has already a parent.");
//            }
//        } else {
//            this.child1 = null;
//        }
//    }
//
//    //we cannot add self or a (distant) parent
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param child DOCUMENT ME!
//     *
//     * @throws CircularReferenceException DOCUMENT ME!
//     * @throws IllegalArgumentException DOCUMENT ME!
//     */
//    public void setLeftChild(BinaryTree child)
//        throws CircularReferenceException {
//        //synchronized (this) {
//        if (child != null) {
//            if ((child != this) && (child.getParent() == null)) {
//                child.setParent(this);
//                this.child1 = child;
//            } else {
//                throw new CircularReferenceException(
//                    "Cannot add BinaryTree child that has already a parent.");
//            }
//        } else {
//            throw new IllegalArgumentException("The child must be not null.");
//        }
//
//        //}
//    }
//
//    //we cannot add self or a (distant) parent
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param child DOCUMENT ME!
//     *
//     * @throws CircularReferenceException DOCUMENT ME!
//     * @throws IllegalArgumentException DOCUMENT ME!
//     */
//    public void setRightChild(BinaryTree child)
//        throws CircularReferenceException {
//        //synchronized (this) {
//        if (child != null) {
//            if ((child != this) && (child.getParent() == null)) {
//                child.setParent(this);
//                this.child2 = child;
//            } else {
//                throw new CircularReferenceException(
//                    "Cannot add BinaryTree child that has already a parent.");
//            }
//        } else {
//            throw new IllegalArgumentException("The child must be not null.");
//        }
//
//        //}
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param child DOCUMENT ME!
//     *
//     * @throws IllegalArgumentException DOCUMENT ME!
//     */
//    public void removeChild(BinaryTree child) {
//        //synchronized (this) {
//        if (child != null) {
//            if (child.equals(child1)) {
//                child.setParent(null);
//                child1 = null;
//            } else if (child.equals(child2)) {
//                child.setParent(null);
//                child2 = null;
//            } else {
//                throw new IllegalArgumentException(
//                    "The child is not a child of this.");
//            }
//        } else {
//            throw new IllegalArgumentException("The child must be not null.");
//        }
//
//        //}
//    }
//
//    //isRoot() == !hasParent();
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public boolean hasParent() {
//        return parent != null;
//    }
//
//    //to set the parent, use element.getParent().removeChild(element);newParent.addChild(element);
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public BinaryTree getParent() {
//        return parent;
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param child DOCUMENT ME!
//     */
//    private void setParent(BinaryTree child) {
//        parent = child;
//    }
//
//    //0 if this element is the root
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public int getDepth() {
//        BinaryTree element;
//        int result;
//
//        element = this;
//        result = 0;
//
//        while (element.getParent() != null) {
//            result += 1;
//            element = element.getParent();
//        }
//
//        return result;
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public BinaryTree getRoot() {
//        BinaryTree element;
//        BinaryTree result;
//
//        element = this;
//        result = null;
//
//        while (element != null) {
//            result = element;
//            element = element.getParent();
//        }
//
//        return result;
//    }
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param tree DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public static BinaryTree getRoot(BinaryTree tree) {
//        BinaryTree element;
//        BinaryTree result;
//
//        result = null;
//
//        if (tree != null) {
//            element = tree;
//
//            while (element != null) {
//                result = element;
//                element = element.getParent();
//            }
//        }
//
//        return result;
//    }
//
//    //the common parent
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param tree1 DOCUMENT ME!
//     * @param tree2 DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public static BinaryTree getCommonRoot(BinaryTree tree1, BinaryTree tree2) {
//        Vector<BinaryTree>  lineage1;
//        Vector<BinaryTree>  lineage2;
//        BinaryTree root1;
//        BinaryTree root2;
//        BinaryTree result;
//        int i;
//
//        if ((tree1 != null) && (tree2 != null)) {
//            root1 = tree1.getRoot();
//            root2 = tree2.getRoot();
//
//            if ((root1 != null) && (root2 != null) && (root1.equals(root2))) {
//                lineage1 = getLineage(root1, tree1);
//                lineage2 = getLineage(root2, tree2);
//                i = 1;
//
//                while ((i < Math.min(lineage1.size(), lineage2.size())) &&
//                        lineage1.get(i).equals(lineage2.get(i))) {
//                    i++;
//                }
//
//                result = (BinaryTree) lineage1.get(i - 1);
//            } else {
//                result = null;
//            }
//        } else {
//            result = null;
//        }
//
//        return result;
//    }
//
//    //get all the children from tree1 to tree2
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param tree1 DOCUMENT ME!
//     * @param tree2 DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public static Vector<BinaryTree> getLineage(BinaryTree tree1, BinaryTree tree2) {
//        BinaryTree element;
//        Vector <BinaryTree> result;
//        boolean found;
//
//        if ((tree1 != null) && (tree2 != null)) {
//            result = new Vector<BinaryTree> ();
//            element = tree2;
//
//            while ((element != null) && !element.equals(tree1)) {
//                result.add(element);
//                element = element.getParent();
//            }
//
//            if (element == null) {
//                result = new Vector<BinaryTree> ();
//            } else {
//                result.add(element);
//            }
//        } else {
//            result = null;
//        }
//
//        return result;
//    }
//
//    //extracts the shortest and simpliest BinaryTree between two individuals
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param tree1 DOCUMENT ME!
//     * @param tree2 DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public static BinaryTree extractBinaryTree(BinaryTree tree1,
//        BinaryTree tree2) {
//        Vector<BinaryTree>  lineage1;
//        Vector<BinaryTree>  lineage2;
//        BinaryTree currentBinaryTree1;
//        BinaryTree nextBinaryTree1;
//        BinaryTree currentBinaryTree2;
//        BinaryTree nextBinaryTree2;
//        BinaryTree result;
//        int i;
//
//        if ((tree1 != null) && (tree2 != null)) {
//            lineage1 = getLineage(tree1.getRoot(), tree1);
//            lineage2 = getLineage(tree2.getRoot(), tree2);
//
//            if (lineage1.get(0).equals(lineage2.get(0))) {
//                result = null;
//
//                try {
//                    //find highest common element
//                    i = 1;
//
//                    if (lineage1.size() < lineage2.size()) {
//                        while ((i < lineage1.size()) &&
//                                lineage1.get(i).equals(lineage2.get(i))) {
//                            i++;
//                        }
//
//                        //build a new simple tree
//                        result = new BinaryTree(((BinaryTree) lineage1.get(i -
//                                    1)).getContents());
//                        currentBinaryTree1 = result;
//                        currentBinaryTree2 = result;
//
//                        while (i < lineage1.size()) {
//                            nextBinaryTree1 = new BinaryTree(((BinaryTree) lineage1.get(
//                                        i)).getContents());
//                            currentBinaryTree1.setLeftChild(nextBinaryTree1);
//                            currentBinaryTree1 = nextBinaryTree1;
//                            nextBinaryTree2 = new BinaryTree(((BinaryTree) lineage2.get(
//                                        i)).getContents());
//                            currentBinaryTree2.setRightChild(nextBinaryTree2);
//                            currentBinaryTree2 = nextBinaryTree2;
//                            i++;
//                        }
//
//                        while (i < lineage2.size()) {
//                            nextBinaryTree2 = new BinaryTree(((BinaryTree) lineage2.get(
//                                        i)).getContents());
//                            currentBinaryTree2.setRightChild(nextBinaryTree2);
//                            currentBinaryTree2 = nextBinaryTree2;
//                            i++;
//                        }
//                    } else {
//                        while ((i < lineage1.size()) &&
//                                lineage1.get(i).equals(lineage2.get(i))) {
//                            i++;
//                        }
//
//                        //build a new simple tree
//                        result = new BinaryTree(((BinaryTree) lineage1.get(i -
//                                    1)).getContents());
//                        currentBinaryTree1 = result;
//                        currentBinaryTree2 = result;
//
//                        while (i < lineage2.size()) {
//                            nextBinaryTree1 = new BinaryTree(((BinaryTree) lineage1.get(
//                                        i)).getContents());
//                            currentBinaryTree1.setLeftChild(nextBinaryTree1);
//                            currentBinaryTree1 = nextBinaryTree1;
//                            nextBinaryTree2 = new BinaryTree(((BinaryTree) lineage2.get(
//                                        i)).getContents());
//                            currentBinaryTree2.setRightChild(nextBinaryTree2);
//                            currentBinaryTree2 = nextBinaryTree2;
//                            i++;
//                        }
//
//                        while (i < lineage1.size()) {
//                            nextBinaryTree1 = new BinaryTree(((BinaryTree) lineage1.get(
//                                        i)).getContents());
//                            currentBinaryTree1.setLeftChild(nextBinaryTree1);
//                            currentBinaryTree1 = nextBinaryTree1;
//                            i++;
//                        }
//                    }
//                } catch (CircularReferenceException e) {
//                }
//            } else {
//                result = null;
//            }
//        } else {
//            result = null;
//        }
//
//        return result;
//    }
//
//    //return all children and grandchildren, etc...
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public Set<BinaryTree> getAllChildren() {
//        Iterator<BinaryTree>  iterator;
//        Set<BinaryTree>  result;
//
//        iterator = getChildren().iterator();
//        result = new HashSet<BinaryTree>(getChildren());
//
//        while (iterator.hasNext()) {
//            result.addAll(iterator.next().getAllChildren());
//        }
//
//        return result;
//    }
//
//    //perhaps this implementation is faster
//
//    /*
//    public Set getAllChildren() {
//
//        Set currentChildren;
//        HashSet newChildren;
//        Iterator iterator;
//        HashSet result;
//
//        currentChildren = getChildren();
//        newChildren = new HashSet();
//        result = new HashSet();
//        result.addAll(currentChildren);
//        iterator = currentChildren.iterator();
//        while (iterator.hasNext()) {
//            newChildren.addAll(((BinaryTree)iterator.next()).getChildren());
//        }
//        while (newChildren.size()!=0) {
//            result.addAll(newChildren);
//            currentChildren=newChildren;
//            newChildren = new HashSet();
//            iterator = currentChildren.iterator();
//            while (iterator.hasNext()) {
//                newChildren.addAll(((BinaryTree)iterator.next()).getChildren());
//            }
//        }
//
//        return result;
//
//    }
//     */
//    public boolean hasDistantChild(BinaryTree child) {
//        Set<BinaryTree>  currentChildren;
//        Set<BinaryTree>  tempChildren;
//        Set<BinaryTree>  newChildren;
//        Iterator<BinaryTree>  iterator;
//        boolean result;
//
//        if (child != null) {
//            currentChildren = getChildren();
//            newChildren = Collections.EMPTY_SET;
//            result = getChildren().contains(child);
//            iterator = currentChildren.iterator();
//
//            while (iterator.hasNext() && !result) {
//                tempChildren = iterator.next().getChildren();
//                result = tempChildren.contains(child);
//                newChildren.addAll(tempChildren);
//            }
//
//            while ((newChildren.size() != 0) && !result) {
//                currentChildren = newChildren;
//                newChildren = new HashSet<BinaryTree> ();
//                iterator = currentChildren.iterator();
//
//                while (iterator.hasNext() && !result) {
//                    tempChildren = iterator.next().getChildren();
//                    result = tempChildren.contains(child);
//                    newChildren.addAll(tempChildren);
//                }
//            }
//
//            return result;
//        } else {
//            throw new IllegalArgumentException("The child must be not null.");
//        }
//    }
//
//    //complete equality, of contents and relations
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @param o DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public boolean equals(Object o) {
//        BinaryTree tree;
//
//        if ((o != null) && (o instanceof BinaryTree)) {
//            tree = (BinaryTree) o;
//
//            return getContents().equals(tree.getContents()) &&
//            getLeftChild().equals(tree.getLeftChild()) &&
//            getRightChild().equals(tree.getRightChild()) &&
//            getParent().equals(tree.getParent());
//        } else {
//            return false;
//        }
//    }
//
//    //builds a shallow BinaryTree (new BinaryTree elements but same contents) from this element and under
//
//    /**
//     * DOCUMENT ME!
//     *
//     * @return DOCUMENT ME!
//     */
//    public Object clone() {
//        BinaryTree result;
//
//        //result = (BinaryTree) super.clone();
//        result = new BinaryTree(this.getContents());
//
//        //result.obj = getContents();
//        try {
//            result.setChildren(this.getLeftChild(), this.getRightChild());
//        } catch (CircularReferenceException e) {
//        }
//
//        return result;
//    }
//
//    public static void main(String arg){
//
//        BinaryTree bt=new BinaryTree();
//
//    }
//}
//
//
