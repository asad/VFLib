/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.mcs.algorithm.vflib.validator;

import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IAtomMatcher;
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
public class VFAtomMatcher implements IAtomMatcher {

    private String symbol;
    private int maximumNeighbors;

    public VFAtomMatcher() {
        symbol = null;
        maximumNeighbors = -1;
    }

    /**
     * 
     * @param atom
     */
    public VFAtomMatcher(IAtom atom) {
        this();

        this.symbol = atom.getSymbol();
    }

    /**
     *
     * @param template
     * @param blockedPositions
     */
    public VFAtomMatcher(IAtom template, int blockedPositions) {
        this(template);
        Integer hCount = template.getHydrogenCount();
        if (hCount != null) {
            this.maximumNeighbors = template.getFormalNeighbourCount() + hCount;
        } else {
            this.maximumNeighbors = template.getFormalNeighbourCount();
        }
         this.maximumNeighbors -= blockedPositions;
    }

    /**
     *
     * @param atom
     * @return
     */
    @Override
    public boolean matches(IAtom atom) {
        if (!matchSymbol(atom)) {
            return false;
        }

        if (!matchMaximumNeighbors(atom)) {
            return false;
        }

        return true;
    }

    public void setMaximumNeighbors(int maximum) {
        this.maximumNeighbors = maximum;
    }

    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    private boolean matchSymbol(IAtom atom) {
        if (symbol == null) {
            return true;
        }

        return symbol.equalsIgnoreCase(atom.getSymbol());
    }

    private boolean matchMaximumNeighbors(IAtom atom) {
        if (maximumNeighbors == -1) {
            return true;
        }

        return atom.getFormalNeighbourCount() <= maximumNeighbors;
    }
}
//
//    private String symbol;
//    private int maximumNeighbors;
//    private int minimumNeighbors;
//    private int minimumValence;
//    private int maximumValence;
//
//    public VFAtomMatcher() {
//        symbol = null;
//        maximumNeighbors = -1;
//        minimumNeighbors = -1;
//        minimumValence = -1;
//        maximumValence = -1;
//    }
//
//    public VFAtomMatcher(IAtom atom) {
//        this();
//
//        this.symbol = atom.getSymbol();
//        this.minimumNeighbors = atom.getFormalNeighbourCount();
//        Integer hCount = atom.getHydrogenCount();
//        if (hCount != null) {
//            this.minimumValence = atom.getFormalNeighbourCount() + atom.getHydrogenCount();
//        } else {
//            this.minimumValence = atom.getFormalNeighbourCount();
//        }
//
////        System.out.println("symbol:" + symbol);
////        System.out.println("minimumNeighbors:" + minimumNeighbors);
////        System.out.println("minimumValence:" + minimumValence);
//    }
//
//    /**
//     *
//     * @param atom
//     * @return
//     */
//    
//    public boolean matches(IAtom atom) {
//        if (!matchSymbol(atom)) {
//            return false;
//        }
//
//        if (!matchMaximumNeighbors(atom)) {
//            return false;
//        }
//
//        if (!matchMinimumNeighbors(atom)) {
//            return false;
//        }
//
//        if (!matchMinimumValence(atom)) {
//            return false;
//        }
//
//        if (!matchMaximumValence(atom)) {
//            return false;
//        }
//
//        return true;
//    }
//
//    public void setMinimumValence(int minimum) {
//        if (minimum > maximumValence && maximumValence != -1) {
//            throw new IllegalStateException("Minimum " + minimum + " exceeds maximum");
//        }
//        this.minimumValence = minimum;
//    }
//
//    public void setMaximumValence(int maximum) {
//        if (maximum < minimumValence) {
//            throw new IllegalStateException("Maximum " + maximum + " less than minimum");
//        }
//        this.maximumValence = maximum;
//    }
//
//    public void setMaximumNeighbors(int maximum) {
//        if (maximum < minimumNeighbors) {
//            throw new IllegalStateException("Maximum " + maximum + " exceeds minimum " + minimumNeighbors);
//        }
//
//        this.maximumNeighbors = maximum;
//    }
//
//    public void setMinimumNeighbors(int minimum) {
//        if (minimum > maximumNeighbors && maximumNeighbors != -1) {
//            throw new IllegalStateException("Minimum " + minimum + " exceeds maximum " + maximumNeighbors);
//        }
//
//        this.minimumNeighbors = minimum;
//    }
//
//    public void setSymbol(String symbol) {
//        this.symbol = symbol;
//    }
//
//    private boolean matchSymbol(IAtom atom) {
//        if (symbol == null) {
//            return true;
//        }
//
//        return symbol.equalsIgnoreCase(atom.getSymbol());
//    }
//
//    private boolean matchMaximumNeighbors(IAtom atom) {
//        if (maximumNeighbors == -1) {
//            return true;
//        }
//
//        return atom.getFormalNeighbourCount() <= maximumNeighbors;
//    }
//
//    private boolean matchMinimumNeighbors(IAtom atom) {
//        if (minimumNeighbors == -1) {
//            return true;
//        }
//
//        return atom.getFormalNeighbourCount() >= minimumNeighbors;
//    }
//
//    private boolean matchMinimumValence(IAtom atom) {
//        if (minimumValence == -1) {
//            return true;
//        }
//
//        Integer hCount = atom.getHydrogenCount();
//        if (hCount != null) {
//            return atom.getFormalNeighbourCount() + hCount >= minimumValence;
//        } else {
//            return atom.getFormalNeighbourCount() >= minimumValence;
//        }
//
//    }
//
//    private boolean matchMaximumValence(IAtom atom) {
//        if (maximumValence == -1) {
//            return true;
//        }
//
//        return atom.getFormalNeighbourCount() + atom.getHydrogenCount() <= maximumValence;
//    }
//}
//
