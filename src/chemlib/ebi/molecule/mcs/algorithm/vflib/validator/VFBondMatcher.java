/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.molecule.mcs.algorithm.vflib.validator;

import chemlib.ebi.molecule.mcs.global.BondType;
import chemlib.ebi.molecule.mcs.algorithm.vflib.interfaces.IBondMatcher;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;


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
public class VFBondMatcher implements IBondMatcher {

    private IBond queryBond;
    private int unsaturation;
    
    public VFBondMatcher() {
        this.queryBond = null;
        this.unsaturation = -1;
    }

    public VFBondMatcher(IBond queryBond) {
        this.queryBond = queryBond;
        this.unsaturation = getUnsaturation(queryBond);
    }

    /**
     * 
     * @param targetBond
     * @return
     */
    @Override
    public boolean matches(IBond targetBond) {

        boolean bondMatch = BondType.getInstance().getBondSensitiveFlag();

        if (!bondMatch) {
            return true;
        }

        if ((queryBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC)) && (queryBond.getOrder().equals(targetBond.getOrder()))) {
            return true;
        } else if (queryBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
            return true;

        }

       if (this.unsaturation == getUnsaturation(targetBond)) {
            return true;
        }

        return false;
    }

    private int getUnsaturation(IBond bond) {
        return getUnsaturation(bond.getAtom(0)) + getUnsaturation(bond.getAtom(1));
    }

    private int getUnsaturation(IAtom atom) {
        return atom.getValency() - atom.getFormalNeighbourCount();
    }
}
