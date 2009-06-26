/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.vflib.query;

import chemlib.ebi.vflib.interfaces.IAtomMatcher;
import chemlib.ebi.vflib.interfaces.IQuery;
import chemlib.ebi.vflib.interfaces.IQueryCompiler;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
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
public class TemplateCompiler implements IQueryCompiler {

//    private Reducer reducer;
//    private Map<IAtom, Integer> reductions;
    private IAtomContainer molecule;

    public TemplateCompiler() {
//        reducer = new Reducer();
//        reductions = new HashMap<IAtom, Integer>();
    }

    public static IQuery compile(IAtomContainer molecule) {
        TemplateCompiler compiler = new TemplateCompiler();

        compiler.setMolecule(molecule);

        return compiler.compile();
    }

    public void setMolecule(IAtomContainer molecule) {
        this.molecule = molecule;
    }

    public IAtomContainer getMolecule() {
        return molecule;
    }

    @Override
    public IQuery compile() {
//        try {
//            IAtomContainer copy = (IAtomContainer) queryMolecule.clone();
//        reductions.clear();
        //Remove Hydrogen if Necesarry Asad
//            reducer.reduce(copy, reductions);
        return build(molecule);
//        } catch (CloneNotSupportedException ex) {
//            Logger.getLogger(TemplateCompiler.class.getName()).log(Level.SEVERE, null, ex);
//        }

//        return null;
    }

    private IQuery build(IAtomContainer queryMolecule) {
        VFQuery result = new VFQuery();

        for (int i = 0; i < queryMolecule.getAtomCount(); i++) {
            IAtom atom = queryMolecule.getAtom(i);
            IAtomMatcher matcher = createMatcher(atom);

            if (matcher != null) {
                result.addNode(matcher, atom);
            }
        }

        for (int i = 0; i < queryMolecule.getBondCount(); i++) {
            IBond bond = queryMolecule.getBond(i);

            IAtom s = bond.getAtom(0);
            IAtom t = bond.getAtom(1);

            int IndexI = queryMolecule.getAtomNumber(s);
            int IndexJ = queryMolecule.getAtomNumber(t);
            result.connect(result.getNode(IndexI), result.getNode(IndexJ), new VFBondMatcher(bond));
        }

        return result;
    }

    private IAtomMatcher createMatcher(IAtom atom) {
        return new VFAtomMatcher(atom);
    }
}
