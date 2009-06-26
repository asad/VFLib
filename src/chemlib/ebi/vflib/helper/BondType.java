/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemlib.ebi.vflib.helper;

/**
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BondType {

    private static BondType INSTANCE = null;
    private boolean bondtype = false;

    /**
     *
     * @return
     */
    public static synchronized BondType getInstance() {
        if (INSTANCE == null) {

            // it's ok, we can call this constructor
            INSTANCE = new BondType();
        }

        return INSTANCE;
    }

    protected BondType() {
    }

    /**
     *
     * @param isBondSensetive (set true if bondsensetive else false)
     */
    public void setBondSensitiveFlag(boolean isBondSensetive) {

        this.bondtype = isBondSensetive;
    }

    /**
     *
     * @return true if bondsensetive else false
     */
    public boolean getBondSensitiveFlag() {

        return bondtype;
    }
}
