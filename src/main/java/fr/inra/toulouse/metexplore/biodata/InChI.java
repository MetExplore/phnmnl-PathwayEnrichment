package fr.inra.toulouse.metexplore.biodata;

import metabolome.util.MassType;
import metabolome.util.MolecularFormula;
import parsebionet.utils.chemicalStructures.ChemicalStructure;
import parsebionet.utils.chemicalStructures.InchiLayer;

import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class InChI extends ChemicalStructure {
    public String wrongLayer = "";
    public String inchiString;
    public int version;
    public boolean standard;
    public String formulaLayer;
    public InchiLayer connectivity;
    public InchiLayer hLayer;
    public InchiLayer chargeLayer;
    public InchiLayer protonationLayer;
    public InchiLayer dbStereoLayer;
    public InchiLayer tetraStereoLayer;
    public InchiLayer isotopicLayer;
    public InchiLayer fixedLayer;
    public InchiLayer reconnectedLayer;

    public InChI() {
        this.validity = false;
    }

    public InChI(String inchiStringValue, String[] inchiLayers){
        this.inchiString = inchiStringValue;
        this.setBoolean(inchiLayers);
    }

    public void setBoolean (String[] inchiLayers){

        Boolean[] boolInchiLayers = {false,false,false,false,false,false,false,false,false};
        for (String i : inchiLayers) {
            if (i != null){
                switch (i) {
                    case "c":
                        boolInchiLayers[0] = true;
                        break;
                    case "h":
                        boolInchiLayers[1] = true;
                        break;
                    case "q":
                        boolInchiLayers[2] = true;
                        break;
                    case "p":
                        boolInchiLayers[3] = true;
                        break;
                    case "b":
                        boolInchiLayers[4] = true;
                        break;
                    case "t":
                        boolInchiLayers[5] = true;
                        break;
                    case "i":
                        boolInchiLayers[6] = true;
                        break;
                    case "f":
                        boolInchiLayers[7] = true;
                        break;
                    case "r":
                        boolInchiLayers[8] = true;
                        break;
                }
            }
        }
        // System.out.println(Arrays.toString(boolInchiLayers));
        this.setLayers(boolInchiLayers);
    }

    private void setLayers(Boolean[] boolInchiLayers) {

        this.setValidInchi(true);
        String trunckedInchi = this.getInchiString();

        if (!trunckedInchi.startsWith("InChI=")) {
            this.setValidInchi(false);
            this.wrongLayer = "InChI=";
        } else {
            trunckedInchi = trunckedInchi.substring(6);

            try {
                this.setVersion(Integer.parseInt(trunckedInchi.substring(0, 1)));
            } catch (Exception var4) {
                this.setValidInchi(false);
                this.wrongLayer="version";
                return;
            }

            if (trunckedInchi.charAt(1) == 'S') {
                this.setStandard(true);
                trunckedInchi = trunckedInchi.substring(2);
            } else {
                this.setStandard(false);
                trunckedInchi = trunckedInchi.substring(1);
                this.wrongLayer="S";
            }

            String[] tmp;
            InchiLayer connect;

            tmp = trunckedInchi.split("/r");
            if (boolInchiLayers[8]) {
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    this.wrongLayer="reconnected (r)";
                    return;
                }
                connect = new InchiLayer('r', tmp[1]);
                this.setReconnectedLayer(connect);
            }
            trunckedInchi = tmp[0];

            tmp = trunckedInchi.split("/f");
            if (boolInchiLayers[7]) {
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    this.wrongLayer="fixed (f)";
                    return;
                }
                connect = new InchiLayer('f', tmp[1]);
                this.setFixedLayer(connect);
            }
            trunckedInchi = tmp[0];

            tmp = trunckedInchi.split("/i");
            if (boolInchiLayers[6]) {
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    this.wrongLayer="isotopic (i)";
                    return;
                }
                connect = new InchiLayer('i', tmp[1]);
                this.setIsotopicLayer(connect);
            }
            trunckedInchi = tmp[0];

            tmp = trunckedInchi.split("/t");
            if (boolInchiLayers[5]) {
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    this.wrongLayer="tetraStereo (t)";
                    return;
                }
                connect = new InchiLayer('t', tmp[1]);
                this.setTetraStereoLayer(connect);
            }
            trunckedInchi = tmp[0];

            tmp = trunckedInchi.split("/b");
            if (boolInchiLayers[4]) {
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    this.wrongLayer="dbStereo (b)";
                    return;
                }
                connect = new InchiLayer('b', tmp[1]);
                this.setDbStereoLayer(connect);
            }
            trunckedInchi = tmp[0];

            tmp = trunckedInchi.split("/p");
            if (boolInchiLayers[3]) {
                if (tmp.length != 2 || tmp[1].contains("/")) {
                    this.setValidInchi(false);
                    this.wrongLayer="protonation (p)";
                    return;
                }
                connect = new InchiLayer('p', tmp[1]);
                this.setProtonationLayer(connect);
            }
            trunckedInchi = tmp[0];


            tmp = trunckedInchi.split("/q");
            if (boolInchiLayers[2]) {
                if (tmp.length != 2 || tmp[1].contains("/")) {
                    this.setValidInchi(false);
                    this.wrongLayer="charge (q)";
                    return;
                }
                connect = new InchiLayer('q', tmp[1]);
                this.setChargeLayer(connect);
            }
            trunckedInchi = tmp[0];


            tmp = trunckedInchi.split("/h");
            if (boolInchiLayers[1]) {
                if (tmp.length != 2 || tmp[1].contains("/")) {
                    this.setValidInchi(false);
                    this.wrongLayer="h";
                    return;
                }
                connect = new InchiLayer('h', tmp[1]);
                this.sethLayer(connect);
            }
            trunckedInchi = tmp[0];

            tmp = trunckedInchi.split("/c");
            if (boolInchiLayers[0]) {
                if (tmp.length != 2 || tmp[1].contains("/")) {
                    this.setValidInchi(false);
                    this.wrongLayer="connectivity (c)";
                    return;
                }
                connect = new InchiLayer('c', tmp[1]);
                this.setConnectivity(connect);
            }
            trunckedInchi = tmp[0];

            if (!trunckedInchi.equals("/")) {
                String formula = trunckedInchi.replaceAll("/", "");
                this.setFormula(formula);
            }
        }
    }

    public String getRealFormula() {
        if (this.getProtonationLayer() != null) {
            String protonState = this.getProtonationLayer().getValue();
            String inchiFormula = this.formulaLayer;
            String realFormula = inchiFormula;
            String sign = "";
            String value = "0";
            Matcher m = Pattern.compile("^([+-])(\\d+)$").matcher(protonState);
            if (m.matches()) {
                sign = m.group(1);
                value = m.group(2);
            }

            m = Pattern.compile(".*H(\\d*).*").matcher(inchiFormula);
            if (m.matches()) {
                String hnumber;
                if (m.group(1).isEmpty()) {
                    hnumber = "1";
                } else {
                    hnumber = m.group(1);
                }

                int realHNumber;
                if (sign.equals("+")) {
                    realHNumber = Integer.parseInt(hnumber) + Integer.parseInt(value);
                } else {
                    realHNumber = Integer.parseInt(hnumber) - Integer.parseInt(value);
                }

                if (realHNumber == 0) {
                    realFormula = inchiFormula.replaceAll("H" + m.group(1), "");
                } else {
                    realFormula = inchiFormula.replaceAll("H" + m.group(1), "H" + realHNumber);
                }
            } else if (sign.equals("+") && !value.equals("1")) {
                realFormula = inchiFormula + "H" + value;
            } else if (sign.equals("+")) {
                realFormula = inchiFormula + "H";
            }

            return realFormula;
        } else {
            return this.getFormulaLayer();
        }
    }

    public int getAbsoluteCharge() {
        int charge = 0;
        int protonation = 0;
        Matcher m;
        String chargeState;
        String chargeSign;
        String chargeValue;
        if (this.getProtonationLayer() != null) {
            chargeState = this.getProtonationLayer().getValue();
            chargeSign = "+";
            chargeValue = "0";
            m = Pattern.compile("^([+-])(\\d+)$").matcher(chargeState);
            if (m.matches()) {
                chargeSign = m.group(1);
                chargeValue = m.group(2);
            }

            if (chargeSign.equals("+")) {
                protonation = Integer.parseInt(chargeValue);
            } else {
                protonation = -Integer.parseInt(chargeValue);
            }
        }

        if (this.getChargeLayer() != null) {
            chargeState = this.getChargeLayer().getValue();
            chargeSign = "+";
            chargeValue = "0";
            m = Pattern.compile("^([+-])(\\d+)$").matcher(chargeState);
            if (m.matches()) {
                chargeSign = m.group(1);
                chargeValue = m.group(2);
            }

            if (chargeSign.equals("+")) {
                charge = Integer.parseInt(chargeValue);
            } else {
                charge = -Integer.parseInt(chargeValue);
            }
        }

        return charge + protonation;
    }

    public String computeChargedExactMass() {
        MolecularFormula molecularFormula = null;
        int charge = this.getAbsoluteCharge();

        try {
            molecularFormula = new MolecularFormula(this.getRealFormula());
        } catch (RuntimeException var9) {
            ;
        }

        if (molecularFormula != null) {
            double mass = molecularFormula.getMass(MassType.MONOISOTOPIC);
            mass *= 1.0E7D;
            mass = (double)Math.round(mass);
            mass /= 1.0E7D;
            double e = 5.486E-4D;
            double exact = mass - (double)charge * e;
            return Double.toString(exact);
        } else {
            return "0.0000000";
        }
    }

    public String computeAverageMass() {
        MolecularFormula molecularFormula = null;

        try {
            molecularFormula = new MolecularFormula(this.getFormulaLayer());
        } catch (RuntimeException var4) {
            ;
        }

        if (molecularFormula != null) {
            double mass = molecularFormula.getMass(MassType.MOLECULAR);
            mass *= 1.0E7D;
            mass = (double)Math.round(mass);
            mass /= 1.0E7D;
            return Double.toString(mass);
        } else {
            return "0.0000000";
        }
    }

    public void displayLayers() {
        System.out.println("initial String: " + this.inchiString);
        System.out.println("Version: " + this.getVersion());
        if (this.isStandard()) {
            System.out.println("Standard: True");
        } else {
            System.out.println("Standard: False");
        }

        System.out.println("Formula: " + this.getFormulaLayer());
        if (this.getConnectivity() != null) {
            System.out.println("connectivity: " + this.getConnectivity().toString());
        }

        if (this.gethLayer() != null) {
            System.out.println("hydrogen layer: " + this.gethLayer().toString());
        }

        if (this.getChargeLayer() != null) {
            System.out.println("charge layer: " + this.getChargeLayer().toString());
        }

        if (this.getProtonationLayer() != null) {
            System.out.println("protonation layer: " + this.getProtonationLayer().toString());
        }

        if (this.getDbStereoLayer() != null) {
            System.out.println("double bond layer: " + this.getDbStereoLayer().toString());
        }

        if (this.getTetraStereoLayer() != null) {
            System.out.println("tetrahedral stereo: " + this.getTetraStereoLayer().toString());
        }

        if (this.getIsotopicLayer() != null) {
            System.out.println("isotopic: " + this.getIsotopicLayer().toString());
        }

        if (this.getFixedLayer() != null) {
            System.out.println("fixed H: " + this.getFixedLayer().toString());
        }

        if (this.getReconnectedLayer() != null) {
            System.out.println("fixed H: " + this.getReconnectedLayer().toString());
        }

    }

    public boolean equals(InChI compared) {
        return this.equals(compared, "s");
    }

    public boolean equals(InChI compared, String Layers) {
        if (this.getVersion() != compared.getVersion()) {
            return false;
        } else if (!this.getFormulaLayer().equals(compared.getFormulaLayer())) {
            return false;
        } else if (this.getConnectivity() == null && compared.getConnectivity() != null || this.getConnectivity() != null && compared.getConnectivity() == null) {
            return false;
        } else if ((this.getConnectivity() != null || compared.getConnectivity() != null) && !this.getConnectivity().equals(compared.getConnectivity())) {
            return false;
        } else if ((this.gethLayer() != null || compared.gethLayer() == null) && (this.gethLayer() == null || compared.gethLayer() != null)) {
            if ((this.gethLayer() != null || compared.gethLayer() != null) && !this.gethLayer().equals(compared.gethLayer())) {
                return false;
            } else if (this.getProtonationLayer() == null && compared.getProtonationLayer() != null || this.getProtonationLayer() != null && compared.getProtonationLayer() == null) {
                return false;
            } else if ((this.getProtonationLayer() != null || compared.getProtonationLayer() != null) && !this.getProtonationLayer().equals(compared.getProtonationLayer())) {
                return false;
            } else if ((this.getChargeLayer() != null || compared.getChargeLayer() == null) && (this.getChargeLayer() == null || compared.getChargeLayer() != null)) {
                if ((this.getChargeLayer() != null || compared.getChargeLayer() != null) && !this.getChargeLayer().equals(compared.getChargeLayer())) {
                    return false;
                } else {
                    if (!Layers.isEmpty()) {
                        Set<Character> mySet = new HashSet();
                        char[] var4 = Layers.toCharArray();
                        int var5 = var4.length;

                        for(int var6 = 0; var6 < var5; ++var6) {
                            char c = var4[var6];
                            mySet.add(c);
                        }

                        if (mySet.contains('s')) {
                            if (this.getDbStereoLayer() == null && compared.getDbStereoLayer() != null || this.getDbStereoLayer() != null && compared.getDbStereoLayer() == null) {
                                return false;
                            }

                            if ((this.getDbStereoLayer() != null || compared.getDbStereoLayer() != null) && !this.getDbStereoLayer().equals(compared.getDbStereoLayer())) {
                                return false;
                            }

                            if (this.getTetraStereoLayer() == null && compared.getTetraStereoLayer() != null || this.getTetraStereoLayer() != null && compared.getTetraStereoLayer() == null) {
                                return false;
                            }

                            if ((this.getTetraStereoLayer() != null || compared.getTetraStereoLayer() != null) && !this.getTetraStereoLayer().equals(compared.getTetraStereoLayer())) {
                                return false;
                            }
                        }

                        if (mySet.contains('i')) {
                            if (this.getIsotopicLayer() == null && compared.getIsotopicLayer() != null || this.getIsotopicLayer() != null && compared.getIsotopicLayer() == null) {
                                return false;
                            }

                            if ((this.getIsotopicLayer() != null || compared.getIsotopicLayer() != null) && !this.getIsotopicLayer().equals(compared.getIsotopicLayer())) {
                                return false;
                            }
                        }

                        if (mySet.contains('f')) {
                            if (this.getFixedLayer() == null && compared.getFixedLayer() != null || this.getFixedLayer() != null && compared.getFixedLayer() == null) {
                                return false;
                            }

                            if ((this.getFixedLayer() != null || compared.getFixedLayer() != null) && !this.getFixedLayer().equals(compared.getFixedLayer())) {
                                return false;
                            }
                        }
                    }

                    return true;
                }
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    public boolean isValid() {
        return this.isValidInchi();
    }

    public String toString() {
        return this.inchiString;
    }

    public boolean isValidInchi() {
        return this.validity;
    }

    public String getInchiString() {
        return this.inchiString;
    }

    public int getVersion() {
        return this.version;
    }

    public boolean isStandard() {
        return this.standard;
    }

    public String getFormulaLayer() {
        return this.formulaLayer;
    }

    public InchiLayer getConnectivity() {
        return this.connectivity;
    }

    public InchiLayer gethLayer() {
        return this.hLayer;
    }

    public InchiLayer getChargeLayer() {
        return this.chargeLayer;
    }

    public InchiLayer getProtonationLayer() {
        return this.protonationLayer;
    }

    public InchiLayer getDbStereoLayer() {
        return this.dbStereoLayer;
    }

    public InchiLayer getTetraStereoLayer() {
        return this.tetraStereoLayer;
    }

    public InchiLayer getIsotopicLayer() {
        return this.isotopicLayer;
    }

    public InchiLayer getFixedLayer() {
        return this.fixedLayer;
    }

    public InchiLayer getReconnectedLayer() {
        return this.reconnectedLayer;
    }

    public void setValidInchi(boolean validInchi) {
        this.validity = validInchi;
    }

    public void setInchiString(String inchiString) {
        this.inchiString = inchiString;
    }

    public void setVersion(int version) {
        this.version = version;
    }

    public void setStandard(boolean standard) {
        this.standard = standard;
    }

    public void setFormula(String formula) {
        this.formulaLayer = formula;
    }

    public void setConnectivity(InchiLayer connectivity) {
        this.connectivity = connectivity;
    }

    public void sethLayer(InchiLayer hLayer) {
        this.hLayer = hLayer;
    }

    public void setChargeLayer(InchiLayer chargeLayer) {
        this.chargeLayer = chargeLayer;
    }

    public void setProtonationLayer(InchiLayer protonationLayer) {
        this.protonationLayer = protonationLayer;
    }

    public void setDbStereoLayer(InchiLayer dbStereoLayer) {
        this.dbStereoLayer = dbStereoLayer;
    }

    public void setTetraStereoLayer(InchiLayer tetraStereoLayer) {
        this.tetraStereoLayer = tetraStereoLayer;
    }

    public void setIsotopicLayer(InchiLayer isotopicLayer) {
        this.isotopicLayer = isotopicLayer;
    }

    public void setFixedLayer(InchiLayer fixedLayer) {
        this.fixedLayer = fixedLayer;
    }

    public void setReconnectedLayer(InchiLayer reconnectedLayer) {
        this.reconnectedLayer = reconnectedLayer;
    }
}