package fr.inra.toulouse.metexplore;

import parsebionet.utils.chemicalStructures.InChI;
import parsebionet.utils.chemicalStructures.InchiLayer;

public class InChI4Galaxy extends InChI {

    public InChI4Galaxy() {
        this.validity = false;
    }

    public InChI4Galaxy(String inchiStringValue, String[] inchiLayers){
        this.inchiString = inchiStringValue;
        this.setBoolean(inchiLayers);
    }

    public InChI4Galaxy(String inchiStringValue) {
        this.inchiString = inchiStringValue;
        String[] inchiLayers = {"c","h"};
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
        this.setLayers(boolInchiLayers);
    }

    private void setLayers(Boolean[] boolInchiLayers) {

        String trunckedInchi = this.getInchiString();

        if (!trunckedInchi.startsWith("InChI=")) {
            this.setValidInchi(false);
        } else {
            trunckedInchi = trunckedInchi.substring(6);

            try {
                this.setVersion(Integer.parseInt(trunckedInchi.substring(0, 1)));
            } catch (Exception var4) {
                this.setValidInchi(false);
                return;
            }

            if (trunckedInchi.charAt(1) == 'S') {
                this.setStandard(true);
                trunckedInchi = trunckedInchi.substring(2);
            } else {
                this.setStandard(false);
                trunckedInchi = trunckedInchi.substring(1);
            }

            String[] tmp;
            InchiLayer connect;
            if (trunckedInchi.contains("/r")) {
                tmp = trunckedInchi.split("/r");
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    return;
                }

                if (boolInchiLayers[8]) {
                    connect = new InchiLayer('r', tmp[1]);
                    this.setReconnectedLayer(connect);
                }
                trunckedInchi = tmp[0];
            }

            if (trunckedInchi.contains("/f")) {
                tmp = trunckedInchi.split("/f");
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    return;
                }

                if (boolInchiLayers[7]) {
                    connect = new InchiLayer('f', tmp[1]);
                    this.setFixedLayer(connect);
                }
                trunckedInchi = tmp[0];
            }

            if (trunckedInchi.contains("/i")) {
                tmp = trunckedInchi.split("/i");
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    return;
                }

                if (boolInchiLayers[6]) {
                    connect = new InchiLayer('i', tmp[1]);
                    this.setIsotopicLayer(connect);
                }
                trunckedInchi = tmp[0];
            }

            if (trunckedInchi.contains("/t")) {
                tmp = trunckedInchi.split("/t");
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    return;
                }

                if (boolInchiLayers[5]) {
                    connect = new InchiLayer('t', tmp[1]);
                    this.setTetraStereoLayer(connect);
                }
                trunckedInchi = tmp[0];
            }

            if (trunckedInchi.contains("/b")) {
                tmp = trunckedInchi.split("/b");
                if (tmp.length != 2) {
                    this.setValidInchi(false);
                    return;
                }

                if (boolInchiLayers[4]) {
                    connect = new InchiLayer('b', tmp[1]);
                    this.setDbStereoLayer(connect);
                }
                trunckedInchi = tmp[0];
            }
            if (trunckedInchi.contains("/p")) {
                tmp = trunckedInchi.split("/p");
                if (tmp.length != 2 || tmp[1].contains("/")) {
                    this.setValidInchi(false);
                    return;
                }

                if (boolInchiLayers[3]) {
                    connect = new InchiLayer('p', tmp[1]);
                    this.setProtonationLayer(connect);
                }
                trunckedInchi = tmp[0];

            }

            if (trunckedInchi.contains("/q")) {
                tmp = trunckedInchi.split("/q");
                if (tmp.length != 2 || tmp[1].contains("/")) {
                    this.setValidInchi(false);
                    return;
                }


                if (boolInchiLayers[2]) {
                    connect = new InchiLayer('q', tmp[1]);
                    this.setChargeLayer(connect);
                }
                trunckedInchi = tmp[0];

            }

            if (trunckedInchi.contains("/h")) {
                tmp = trunckedInchi.split("/h");
                if (tmp.length != 2 || tmp[1].contains("/")) {
                    this.setValidInchi(false);
                    return;
                }


                if (boolInchiLayers[1]) {
                    connect = new InchiLayer('h', tmp[1]);
                    this.sethLayer(connect);
                }
                trunckedInchi = tmp[0];
            }

            if (trunckedInchi.contains("/c")) {
                tmp = trunckedInchi.split("/c");
                if (tmp.length != 2 || tmp[1].contains("/")) {
                    this.setValidInchi(false);
                    return;
                }

                if (boolInchiLayers[0]) {
                    connect = new InchiLayer('c', tmp[1]);
                    this.setConnectivity(connect);
                }
                trunckedInchi = tmp[0];
            }

            if (!trunckedInchi.equals("/")) {
                String formula = trunckedInchi.replaceAll("/", "");
                this.setFormula(formula);
            }
        }
    }

}
