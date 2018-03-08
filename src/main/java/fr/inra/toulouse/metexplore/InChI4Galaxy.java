package fr.inra.toulouse.metexplore;

import parsebionet.utils.chemicalStructures.InChI;
import parsebionet.utils.chemicalStructures.InchiLayer;

import java.util.Arrays;

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
       // System.out.println(Arrays.toString(boolInchiLayers));
        this.setLayers(boolInchiLayers);
    }

    private void setLayers(Boolean[] boolInchiLayers) {

        this.setValidInchi(true);
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

                tmp = trunckedInchi.split("/r");
                if (boolInchiLayers[8]) {
                    if (tmp.length != 2) {
                        this.setValidInchi(false);
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

}
