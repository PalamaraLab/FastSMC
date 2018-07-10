//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.


package ASMCprepareDecoding;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.net.URLDecoder;
import java.util.ArrayList;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.MutuallyExclusiveGroup;
import net.sourceforge.argparse4j.inf.Namespace;

public class ASMCprepareDecoding {

    final static String VERSION = "1.0";
    final static String VERSION_DATE = "July 1, 2018";
    final static String YEAR = "2018";
    final static String LICENSE = "GPL v3";
    final static String WEBSITE = "www.palamaralab.org/software/ASMC";
    final static String PROGRAM = "ASMC";

    public static void main(String[] args) throws IOException {

        Utils.print("");
        Utils.print(" █████╗   ███████╗  ███╗   ███╗   ██████╗");
        Utils.print("██╔══██╗  ██╔════╝  ████╗ ████║  ██╔════╝");
        Utils.print("███████║  ███████╗  ██╔████╔██║  ██║     ");
        Utils.print("██╔══██║  ╚════██║  ██║╚██╔╝██║  ██║     ");
        Utils.print("██║  ██║  ███████║  ██║ ╚═╝ ██║  ╚██████╗");
        Utils.print("╚═╝  ╚═╝  ╚══════╝  ╚═╝     ╚═╝   ╚═════╝");

        Utils.print("\nAscertained Sequentially Markovian Coalescent (ASMC) - Precompute decoding quantities v." + VERSION);
        Utils.print("GNU GPL v3, Copyright (C) " + YEAR + " Pier Palamara");
        Utils.print("Manual: " + WEBSITE + "\n");

        ArgumentParser parser = ArgumentParsers.newArgumentParser("ASMCprepareDecoding");
        parser.defaultHelp(true);

        parser.description("Ascertained Sequentially Markovian Coalescent (ASMC) - Precompute decoding quantities.");

        // mandatory: discretization (number of intervals or file)
        // mandatory: frequencies (frq file or haps file to read frequencies from)
        // mandatory: output file root
        // optional: demography (default: EU)
        // optional: CSFS file
        // optional: maximum CSFS size
        // optional: mutation rate
        //**********************************************************************
        // MANDATORY ARGUMENTS *************************************************
        //**********************************************************************
        // discretization. Either file with one time per line, or quantiles of coalescent or erlang distribution
        MutuallyExclusiveGroup groupDiscretization = parser.addMutuallyExclusiveGroup()
                .required(true);
        groupDiscretization.addArgument("-d", "--discretization")
                .metavar("Discretization")
                .type(String.class)
                .help("File with vector of time discretization intervals.");
        groupDiscretization.addArgument("-qCoal", "--coalescentQuantiles")
                .metavar("numberOfIntervals")
                .type(Integer.class)
                .help("Desired number of discretization intervals (quantiles from the pairwise coalescent distribution).");
//        groupDiscretization.addArgument("-qMut", "--mutationQuantiles")
//                .metavar("numberOfIntervals")
//                .type(Integer.class)
//                .help("Infer desired number of discretization intervals from the distribution of mutation age.");

        // to read array allele frequencies
        // root for data files to read frequencies
        MutuallyExclusiveGroup groupFrequency = parser.addMutuallyExclusiveGroup("Frequency")
                .required(true);
        groupFrequency.addArgument("-f", "--fileRoot")
                .metavar("FileRoot")
                .type(String.class)
                .help("Root for name of hap/samples files from which allele frequencies for array data should be read.");
        // frequency file for ascertainment (optional, use haps otherwise)
        groupFrequency.addArgument("-F", "--freqFile")
                .metavar("freqFile")
                .type(String.class)
                .help("Plink .frq file containing minor allele frequencies for array data.");

        // output root
        parser.addArgument("-o", "--out")
                .metavar("OutputFileRoot")
                .required(true)
                .type(String.class)
                .help("Output file root.");

        //**********************************************************************
        // OPTIONAL ARGUMENTS **************************************************
        //**********************************************************************
        // demographic model
        parser.addArgument("-D", "--demography")
                .metavar("Demography")
                .type(String.class)
                .help("File with demographic model. If not specified, the default EU model will be used.");

        // can input precomputed CSFS
        parser.addArgument("-C", "--CSFS")
                .metavar("CSFS")
                .type(String.class)
                .help("File with precomputed CSFS. If not specified this program will try to compute it internally.");

        // bound on number of total samples for CSFS if not those in haps
        parser.addArgument("-n", "--samples")
                .metavar("Samples")
                .type(Integer.class)
                .setDefault(300)
                .help("Number of total samples to be used for emission probability. If not specified will use min(300, sample size).");

        // mutation rate
        parser.addArgument("-mu", "--mut")
                .metavar("MutationRate")
                .type(Double.class)
                .setDefault(1.65e-8)
                .help("Mutation rate assumed by demographic model. If not specified, will assume mu = 1.65E-8.");

        //**********************************************************************
        // ********************* PARSE ARGUMENTS *******************************
        //**********************************************************************
        Namespace ns = null;
        try {
            ns = parser.parseArgs(args);
        } catch (ArgumentParserException e) {
            parser.handleError(e);
            System.exit(1);
//            Utils.exit("Error parsing input arguments.");
        }

        // PARSE OUTPUT INFO ***************************************************
        String outputFileRoot = ns.get("out");
        Utils.print("Will output files to: " + outputFileRoot + ".*");

        // PARSE DEMOGRAPHIC MODEL *********************************************
        ArrayList<Double> arrayTime = new ArrayList<Double>();
        ArrayList<Double> arraySize = new ArrayList<Double>();
        String demographicFile = ns.get("demography");
        if (!Utils.fileExists(demographicFile)) {
            // USE DEFAULT
            arrayTime = Transition.EUtime;
            arraySize = Transition.EUsize;
            Utils.print("Did not input a demographic model, using default EU model.");
//            Utils.print("Time: " + arrayTime);
//            Utils.print("Size: " + arraySize);
        } else {
            Utils.print("Will read demographic model from " + demographicFile + " ...");
            BufferedReader bDemo = Utils.openFile(demographicFile);
            String line;
            while ((line = bDemo.readLine()) != null) {
                String[] split = line.split("\\s+");
                double time = Double.parseDouble(split[0]);
                double size = Double.parseDouble(split[1]);
                arrayTime.add(time);
                arraySize.add(size);
            }
            arrayTime.add(Double.POSITIVE_INFINITY);
            arraySize.add(arraySize.get(arraySize.size() - 1));
//            Utils.print("Time vector was input from " + demographicFile + ": " + arrayTime);
//            Utils.print("Size vector was input from " + demographicFile + ": " + arraySize);
        }

        // PARSE DISCRETIZAITON ************************************************
        ArrayList<Double> arrayDiscretization = new ArrayList<Double>();
        String discretizationFile = ns.get("discretization");
        Integer coalescentIntervals = ns.get("coalescentQuantiles");
        Integer mutationAgeIntervals = ns.get("mutationQuantiles");
        if (discretizationFile != null) {
            Utils.print("Will read discretization intervals from " + discretizationFile + " ...");
            BufferedReader bDisc = Utils.openFile(discretizationFile);
            String line;
            while ((line = bDisc.readLine()) != null) {
                String[] split = line.split("\\s+");
                double disc = Double.parseDouble(split[0]);
                arrayDiscretization.add(disc);
            }
            arrayDiscretization.add(Double.POSITIVE_INFINITY);
//            Utils.print("Discretization vector was input: " + arrayDiscretization);
        } else {
            if (coalescentIntervals != null) {
                arrayDiscretization = Transition.getTimeExponentialQuantiles(coalescentIntervals, arrayTime, arraySize);
                arrayDiscretization.add(Double.POSITIVE_INFINITY);
//                Utils.print(coalescentIntervals + " discretization intervals from coalescent distribution: " + arrayDiscretization);
                Utils.print("Using " + coalescentIntervals + " discretization intervals from coalescent distribution.");
            } else if (mutationAgeIntervals != null) {
                arrayDiscretization = Transition.getTimeErlangQuantiles(mutationAgeIntervals, arrayTime, arraySize);
                arrayDiscretization.add(Double.POSITIVE_INFINITY);
//                Utils.print(mutationAgeIntervals + " discretization intervals from coalescent distribution: " + arrayDiscretization);
                Utils.print("Using " + mutationAgeIntervals + " discretization intervals from coalescent distribution.");
            }
        }

        // PARSE FREQUENCIES ***************************************************
        String fileRoot = ns.get("fileRoot");
        String freqFile = ns.get("freqFile");
        Data data = null;
        if (fileRoot != null) {
            Utils.print("Files will be read from: " + fileRoot + "*");
        }
        if (freqFile != null) {
            if (!Utils.fileExists(freqFile)) {
                Utils.exit("Could not open " + freqFile);
            }
            Utils.print("Will load minor allele frequencies from " + freqFile + " ...");
            data = new Data();
            data.addFreq(freqFile);
        } else {
            // user must have specified --fileRoot
            Utils.print("Will compute allele frequencies in haps files with root " + fileRoot + ".");
            if (fileRoot == null) {
                Utils.exit("Did not specify a valid file root.");
            }
            data = new Data(fileRoot);
        }

        // PARSE EMISSION OPTIONS **********************************************
        double mutRate = ns.get("mut");
        Utils.print("Will use mutation rate mu = " + mutRate + ".");

        // PARSE MAX CSFS SAMPLES **********************************************
        int samples = ns.get("samples");
        samples = (int) Math.min(samples, data.haploidSampleSize);
        Utils.print("Number of samples in CSFS calculations: " + samples + ".");

        // BUILD TRANSITION ****************************************************
        double[] timeArray = new double[arrayTime.size()];
        double[] sizeArray = new double[arraySize.size()];
        double[] discArray = new double[arrayDiscretization.size()];
        for (int i = 0; i < arrayTime.size(); i++) {
            timeArray[i] = arrayTime.get(i);
        }
        for (int i = 0; i < arraySize.size(); i++) {
            sizeArray[i] = arraySize.get(i);
        }
        for (int i = 0; i < arrayDiscretization.size(); i++) {
            discArray[i] = arrayDiscretization.get(i);
        }
        Transition transition = new Transition(timeArray, sizeArray, discArray, Transition.TransitionType.CSC);

//        // OUTPUT EXPECTED TIMES ***********************************************
//        Utils.printf("ExpectedTimes:");
//        for (int i = 0; i < transition.expectedTimes.length; i++) {
//            Utils.printf("\t" + transition.expectedTimes[i]);
//        }
        // PARSE CSFS **********************************************************
        String CSFSFile = ns.get("CSFS");
        CSFS csfs = null;
        // check that CSFS exists and contains all required entries.
        if (Utils.fileExists(CSFSFile)) {
            Utils.print("Will load precomputed CSFS from " + freqFile + " ...");
            csfs = CSFS.loadFromFile(CSFSFile);
            Utils.printf("Verifying CSFS loaded from " + CSFSFile + " ... ");
            if (!csfs.verify(arrayTime, arraySize, mutRate, samples, arrayDiscretization)) {
                Utils.print("Will compute new CSFS.");
                csfs = null;
            } else {
                Utils.print("Verified " + csfs.CSFS.size() + " CSFS entries.");
            }
        }
        if (csfs == null) {
            File jarFile = new File(ASMCprepareDecoding.class.getProtectionDomain().getCodeSource().getLocation().getPath());
            String decodedPath = URLDecoder.decode(jarFile.getParent(), "UTF-8");
            Utils.print("Generating and saving new CSFS.");
            if (demographicFile == null) {
                demographicFile = decodedPath + "/" + "CSFS_prepare.demo";
                Utils.print("Writing default EU demographic model in " + demographicFile);
                BufferedWriter bw = Utils.openFileForWriting(demographicFile);
                for (int i = 0; i < Transition.EUsize_array.length - 1; i++) {
                    bw.write(Transition.EUtime_array[i] + "\t" + Transition.EUsize_array[i] + "\n");
                }
                bw.close();
            }
            if (discretizationFile == null) {
                discretizationFile = decodedPath + "/" + "CSFS_prepare.disc";
                Utils.print("Writing discretization intervals in " + discretizationFile);
                BufferedWriter bw = Utils.openFileForWriting(discretizationFile);
                for (int i = 0; i < arrayDiscretization.size() - 1; i++) {
                    bw.write(arrayDiscretization.get(i) + "\n");
                }
                bw.close();
            }
            // generate and save.
            String[] command = new String[6];
            command[0] = "sh";
            command[1] = decodedPath + "/" + "getCSFS.sh";
            command[2] = demographicFile;
            command[3] = discretizationFile;
            command[4] = samples + "";
            command[5] = outputFileRoot;
            Utils.printf("Running: " + command[0]);
            for (int i = 1; i < command.length; i++) {
                Utils.printf(" " + command[i]);
            }
            Utils.print();
            try {
                ProcessBuilder pb = new ProcessBuilder(command);
                pb.inheritIO();
                Process p = pb.start();     // Start the process.
                p.waitFor();                // Wait for the process to finish.
                Utils.print("CSFS constructed.");
                CSFSFile = outputFileRoot + ".csfs";
                if (!Utils.fileExists(CSFSFile)) {
                    Utils.exit("Something went wrong. CSFS constructed, but CSFS file is not in " + CSFSFile + ".");
                }
            } catch (Exception e) {
                Utils.exit("Something went wrong when trying to build the CSFS. Verify above command, or build CSFS externally using the getCSFS.sh script and use --CSFS flag.");
//                e.printStackTrace();
            }
            csfs = CSFS.loadFromFile(CSFSFile);
            if (!csfs.verify(arrayTime, arraySize, mutRate, samples, arrayDiscretization)) {
                Utils.exit("Something went wrong when constructing the CSFS.");
            } else {
                Utils.print("Verified " + csfs.CSFS.size() + " CSFS entries.");
            }
        }
        csfs.fixAscertainment(data, samples, transition);

        // BUILD DECODING QUANTITIES *******************************************
        Utils.print("");
        Utils.print("Building decoding quantities...");
        DecodingQuantities decodingQuantities = new DecodingQuantities(csfs, transition, mutRate);
        decodingQuantities.saveDecodingQuantities(outputFileRoot);
        BufferedWriter bw = Utils.openFileForWriting(outputFileRoot + ".intervalsInfo");
        for (int i = 0; i < decodingQuantities.expectedTimes.length; i++) {
            bw.write(decodingQuantities.discretization[i] + "\t" + decodingQuantities.expectedTimes[i] + "\t" + decodingQuantities.discretization[i + 1] + "\n");
        }
        bw.close();
        Utils.print("Done\n");

    }

}
