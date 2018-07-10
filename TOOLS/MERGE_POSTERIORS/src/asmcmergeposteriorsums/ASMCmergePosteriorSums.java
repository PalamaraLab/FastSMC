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


package asmcmergeposteriorsums;

import java.io.IOException;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.MutuallyExclusiveGroup;
import net.sourceforge.argparse4j.inf.Namespace;

public class ASMCmergePosteriorSums {

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

        Utils.print("\nAscertained Sequentially Markovian Coalescent (ASMC) - merge posterior sums v." + VERSION + ", " + VERSION_DATE);
        Utils.print("GNU GPL v3, Copyright (C) " + YEAR + " Pier Palamara");
        Utils.print("Manual: " + WEBSITE + "\n");

        ArgumentParser parser = ArgumentParsers.newArgumentParser("MergePosteriorSums");

        parser.defaultHelp(true);

        parser.description("Merge ASMC posterior sums.");

        MutuallyExclusiveGroup input = parser.addMutuallyExclusiveGroup()
                .required(true);
        input.addArgument("--fileRoot")
                .metavar("root")
                .type(String.class)
                .help("Root of files containing posterior decodings");
        input.addArgument("--fileList")
                .metavar("root")
                .type(String.class)
                .help("File containing a list of file roots containing posterior decodings");

        parser.addArgument("--numJobs")
                .metavar("n")
                .type(Integer.class)
                .help("Total number of jobs that were run");

        parser.addArgument("--computeExpCoal")
                .metavar("intervalsInfoFile")
                .type(String.class)
                .help("Compute expected coalescent times. Requires a file containing intervals info");

        parser.addArgument("--norm")
                .metavar("normalize")
                .required(false)
                .action(Arguments.storeTrue())
                .help("Normalize output");

        parser.addArgument("--out")
                .metavar("fileName")
                .type(String.class)
                .help("Output file");

        /**
         * ********************* PARSE ARGUMENTS ***********************
         */
        Namespace ns = null;
        try {
            ns = parser.parseArgs(args);
        } catch (ArgumentParserException e) {
            parser.handleError(e);
            System.exit(1);
        }

        if (ns.get("fileRoot") != null && ns.get("numJobs") == null) {
            Utils.exit("fileRoot requires numJobs");
        }

        String intervalInfoFile = ns.get("computeExpectedCoal");
        boolean getExpectedCoal = intervalInfoFile != null;
        String outFileRoot = ns.get("out");
        boolean norm = ns.get("norm");

        PosteriorMerger pm = null;

        if (ns.get("fileRoot") != null) {
            String fileRoot = ns.get("fileRoot");
            if (outFileRoot == null) {
                outFileRoot = fileRoot;
            }
            int jobs = ns.get("numJobs");
            pm = new PosteriorMerger(fileRoot, jobs, outFileRoot, norm);
        }
        if (ns.get("fileList") != null) {
            String fileList = ns.get("fileList");
            if (outFileRoot == null) {
                outFileRoot = fileList;
            }
            pm = new PosteriorMerger(fileList, outFileRoot, norm);
        }

        pm.writeMatricesToFile(outFileRoot);

        if (getExpectedCoal) {
            pm.computeCoalescentTimes(intervalInfoFile);
            pm.writeExpCoalToFile(outFileRoot);
        }
    }

}
