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
import java.io.IOException;
import java.util.ArrayList;

public class Data {

    // holding SNP info
    ArrayList<Double> allSNPsFreq = new ArrayList<Double>();
    ArrayList<Integer> allSNPsMinorAlleles = new ArrayList<Integer>();
    ArrayList<Integer> allSNPsAlleleCounts = new ArrayList<Integer>();

    // some variables about the data
    String chr;
    int haploidSampleSize = 0;
    
    public Data() {}

    // constructor: no freq file was specified. Will look for one, or use haps file.
    public Data(String hapsFileRoot) throws IOException {
        // look for frq file
        if (Utils.fileExists(hapsFileRoot + ".frq.gz")) {
            readMinorAlleleFrequencies(hapsFileRoot + ".frq.gz");
        } else if (Utils.fileExists(hapsFileRoot + ".frq")) {
            readMinorAlleleFrequencies(hapsFileRoot + ".frq");
        } else {
            // or compute it from haps
            computeMinorAlleleFrequenciesFromHaps(hapsFileRoot);
        }
    }

    // constructor. Allele frequencies are specified in freq file. No need to compute them.
    public void addFreq(String freqFile) throws IOException {
        // look for frq file
        readMinorAlleleFrequencies(freqFile);

    }

    // read minor allele frequencies using freq file
    private void readMinorAlleleFrequencies(String freqFile) {
        try {
            BufferedReader br = Utils.openFile(freqFile);
            // skip header
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                String[] strSplit = line.split("\\s+");
                double freq = Double.parseDouble(strSplit[5]);
                Double popSize = Double.parseDouble(strSplit[6]);
                if (popSize.intValue() > this.haploidSampleSize) {
                    this.haploidSampleSize = popSize.intValue();
                }
                Double minorAlleles = popSize * freq;
                allSNPsFreq.add(freq);
                allSNPsMinorAlleles.add(minorAlleles.intValue());
                allSNPsAlleleCounts.add(popSize.intValue());
            }
        } catch (IOException ex) {
            Utils.exit("Could not read frq file.");
        }
    }

    // read minor allele frequencies using haps file
    private void computeMinorAlleleFrequenciesFromHaps(String hapsFileRoot) throws IOException {
        BufferedReader hapsBr = null;
        if (Utils.fileExists(hapsFileRoot + ".hap.gz")) {
            hapsBr = Utils.openFile(hapsFileRoot + ".hap.gz");
        } else if (Utils.fileExists(hapsFileRoot + ".hap")) {
            hapsBr = Utils.openFile(hapsFileRoot + ".hap");
        } else if (Utils.fileExists(hapsFileRoot + ".haps.gz")) {
            hapsBr = Utils.openFile(hapsFileRoot + ".haps.gz");
        } else if (Utils.fileExists(hapsFileRoot + ".haps")) {
            hapsBr = Utils.openFile(hapsFileRoot + ".haps");
        } else {
            Utils.exit("Could not find hap file in " + hapsFileRoot + ".hap.gz, " + hapsFileRoot + ".hap, " + ".haps.gz, or " + hapsFileRoot + ".haps");
        }
        String line;
        int pos = 0, monomorphic = 0;
        while ((line = hapsBr.readLine()) != null) {
            int DAcount = 0;
            String[] splitStr = line.split("\\s+");
            for (int i = 5; i < splitStr.length; i++) {
                if (splitStr[i].compareToIgnoreCase("1") == 0) {
                    DAcount++;
                }
            }
            int samples = splitStr.length - 5;
            if (samples > this.haploidSampleSize) {
                this.haploidSampleSize = samples;
            }
            if (samples % 2 != 0) {
                Utils.exit("Haps file contains a line with odd haploid sample size.");
            }
            if (DAcount > samples / 2) {
                DAcount = samples - DAcount;
            }
            double DAfreq = DAcount/(double)samples;
            double MAfreq = Math.min(DAfreq, 1-DAfreq);
            allSNPsFreq.add(MAfreq);
            allSNPsMinorAlleles.add(DAcount);
            allSNPsAlleleCounts.add(samples);
            if (DAcount == 0) {
                monomorphic++;
            }
            pos++;
        }
        hapsBr.close();
        Utils.print("Computed frequencies for " + haploidSampleSize + " haploid samples and " + pos + " markers, " + monomorphic + " of which are monomorphic.");
    }

}
