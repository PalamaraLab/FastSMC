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

import java.util.HashMap;
import org.apache.commons.math3.distribution.HypergeometricDistribution;

public class ArraySpectrum {

    // spectrum does not include probability of monomorphic alleles due to no variation or subsampling
    double[] spectrum;
    // this is computed externally given population spectrum from demographic model
    double[] samplingFactors;
    // probability of monomorphic is stored separately
    double monomorphic;

    // in: data - contains allele frequencies, either computed from haps file, or from freq file.
    // in: samples - is bound on number of samples in csfs
    //
    // builds a spectrum for array data (using hypergeometric subsampling given
    // frequencies if samples is less than what was used to compute frequencies)
    //
    public ArraySpectrum(Data data, int samples) {

        HashMap<Double, HypergeometricDistribution> distributions = new HashMap<Double, HypergeometricDistribution>();
        HashMap<Double, Integer> distCounts = new HashMap<Double, Integer>();

        // get hypergeometric for each allele
        int monoMorphic = 0;
        for (int i = 0; i < data.allSNPsAlleleCounts.size(); i++) {
            int popSize = data.allSNPsAlleleCounts.get(i);
            int minorAlleles = data.allSNPsMinorAlleles.get(i);
            double freq = data.allSNPsFreq.get(i);
            if (minorAlleles == 0) {
                monoMorphic++;
                continue;
            }
            if (!distributions.containsKey(freq)) {
                HypergeometricDistribution dist = new HypergeometricDistribution(popSize, minorAlleles, samples);
                distributions.put(freq, dist);
                distCounts.put(freq, 1);
            } else {
                distCounts.put(freq, distCounts.get(freq) + 1);
            }
        }

        // sum all hypergeometrics to get spectrum in subsample (exclude [0] and [samples], which are monomorphic.
        spectrum = new double[samples + 1]; // include monomorphic at 0 and samples
        for (double d : distributions.keySet()) {
            int c = distCounts.get(d);
            HypergeometricDistribution dist = distributions.get(d);
            for (int i = 0; i <= samples; i++) {
                spectrum[i] += dist.probability(i) * c;
            }
        }
        // the term in 0 will contain monomorphic samples that are either present in the data, or due to subsampling
        spectrum[0] += monoMorphic;

        // normalize including monomorphic at 0 and samples
        spectrum = normalize(spectrum);

        // store monomorphic probability separately, then renormalize excluding monomorphic at 0 and samples
        // add alleles for which all samples are carriers to monomorphic probability
        monomorphic = spectrum[0] + spectrum[samples];
        Utils.print("Probability of a site being monomorphic due to subsampling: " + Math.round(monomorphic*1000)/1000. + ".");
        // renormalize without monomorphic
        spectrum[0] = 0.;
        spectrum[samples] = 0.;
        spectrum = normalize(spectrum);

        // fold to minor allele
        int halfTotal = samples / 2;
        double[] foldedAFS = new double[halfTotal + 1];
        for (int i = 0; i < halfTotal; i++) {
            foldedAFS[i] = spectrum[i] + spectrum[samples - i];
        }
        foldedAFS[halfTotal] = spectrum[halfTotal];
        spectrum = foldedAFS;

    }

    // normalize spectrum to 1
    private double[] normalize(double[] spectrum) {
        double[] normalized = new double[spectrum.length];
        double tot = 0.;
        for (int i = 0; i < spectrum.length; i++) {
            tot += spectrum[i];
        }
        for (int i = 0; i < spectrum.length; i++) {
            normalized[i] = spectrum[i] / tot;
        }
        return normalized;
    }

}
