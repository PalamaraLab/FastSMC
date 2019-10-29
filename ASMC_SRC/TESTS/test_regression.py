import unittest
import numpy as np
from pathlib import Path
from asmc import (
    HMM,
    DecodingQuantities,
    DecodingParams,
    Data,
)


class TestASMCRegression(unittest.TestCase):

    def setUp(self):
        hapsFileRoot = "FILES/EXAMPLE/exampleFile.n300.array"
        decodingQuantFile = "FILES/DECODING_QUANTITIES" \
            "/30-100-2000.decodingQuantities.gz"
        self.sequenceLength = Data.countHapLines(hapsFileRoot)
        self.params = DecodingParams(hapsFileRoot, decodingQuantFile,
                                     doPosteriorSums=True)
        self.decodingQuantities = DecodingQuantities(decodingQuantFile)
        self.data = Data(hapsFileRoot, self.sequenceLength,
                         self.decodingQuantities.CSFSSamples,
                         self.params.foldData, self.params.usingCSFS)
        self.hmm = HMM(self.data, self.decodingQuantities, self.params,
                       not self.params.noBatches, 1)

    def test_regression(self):
        old_sumOverPairs = np.loadtxt(Path(__file__).parent / 'data' /
                                      'regression_test_original.gz')
        self.hmm.decodeAll(self.params.jobs, self.params.jobInd)
        ret = self.hmm.getDecodingReturnValues()
        self.assertEqual(np.allclose(ret.sumOverPairs, old_sumOverPairs), True)
