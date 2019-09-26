import unittest
from asmc import (
    HMM,
    DecodingQuantities,
    DecodingParams,
    Data
)


class TestASMC(unittest.TestCase):

    def setUp(self):
        hapsFileRoot = "FILES/EXAMPLE/exampleFile.n300.array"
        decodingQuantFile = "FILES/DECODING_QUANTITIES" \
            "/30-100-2000.decodingQuantities.gz"
        sequenceLength = Data.countHapLines(hapsFileRoot)
        params = DecodingParams(hapsFileRoot, decodingQuantFile)
        decodingQuantities = DecodingQuantities(decodingQuantFile)
        self.data = Data(hapsFileRoot, sequenceLength,
                         decodingQuantities.CSFSSamples, params.foldData,
                         params.usingCSFS)
        self.hmm = HMM(self.data, decodingQuantities, params,
                       not params.noBatches, 1)

    def test_initialization(self):
        self.assertGreater(len(self.data.individuals), 20)


if __name__ == "__main__":
    unittest.main()
