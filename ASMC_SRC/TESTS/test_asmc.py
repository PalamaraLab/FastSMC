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

    def test_decode_pair(self):
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)
        self.hmm.decodePair(0, 9)
        self.assertEqual(len(self.hmm.getBatchBuffer()), 4)
        for i in range(4):
            self.assertEqual(self.hmm.getBatchBuffer()[0].iName,
                             self.data.individuals[0].name)
            self.assertEqual(self.hmm.getBatchBuffer()[0].jName,
                             self.data.individuals[9].name)
        self.hmm.decodePair(1, 1)
        self.assertEqual(len(self.hmm.getBatchBuffer()), 5)
        self.assertEqual(self.hmm.getBatchBuffer()[4].iName,
                         self.data.individuals[1].name)
        self.assertEqual(self.hmm.getBatchBuffer()[4].jName,
                         self.data.individuals[1].name)

    def test_decode_pairs(self):
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)
        self.hmm.decodePairs([0, 1], [9, 1])
        self.assertEqual(len(self.hmm.getBatchBuffer()), 5)
        for i in range(4):
            self.assertEqual(self.hmm.getBatchBuffer()[0].iName,
                             self.data.individuals[0].name)
            self.assertEqual(self.hmm.getBatchBuffer()[0].jName,
                             self.data.individuals[9].name)
        self.assertEqual(self.hmm.getBatchBuffer()[4].iName,
                         self.data.individuals[1].name)
        self.assertEqual(self.hmm.getBatchBuffer()[4].jName,
                         self.data.individuals[1].name)

    def test_finish_decoding(self):
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)
        self.hmm.decodePair(0, 9)
        self.assertEqual(len(self.hmm.getBatchBuffer()), 4)
        self.hmm.finishDecoding()
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)

    def test_fill_up_buffer(self):
        for i in range(1, (64 // 4) + 1):
            self.hmm.decodePair(0, i)
        # buffer should be empty now
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)


if __name__ == "__main__":
    unittest.main()
