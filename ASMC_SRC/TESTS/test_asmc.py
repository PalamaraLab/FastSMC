import unittest
from asmc import (
    HMM,
    DecodingQuantities,
    DecodingParams,
    Data,
    makePairObs
)


class TestASMC(unittest.TestCase):

    def setUp(self):
        hapsFileRoot = "FILES/EXAMPLE/exampleFile.n300.array"
        decodingQuantFile = "FILES/DECODING_QUANTITIES" \
            "/30-100-2000.decodingQuantities.gz"
        self.sequenceLength = Data.countHapLines(hapsFileRoot)
        self.params = DecodingParams(hapsFileRoot, decodingQuantFile)
        self.decodingQuantities = DecodingQuantities(decodingQuantFile)
        self.data = Data(hapsFileRoot, self.sequenceLength,
                         self.decodingQuantities.CSFSSamples, self.params.foldData,
                         self.params.usingCSFS)
        self.hmm = HMM(self.data, self.decodingQuantities, self.params,
                       not self.params.noBatches, 1)

    def tearDown(self):
        del self.hmm

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

    def test_decode_pair_observation(self):
        self.assertEqual(len(self.decodingQuantities.discretization),
                         len(self.decodingQuantities.expectedTimes) + 1)
        self.assertEqual(self.data.sites, self.sequenceLength)

        for p in [
            makePairObs(self.data.individuals[0], 1, self.data.individuals[0], 2),
            makePairObs(self.data.individuals[0], 1, self.data.individuals[1], 1),
            makePairObs(self.data.individuals[0], 2, self.data.individuals[1], 2)]:
            d = self.hmm.decode(p)
            self.assertEqual(len(d), len(self.decodingQuantities.expectedTimes))
            for i in range(len(d)):
                self.assertEqual(len(d[i]), self.data.sites)

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

class TestASMCDecodingParams(unittest.TestCase):
    def test_no_compress(self):
        hapsFileRoot = "FILES/EXAMPLE/exampleFile.n300.array"
        decodingQuantFile = "FILES/DECODING_QUANTITIES" \
            "/30-100-2000.decodingQuantities.gz"
        sequenceLength = Data.countHapLines(hapsFileRoot)
        params = DecodingParams(hapsFileRoot, decodingQuantFile, compress=True,
            skipCSFSdistance=float('nan'))

        self.assertEqual(params.compress, True)
        self.assertEqual(params.skipCSFSdistance, float('inf'))

        decodingQuantities = DecodingQuantities(decodingQuantFile)
        data = Data(hapsFileRoot, sequenceLength,
                    decodingQuantities.CSFSSamples, params.foldData,
                    params.usingCSFS)
        hmm = HMM(data, decodingQuantities, params, not params.noBatches, 1)

        p = makePairObs(data.individuals[0], 1, data.individuals[0], 2)
        d = hmm.decode(p)
        self.assertEqual(len(d), len(decodingQuantities.expectedTimes))
        for i in range(len(d)):
            self.assertEqual(len(d[i]), data.sites)
        del hmm


if __name__ == "__main__":
    unittest.main()
