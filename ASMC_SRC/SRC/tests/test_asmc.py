import asmc


def test_asmc():
    hapsFileRoot = "FILES/EXAMPLE/exampleFile.n300.array"
    decodingQuantFile = "FILES/DECODING_QUANTITIES" + \
        "/30-100-2000.decodingQuantities.gz"
    decodingQuantities = asmc.DecodingQuantities(decodingQuantFile)
    params = asmc.DecodingParams(hapsFileRoot, decodingQuantFile)
    sequenceLength = asmc.Data.countHapLines(hapsFileRoot)
    data = asmc.Data(hapsFileRoot, sequenceLength,
                     decodingQuantities.CSFSSamples,
                     params.foldData, params.usingCSFS)
    hmm = asmc.HMM(data, decodingQuantities, params, not params.noBatches, False)
    hmm.decodeAll(params.jobs, params.jobInd)
