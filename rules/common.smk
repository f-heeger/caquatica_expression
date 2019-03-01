config: "config.json"

samples = {}
for line in open("samples.tsv"):
    sId, sName, sraId, condNum, culture, medium, stage = line.strip().split("\t")
    samples[sId] = {"ID": sId, "name": sName, "sraId": sraId, "condNum": condNum}
