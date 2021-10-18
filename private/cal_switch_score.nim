## calcualte switch score for snps positions that are high risk

let binSize = 1000
let movingStep = 200

# reduce the search space
# return seq of snp indexes

proc findHighRiskSnps(fullGeno:seq[int]): seq[int] = 