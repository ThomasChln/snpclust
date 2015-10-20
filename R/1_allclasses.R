
# GenotypeDataSubset

setClass("GenotypeDataSubset", 
  contains = 'GenotypeData', 
  representation(snps_idx = 'integer', scans_idx = 'integer'),
  prototype(snps_idx = NULL, scans_idx = NULL)
)

#setOldClass("check.marker")

