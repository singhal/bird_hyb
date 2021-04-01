library(dplyr)

d = read.csv("~/Desktop/metadata/samples.csv",
             stringsAsFactors = F)

# read1, read2, adaptor1, adaptor2, barcode1, barcode2, lineage
dir = '/home/babs/Desktop/birds/Raw-Data-renamed/'

b = read.csv("~/Desktop/metadata/barcodes.csv")
b1 = b %>% select(SampleName, Barcode_Seq_i5, Barcode_Seq_i7) %>% unique()

# some individuals
# have two barcodes :(
d1 = left_join(d, b1)
# just drop one
d2 = d1[!duplicated(d1$SampleName), ]

# read names
d2$read1 = paste0(dir, d2$SampleName, '_R1.fastq.gz')
d2$read2 = paste0(dir, d2$SampleName, '_R2.fastq.gz')

# species
s = read.csv("~/Desktop/metadata/locality_info.csv", stringsAsFactors = F)
s$lineage = paste0(s$Genus, "_", s$Species) 
d3 = left_join(d2, s %>% select(sample_uid, lineage))

ad7s = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG'
ad5s = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

ad7d = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG'
ad5d = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT'

d3$adaptor1 = NA
d3$adaptor2 = NA
for (i in 1:nrow(d3)) {
  if (d3[i, "Barcode_Seq_i5"] == "NULL") {
    d3[i, "adaptor1"] = ad5s
    d3[i, "adaptor2"] = ad7s
  } else {
    d3[i, "adaptor1"] = ad5d
    d3[i, "adaptor2"] = ad7d
  }
}

names(d3)[ which(names(d3) == 'Barcode_Seq_i5') ] = "barcode1"
names(d3)[ which(names(d3) == 'Barcode_Seq_i7') ] = "barcode2"
d3[d3$barcode1 == "NULL", "barcode1"] = NA
write.csv(d3, "~/Desktop/bird_ILS.csv", row.names = F)
