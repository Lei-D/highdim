# NSCLC
load("~/Desktop/Lei/AIMER2/ISMB/real/NSCLC/suffPCR.RData")
print(suffPCR$best)
sum(suffPCR$betahat != 0) # 39 selected features
index = which(suffPCR$betahat != 0)
# read in data
load("~/Desktop/Lei/AIMER/NSCLC/rawdata/Chemores-miRNA-LVS-normalized.RData")
feature.name = row.names(mirdat)
feature.name[index]
# [1] "hsa-miR-376c"   "hsa-miR-320c"   "hsa-miR-299-3p" "hsa-miR-154"   
# [5] "hsa-miR-410"    "hsa-miR-1182"   "hsa-miR-136"    "hsa-miR-379"   
# [9] "hsa-miR-765"    "hsa-miR-610"    "hsa-miR-487a"   "hsa-miR-136*"  
# [13] "hsv1-miR-H1"    "hsa-miR-622"    "hsa-miR-659"    "hsa-miR-376b"  
# [17] "hsa-miR-154*"   "hsa-miR-483-5p" "hsa-miR-337-5p" "hsa-miR-376a*" 
# [21] "hsa-miR-409-5p" "hcmv-miR-US4"   "hsa-miR-376a"   "hsa-miR-1471"  
# [25] "hsa-miR-411"    "bkv-miR-B1-5p"  "hsa-miR-377"    "hsa-miR-601"   
# [29] "hsa-miR-299-5p" "hsa-miR-543"    "hsa-miR-381"    "hsa-miR-329"   
# [33] "hsa-miR-760"    "hsa-miR-409-3p" "hsa-miR-617"    "hsa-miR-758"   
# [37] "hsa-miR-1183"   "hsa-miR-671-5p" "hsa-miR-127-3p"


# AML
load("~/Desktop/Lei/AIMER2/ISMB/real/AML/suffPCR.RData")
print(suffPCR$best)
sum(suffPCR$betahat != 0) # 4 selected features
index = which(suffPCR$betahat != 0)
index
# [1] 4425 6012 6057 6059

# 105800 BPGM 2,3-bisphosphoglycerate mutase Hs.198365 AA678065  669 
# 309674 PIK4CB phosphatidylinositol 4-kinase, catalytic, beta polypeptide Hs.154846 AA282706  5298
# 109185  EST Hs.321434 H96229  
# 318982 LOC90379 hypothetical protein BC002926 Hs.298553 AI301071  90379 


# BreastCancer1
load("~/Desktop/Lei/AIMER2/ISMB/real/BreastCancer1/suffPCR.RData")
print(suffPCR$best)
sum(suffPCR$betahat != 0) # 28 selected features
index = which(suffPCR$betahat != 0)
index
# read in data
gene.name = read.csv("~/Desktop/Lei/AIMER/breastCancer/rawdata/vantVeer2002_gene_names.CSV", head=FALSE)
gene.name[index,1]
# [1] NM_003118      NM_003247      Contig46244_RC NM_004079      Contig55801_RC
# [6] NM_002775      M37033         NM_004369      NM_004385      NM_002985     
# [11] Contig43613_RC Contig43833_RC Contig42919_RC NM_005565      NM_016081     
# [16] Contig52398_RC NM_016187      NM_006889      Contig30260_RC Contig66347   
# [21] NM_000089      NM_000090      NM_000138      Contig25362_RC NM_000393     
# [26] NM_000560      Contig58512_RC NM_001387 


