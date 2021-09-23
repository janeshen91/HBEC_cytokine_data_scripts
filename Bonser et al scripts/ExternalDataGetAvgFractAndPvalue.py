import scanpy as sc

##import data downloaded from covid19cell atlas
viData=sc.read_h5ad("/Users/janeshen/Documents/covid19cellatlas/dataH5AD/vieira19_Bronchi_anonymised.processed.h5ad")
##subset out the 3 celltypes we're focusing on
vidata_sub = viData[viData.obs.CellType.isin(["Basal 1", "Basal 2","Club", "Goblet_1","Goblet 2", "Ciliated 1", "Ciliated 2"]),:]
#read in the files of cell type specific genes
my_file = open("/Users/janeshen/Documents/covid19cellatlas/ITGBspecificGenes.txt", "r")
content_list = my_file.readlines()
print(content_list)
geneList = [x.rstrip('\n') for x in content_list]
##Dotplot
ax = sc.pl.dotplot(vidata_sub, geneList, groupby='BroadCellType')

##print out the average fraction and the pvalues of the genes in geneList
import pandas as pd
gene_ids = vidata_sub.var.index.values
obs = vidata_sub[:,gene_ids].X.toarray()
obs = pd.DataFrame(obs,columns=gene_ids,index=vidata_sub.obs.index)

obs2=obs[geneList]
obs2['BroadCellType']= vidata_sub.obs.BroadCellType
average_obs = obs2.groupby('BroadCellType').mean()

sum3celltypes=average_obs.sum(axis=0)

cellNumbers=vidata_sub.obs.groupby(["BroadCellType"]).apply(len)###print out new!
cellNumbers.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/viBronchi_CTSB_CellNumbers.csv")
obs_bool = obs2.astype(bool)
obs_bool['BroadCellType'] = vidata_sub.obs.BroadCellType

fraction_obs = obs_bool.groupby('BroadCellType').sum()/obs_bool.groupby('BroadCellType').count()
average_obs.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/viBronchi_CTSB_BroadCelltype.csv")
fraction_obs.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/viBronchi_CTSB_BroadCelltype_fraction.csv")

sc.tl.rank_genes_groups(vidata_sub, 'BroadCellType', method='t-test',n_genes=10000)
ranksOfGenesValueVi=vidata_sub.uns.get('rank_genes_groups')
ViBronchiPvalAdj=ranksOfGenesValueVi['pvals_adj']

top100Vi=pd.DataFrame(vidata_sub.uns['rank_genes_groups']['names'])
top100Vi.to_csv("/Users/janeshen/Documents/covid19cellatlas/ViBronchi_CTSB_Log2fcPvalAdjustedTop100.txt")


top100ViPval=pd.DataFrame(ViBronchiPvalAdj)
top100ViPval.to_csv("/Users/janeshen/Documents/covid19cellatlas/ViBronchi_CTSB_Log2fcPvalAdjusted.txt")

#######repeat for vieira nasal data
viNasalData=sc.read_h5ad("/Users/janeshen/Documents/covid19cellatlas/dataH5AD/vieira19_Nasal_anonymised.processed.h5ad")

#ax = sc.pl.dotplot(viNasalData_sub, geneList, groupby='BroadCellType')
viNasalData_sub = viNasalData[viNasalData.obs.CellType.isin(["Basal 1", "Basal 2","Club", "Ciliated 1", "Ciliated 2"]),:]

import pandas as pd
gene_ids = viNasalData_sub.var.index.values
obs = viNasalData_sub[:,gene_ids].X.toarray()
obs = pd.DataFrame(obs,columns=gene_ids,index=viNasalData_sub.obs.index)

obs2=obs[geneList]
obs2['BroadCellType']= viNasalData_sub.obs.BroadCellType
average_obs = obs2.groupby('BroadCellType').mean()

cellNumbers=viNasalData_sub.obs.groupby(["BroadCellType"]).apply(len)###print out new!
cellNumbers.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/viNasal_CTSB_CellNumbers.csv")
obs_bool = obs2.astype(bool)
obs_bool['BroadCellType'] = viNasalData_sub.obs.BroadCellType

fraction_obs = obs_bool.groupby('BroadCellType').sum()/obs_bool.groupby('BroadCellType').count()
average_obs.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/viNasal_CTSB_BroadCelltype.csv")
fraction_obs.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/viNasal_CTSB_BroadCelltype_fraction.csv")

sc.tl.rank_genes_groups(viNasalData_sub, 'BroadCellType', method='t-test',n_genes=10000)
ranksOfGenesValueVi=viNasalData_sub.uns.get('rank_genes_groups')
ViBronchiPvalAdj=ranksOfGenesValueVi['pvals_adj']

top100Vi=pd.DataFrame(viNasalData_sub.uns['rank_genes_groups']['names'])
top100Vi.to_csv("/Users/janeshen/Documents/covid19cellatlas/ViNasal_CTSB_Log2fcPvalAdjustedTop100.txt")


top100ViPval=pd.DataFrame(ViBronchiPvalAdj)
top100ViPval.to_csv("/Users/janeshen/Documents/covid19cellatlas/ViNasal_CTSB_Log2fcPvalAdjusted.txt")


##repeat for lukassen data
lukassenData=sc.read_h5ad("/Users/janeshen/Documents/covid19cellatlas/dataH5AD/lukassen20_airway_orig.processed.h5ad")

lukassenData_sub = lukassenData[lukassenData.obs.CellType.isin(["Basal1", "Basal2","Basal3","Basal_Mitotic","Club", "Goblet","Goblet", "Ciliated1", "Ciliated2","Secretory1","Secretory2","Secretory3"]),:]
ax = sc.pl.dotplot(lukassenData_sub, geneList, groupby='CellType')



gene_ids = lukassenData_sub.var.index.values
obs = lukassenData_sub[:,gene_ids].X.toarray()
obs = pd.DataFrame(obs,columns=gene_ids,index=lukassenData_sub.obs.index)

obs2=obs[geneList]
obs2['CellType']= lukassenData_sub.obs.CellType
average_obs = obs2.groupby('CellType').mean()

average_obs.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/Lukassen_CTSB_Avg.csv")


cellNumbers=lukassenData_sub.obs.groupby(["CellType"]).apply(len)###print out new!
cellNumbers.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/Lukassen_CTSB_CellNumbers.csv")

obs_bool = obs2.astype(bool)
obs_bool['CellType'] = lukassenData_sub.obs.CellType
fraction_obs = obs_bool.groupby('CellType').sum()/obs_bool.groupby('CellType').count()
fraction_obs.T.to_csv("/Users/janeshen/Documents/covid19cellatlas/Lukassen_CTSB_fraction.csv")

sc.tl.rank_genes_groups(lukassenData_sub, 'CellType', method='t-test',n_genes=10000)

ranksOfGenesValueVi=lukassenData_sub.uns.get('rank_genes_groups')
ViBronchiPvalAdj=ranksOfGenesValueVi['pvals_adj']

top100Vi=pd.DataFrame(lukassenData_sub.uns['rank_genes_groups']['names'])
top100Vi.to_csv("/Users/janeshen/Documents/covid19cellatlas/Lukassen_CTSB_PvalAdjustedTop100.txt")


top100ViPval=pd.DataFrame(ViBronchiPvalAdj)
top100ViPval.to_csv("/Users/janeshen/Documents/covid19cellatlas/Lukassen_CTSB_PvalAdjusted.txt")


