############################################################################
#===========================================================================
#
#--------- Transcriptome data to be added to the somalogic paper ---------
#
#===========================================================================
#
#------------ Biomarker subjects 83-83 --------------------------
#
setwd("/Users/seidmuhie/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/10282020_Somalogic_WoundHealing_toScienceWriter/04052022_Manuscript_for_MolecularCell_Submission")
hptsd8383ge = read.csv("OriginalMale-mRNA-PBMC-USACEHR-Agilent-20171031.csv", row.names = 1)
dim(hptsd8383ge)
names(hptsd8383ge[1:10])



#PCA
pca.ma<-prcomp(as.matrix(hptsd8383ge))
pc1<-summary(pca.ma)[6][[1]][2,1]
pc2<-summary(pca.ma)[6][[1]][2,2]
plot(pca.ma$rotation[,1:2],col=chip.col,pch=ifelse(is.ptsd,"y","n"),
     ylab=paste("PC2 (",pc2,")",sep=""),xlab=paste("PC1 (",pc1,")",sep=""),
     main="M PCA after norm",sub="symbol=PTSD y/n, color=chip")


plot(pca.ma$rotation[,1:2],pch=sub(".*_", "",names(hptsd8383ge)), #, value=TRUE), grep("_P", names(hptsd8383ge), value=TRUE)),#col=chip.col,
     ylab=paste("PC2 (",pc2,")",sep=""),xlab=paste("PC1 (",pc1,")",sep=""),
     main="PCA for norm",sub="symbol=PTSD P/N")

xtring = "hello xxx other stuff"

sub(" xxx.*", "", xtring) # subset before xxx
sub(".*xxx ", "", xtring) # subset after xxx

sub(".*_", "",names(hptsd8383ge))


names(hptsd8383ge_ptsd) = sub(".*_", "",names(hptsd8383ge))
View(hptsd8383ge_ptsd)
hptsd8383ge_ptsd_matrix = as.matrix(hptsd8383ge_ptsd)

#fit linear model
design <- model.matrix(~0+factor(names(hptsd8383ge_ptsd))) #
colnames(design) <- c("P","N")
fit <- lmFit(hptsd8383ge_ptsd_matrix, design)
cont.matrix <- makeContrasts(caseVScontrol=P-N, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- classifyTestsF(fit2)
topTable.all<-topTable(fit2,sort="P",n=Inf,adjust="BH")

write.csv(topTable.all, "hPTSD_transcriptome_DEGs.csv", row.names = T)

View(topTable.all[topTable.all$P.Value<0.05,])


#############################################################################################
#
#-------------- GSE64814-GPL6244_Woelk_96samples_molecPsychiatry -------------------------
#
#==========================================================================================

setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/10282020_Somalogic_WoundHealing_toScienceWriter/04052022_Manuscript_for_MolecularCell_Submission/Marin-resilience-study")

marindata_96samples = read.csv("GSE64814-GPL6244_Woelk_96samples_molecPsychiatry.csv", row.names = 1)
dim(marindata_96samples)

marindata_96samples_ptsd = marindata_96samples
names(marindata_96samples_ptsd) = sub("_.*", "",names(marindata_96samples))
View(marindata_96samples_ptsd)
marindata_96samples_ptsd_matrix = as.matrix(marindata_96samples_ptsd)

#fit linear model
design <- model.matrix(~0+factor(names(marindata_96samples_ptsd))) #
colnames(design) <- c("Case","Control")
fit <- lmFit(marindata_96samples_ptsd_matrix, design)
cont.matrix <- makeContrasts(caseVScontrol=Case-Control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- classifyTestsF(fit2)
topTable.all<-topTable(fit2,sort="P",n=Inf,adjust="BH")

write.csv(topTable.all, "GSE64814-GPL6244_Woelk_96samples_DEGs.csv", row.names = T)

View(topTable.all[topTable.all$P.Val<0.05,])

annotation_mrdata = read.csv("GPL6244-24073_GSE64814_annotationFile_simplified_GeneSymbol.csv")
names(annotation_mrdata)

toptable.all_annotated = merge(topTable.all, annotation_mrdata, by.x = 0, by.y = "ID")

View(toptable.all_annotated[toptable.all_annotated$P.Value<0.01,])

write.csv(toptable.all_annotated, "GSE64814-GPL6244_Woelk_96samples_DEGs.csv", row.names = T)

#====================================================================================================
#
#----------------- post-case vs post-control ----------------------------
#


marindata_96samples_ptsd_post = marindata_96samples
names(marindata_96samples_ptsd_post) = sub("_._*._.*", "",names(marindata_96samples))
View(marindata_96samples_ptsd_post)
marindata_96samples_ptsd_post_matrix = as.matrix(marindata_96samples_ptsd_post)

#fit linear model
design <- model.matrix(~0+factor(names(marindata_96samples_ptsd))) #
colnames(design) <- c("Case","Control")
fit <- lmFit(marindata_96samples_ptsd_matrix, design)
cont.matrix <- makeContrasts(caseVScontrol=Case-Control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- classifyTestsF(fit2)
topTable.all<-topTable(fit2,sort="P",n=Inf,adjust="BH")

write.csv(topTable.all, "GSE64814-GPL6244_Woelk_96samples_DEGs.csv", row.names = T)

View(topTable.all[topTable.all$P.Val<0.05,])

annotation_mrdata = read.csv("GPL6244-24073_GSE64814_annotationFile_simplified_GeneSymbol.csv")
names(annotation_mrdata)

toptable.all_annotated = merge(topTable.all, annotation_mrdata, by.x = 0, by.y = "ID")

View(toptable.all_annotated[toptable.all_annotated$P.Value<0.01,])

write.csv(toptable.all_annotated, "GSE64814-GPL6244_Woelk_96samples_DEGs.csv", row.names = T)

#==============================================================================
#
#------------ common between somalogic and postmorteme datasets --------------
#
setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/10282020_Somalogic_WoundHealing_toScienceWriter/04052022_Manuscript_for_MolecularCell_Submission")
commonids = read.csv("commonIDs_SL_PM.csv", row.names = 1)
dlpfcdata = read.csv("postMortem brain data 41593_2020_748_MOESM3_ESM_diPFC_both.csv", row.names = 1)

dlpfc_commonids = merge(commonids, dlpfcdata, by = 0)
colnames(dlpfc_commonids)[1] <- "GeneName"
dlpfc_commonidsv2 = dlpfc_commonids[,-2]
names(dlpfc_commonidsv2)

sl_pathways = read.csv("01_0_hPTSD_SomaScan_Proteins_metaTopTables_p_0.1_pathways_BPs_q_0.05_wide_good_v2.csv")
dim(sl_pathways)

somalogic_dlpfc_pathwaydata = merge(sl_pathways ,dlpfc_commonidsv2, by = "GeneName", all = T)
dim(somalogic_dlpfc_pathwaydata)



commonids_mr = read.csv("commonIDs_SL_MR.csv")
dim(commonids_mr)
names(commonids_mr)

marindata = read.csv("GSE64814-GPL6244_Woelk_96samples_DEGs.csv")
names(marindata)
marindata_commonid = merge(commonids_mr, marindata, by = "GeneName")
dim(marindata_commonid)


length(intersect(somalogic_dlpfc_pathwaydata$GeneName, marindata_commonid$GeneName))

somalogic_dlpfc_marin_data_pathways = merge(somalogic_dlpfc_pathwaydata, marindata_commonid, by ="GeneName", all = T)

dim(somalogic_dlpfc_pathwaydata)
dim(somalogic_dlpfc_marin_data_pathways)
write.csv(somalogic_dlpfc_marin_data_pathways, "hPTSD_somalogic_post-mortem_marin_data_pathways.csv", row.names = F)

#============ somalogic and sbc transcriptome intersection ==================

commonid_sbctranscripts = read.csv("commonIDs_SL_SBCtranscriptome.csv")
dim(commonid_sbctranscripts)
sbc_trasncriptsdata = read.csv("hPTSD_transcriptome_DEGs_v3.csv")
dim(sbc_trasncriptsdata)

sbc_commonid_transcriptdata = merge(commonid_sbctranscripts, sbc_trasncriptsdata, by = "GeneName")
dim(sbc_commonid_transcriptdata)

dim(somalogic_dlpfc_marin_data_pathways)

sl_dlpfc_marin_sbctrans_data_pathways = merge(somalogic_dlpfc_marin_data_pathways, sbc_commonid_transcriptdata, by = "GeneName", all = T)

dim(sl_dlpfc_marin_sbctrans_data_pathways)

write.csv(sl_dlpfc_marin_sbctrans_data_pathways, "hPTSD_somalogic_post-mortem_marin_sbc-transcripts_data_pathways_good.csv", row.names = F)






#####################################################################################

#===========================================================================================================
#
#   ----------      proteins common and unique to active duty and veteran cohorts   -----------------
#


library(dplyr)
library(ComplexHeatmap)
library(seriation)
library(circlize)

rm(list = ls())
options(stringsAsFactors = F)

setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/10282020_Somalogic_WoundHealing_toScienceWriter/04052022_Manuscript_for_MolecularCell_Submission")
# up in active duty

protein_transcripts_logfc_pvalue = read.csv("hPTSD_somalogic_post-mortem_data_pathways_good_noGrady_noFemales.csv")
names(protein_transcripts_logfc_pvalue)
dim(protein_transcripts_logfc_pvalue)


protein_transcripts_logfc_pvalue_inflammation = protein_transcripts_logfc_pvalue[protein_transcripts_logfc_pvalue$Pathway=="Inflammatory response",]
dim(protein_transcripts_logfc_pvalue_inflammation)
inflammatory_response = protein_transcripts_logfc_pvalue_inflammation[,-2]

row.names(inflammatory_response) = inflammatory_response$GeneName
inflammatory_response = inflammatory_response[,-1]
View(inflammatory_response)

# colnames(protein_logfc_pvalue) = c("Training_M_VT", "Female_VT", "Female_AD", "Longitudinal_M_AD", "Testing_M_AD", 
#                                    "Subthreshold_M_AD", "civilian_M", "Validating_M_VT", "Training_M_VT_pval", "Female_VT_pval", "Female_AD_pval", "Longitudinal_M_AD_pval", "Testing_M_AD_pval", 
#                                    "Subthreshold_M_AD_pval", "civilian_M_pval", "Validating_M_VT_val")


#================================================================================================
#
#--------- both up and down across all cohorts (for inflammatory response proteins) -----------------------------------------------
#
# ------  down across cohorts ------------
#
names(inflammatory_response)

inflammatory_response_down = inflammatory_response[!(inflammatory_response[,1]>0|inflammatory_response[,2]>0|inflammatory_response[,3]>0|inflammatory_response[,4]>0|inflammatory_response[,5]>0|inflammatory_response[,16]>0|inflammatory_response[,18]>0|inflammatory_response[,20]>0|inflammatory_response[,22]>0|inflammatory_response[,24]>0|inflammatory_response[,27]>0),]
View(inflammatory_response_down)

#---  up across cohorts

protein_logfc_pvalue_up = protein_logfc_pvalue[!(protein_logfc_pvalue[,1]<0|protein_logfc_pvalue[,2]<0|protein_logfc_pvalue[,3]<0|protein_logfc_pvalue[,4]<0|protein_logfc_pvalue[,5]<0|protein_logfc_pvalue[,6]<0|protein_logfc_pvalue[,7]<0|protein_logfc_pvalue[,8]<0),]
View(protein_logfc_pvalue_up)


#--------- both up and down-regulated across cohorts ---

protein_logfc_pvalue_updown = rbind(protein_logfc_pvalue_up, protein_logfc_pvalue_down)

write.csv(protein_logfc_pvalue_updown, "proteins consistent across all cohorts.csv", row.names = T)

View(protein_logfc_pvalue_updown)

#--------- inflammatory response heatmap ----------------

inflammatory_response_matrix = as.matrix(inflammatory_response)

o1 = seriate(dist(inflammatory_response_matrix[,c(1:5,16,18,20,22,24)]), method = "GW")
o2 = seriate(dist(t(inflammatory_response_matrix[,c(1:5,16,18,20,22,24)])), method = "GW")


small_mat = inflammatory_response_matrix[,c(1:5,16,18,20,22,24)]

small_mat_pv = inflammatory_response_matrix[,c(6:10,17,19,21,23,25)]

#         working code for heatmap 

library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp1 = Heatmap(small_mat, name = "log2FC", col = col_fun,
                show_heatmap_legend = F,
                #row_km = 2, column_km = 2,
                column_names_rot = 45,
                row_km = 2,
                rect_gp = gpar(col = "white", lwd = 2),
                
                layer_fun = function(j, i, x, y, w, h, fill) {
                  # restore_matrix() is explained after this chunk of code
                  ind_mat = restore_matrix(j, i, x, y)
                  for(ir in seq_len(nrow(ind_mat))) {
                    # start from the second column
                    for(ic in seq_len(ncol(ind_mat))) {
                      ind = ind_mat[ir, ic] # previous column
                      v = small_mat_pv[i[ind], j[ind]]
                      
                      grid.points(x[ind], y[ind], 
                                  pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(3, "mm"))
                      
                    }
                  }
                }
                
)

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2 (fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4, "cm"))

lgd2 = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "p-value", type = "points", title_position =  "topleft", 
              grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
              pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")
pd_ldgs = packLegend(lgd1, lgd2)

draw(htmp1, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")




#--------- inflammatory response significant in SBC protein data ----------------
inflammatory_response_sbc_significant = inflammatory_response[inflammatory_response$adj.P.Val.TT.Discovery_validation_topTable_cohort_race_pc1.2_HybConNorScl.csv<0.1,]
dim(inflammatory_response_sbc_significant)

inflammatory_response_sbc_significant_matrix = as.matrix(inflammatory_response_sbc_significant)

o1 = seriate(dist(inflammatory_response_sbc_significant_matrix[,c(1:5,16,18,20,22,24,27)]), method = "GW")
o2 = seriate(dist(t(inflammatory_response_sbc_significant_matrix[,c(1:5,16,18,20,22,24,27)])), method = "GW")


small_mat = inflammatory_response_sbc_significant_matrix[,c(1:5,16,18,20,22,24,27)]

small_mat_pv = inflammatory_response_sbc_significant_matrix[,c(6:10,17,19,21,23,25,28)]

#         working code for heatmap 

library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp1 = Heatmap(small_mat, name = "log2FC", col = col_fun,
                show_heatmap_legend = F,
                #row_km = 2, column_km = 2,
                column_names_rot = 45,
                row_km = 2,
                rect_gp = gpar(col = "white", lwd = 2),
                
                layer_fun = function(j, i, x, y, w, h, fill) {
                  # restore_matrix() is explained after this chunk of code
                  ind_mat = restore_matrix(j, i, x, y)
                  for(ir in seq_len(nrow(ind_mat))) {
                    # start from the second column
                    for(ic in seq_len(ncol(ind_mat))) {
                      ind = ind_mat[ir, ic] # previous column
                      v = small_mat_pv[i[ind], j[ind]]
                      
                      grid.points(x[ind], y[ind], 
                                  pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(3, "mm"))
                      
                    }
                  }
                }
                
)

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2 (fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4, "cm"))

lgd2 = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "p-value", type = "points", title_position =  "topleft", 
              grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
              pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")
pd_ldgs = packLegend(lgd1, lgd2)

draw(htmp1, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")

#==============================================================
#
#------ heatmap for consistent inflammatory proteins ------
#

inflamm_proteins_consistent = c("IL1RAP", "SERPINC1", "OLR1",  
"NOTCH1", "CFP", "CHI3L1", "C4B", "S100A9",
"ANXA1", "STAT3", "CEBPB", "C2", "PRKCD", "MAPK13", "SELE",
"BDNF", "S100A12", "IL1B", "SERPINE1", "CD59", "PRKCZ", "CHST2",
"GIG24", "CCL13", "INS", "TNF", "TEK", "IL6ST", "AIMP1", "F2", 
"NCR3", "SERPINA3", "JAK2", "KLK3", "NPPB", "CRP", "TF", "A2M", "CCL5")
#( "CFI", "CAMK1D", "AZU1", "HIF1A") # these were not consistent throughout

inflammatory_prot_consistent = inflammatory_response[inflamm_proteins_consistent,]
View(inflammatory_prot_consistent)
write.csv(inflammatory_prot_consistent, "heatmap inflammatory-response.csv")
inflammatory_prot_consistent_matrix = as.matrix(inflammatory_prot_consistent)

o1 = seriate(dist(inflammatory_prot_consistent_matrix[,c(1:5,16,18,20,22,24)]), method = "GW")
o2 = seriate(dist(t(inflammatory_prot_consistent_matrix[,c(1:5,16,18,20,22,24)])), method = "GW")


small_mat = inflammatory_prot_consistent_matrix[,c(1:5,16,18,20,22,24)]

small_mat_pv = inflammatory_prot_consistent_matrix[,c(6:10,17,19,21,23,25)]

#         working code for heatmap 

library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp1 = Heatmap(small_mat, name = "log2FC", col = col_fun,
                show_heatmap_legend = F,
                #row_km = 2, column_km = 2,
                column_names_rot = 45,
                row_km = 2,
                rect_gp = gpar(col = "white", lwd = 2),
                
                layer_fun = function(j, i, x, y, w, h, fill) {
                  # restore_matrix() is explained after this chunk of code
                  ind_mat = restore_matrix(j, i, x, y)
                  for(ir in seq_len(nrow(ind_mat))) {
                    # start from the second column
                    for(ic in seq_len(ncol(ind_mat))) {
                      ind = ind_mat[ir, ic] # previous column
                      v = small_mat_pv[i[ind], j[ind]]
                      
                      grid.points(x[ind], y[ind], 
                                  pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(3, "mm"))
                      
                    }
                  }
                }
                
)

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2 (fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4, "cm"))

lgd2 = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "p-value", type = "points", title_position =  "topleft", 
              grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
              pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")
pd_ldgs = packLegend(lgd1, lgd2)

draw(htmp1, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")






#==============================================================
#
#------ heatmap for oxidative stress proteins ------
#

protein_transcripts_logfc_pvalue_oxidativestress = protein_transcripts_logfc_pvalue[protein_transcripts_logfc_pvalue$Pathway=="response to oxidative stress",]
dim(protein_transcripts_logfc_pvalue_oxidativestress)
oxidative_stress = protein_transcripts_logfc_pvalue_oxidativestress[,-2]

row.names(oxidative_stress) = oxidative_stress$GeneName
oxidative_stress = oxidative_stress[,-1]
View(oxidative_stress)

rownames(oxidative_stress)

oxidativestress_proteins_consistent = c("ADAM9", "ADRBK1", "BCL2", "CAT", "CCL5", "CLU",      
                                        "HIF1A", "HMOX2", "JAK2", "MAP2K1", "MPO", "OLR1",
                                        "PTGS2", "SERPINE1", "SNCA", "SOD1", "SOD2", "TF")

oxidativestress_proteins_consistentdata = oxidative_stress[oxidativestress_proteins_consistent,]
View(oxidativestress_proteins_consistentdata)
write.csv(oxidativestress_proteins_consistentdata, "heatmap reponse to oxidative stress.csv")

oxidativestress_proteins_consistentdata_matrix = as.matrix(oxidativestress_proteins_consistentdata)

o1 = seriate(dist(oxidativestress_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)]), method = "GW")
o2 = seriate(dist(t(oxidativestress_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)])), method = "GW")


small_mat = oxidativestress_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)]

small_mat_pv = oxidativestress_proteins_consistentdata_matrix[,c(6:10,17,19,21,23)]

#         working code for heatmap 

library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp1 = Heatmap(small_mat, name = "log2FC", col = col_fun,
                show_heatmap_legend = F,
                #row_km = 2, column_km = 2,
                column_names_rot = 45,
                row_km = 2,
                rect_gp = gpar(col = "white", lwd = 2),
                
                layer_fun = function(j, i, x, y, w, h, fill) {
                  # restore_matrix() is explained after this chunk of code
                  ind_mat = restore_matrix(j, i, x, y)
                  for(ir in seq_len(nrow(ind_mat))) {
                    # start from the second column
                    for(ic in seq_len(ncol(ind_mat))) {
                      ind = ind_mat[ir, ic] # previous column
                      v = small_mat_pv[i[ind], j[ind]]
                      
                      grid.points(x[ind], y[ind], 
                                  pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(3, "mm"))
                      
                    }
                  }
                }
                
)

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2 (fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4, "cm"))

lgd2 = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "p-value", type = "points", title_position =  "topleft", 
              grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
              pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")
pd_ldgs = packLegend(lgd1, lgd2)

draw(htmp1, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")



#======================================================================
#
#------ heatmap for vasculature development (angiogenesis) ------
#

protein_transcripts_logfc_pvalue_angiogenesis = protein_transcripts_logfc_pvalue[protein_transcripts_logfc_pvalue$Pathway=="vasculature development",]
dim(protein_transcripts_logfc_pvalue_angiogenesis)
angiogenesis = protein_transcripts_logfc_pvalue_angiogenesis[,-2]

row.names(angiogenesis) = angiogenesis$GeneName
angiogenesis = angiogenesis[,-1]
View(angiogenesis)

rownames(angiogenesis)

angiogenesis_proteins_consistent = c("AIMP1", "BGN", "COL18A1", "COL8A1",      
                                     "CDH5", "CX3CL1", "EFNB2", "ENG", "EPHA2", 
                                     "ESM1", "FLT4", "FN1", "GHRL", "GPI",   
                                     "KDR", "KLK3", "LRP1", "MAP2K1", "MED1",
                                     "MMP2", "NOTCH1", "NPPB", "NRCAM",  
                                     "NRXN3", "NTRK2", "PLAT", "PTGS2", "SERPINE1",   
                                     "SPINT1", "TEK", "TGFBI", "TGFBR3", "TIE1")

angiogenesis_proteins_consistentdata = angiogenesis[angiogenesis_proteins_consistent,]
dim(angiogenesis_proteins_consistentdata)
write.csv(angiogenesis_proteins_consistentdata, "heatmap angiogenesis.csv")


angiogenesis_proteins_consistentdata_matrix = as.matrix(angiogenesis_proteins_consistentdata)

o1 = seriate(dist(angiogenesis_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)]), method = "GW")
o2 = seriate(dist(t(angiogenesis_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)])), method = "GW")


small_mat = angiogenesis_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)]

small_mat_pv = angiogenesis_proteins_consistentdata_matrix[,c(6:10,17,19,21,23)]

#         working code for heatmap 

library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp1 = Heatmap(small_mat, name = "log2FC", col = col_fun,
                show_heatmap_legend = F,
                #row_km = 2, column_km = 2,
                column_names_rot = 45,
                row_km = 2,
                rect_gp = gpar(col = "white", lwd = 2),
                
                layer_fun = function(j, i, x, y, w, h, fill) {
                  # restore_matrix() is explained after this chunk of code
                  ind_mat = restore_matrix(j, i, x, y)
                  for(ir in seq_len(nrow(ind_mat))) {
                    # start from the second column
                    for(ic in seq_len(ncol(ind_mat))) {
                      ind = ind_mat[ir, ic] # previous column
                      v = small_mat_pv[i[ind], j[ind]]
                      
                      grid.points(x[ind], y[ind], 
                                  pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(3, "mm"))
                      
                    }
                  }
                }
                
)

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2 (fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4, "cm"))

lgd2 = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "p-value", type = "points", title_position =  "topleft", 
              grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
              pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")
pd_ldgs = packLegend(lgd1, lgd2)

draw(htmp1, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")


#======================================================================
#
#------ heatmap for heart development  ------
#

protein_transcripts_logfc_pvalue_heartdevelopment = protein_transcripts_logfc_pvalue[protein_transcripts_logfc_pvalue$Pathway=="heart development",]
dim(protein_transcripts_logfc_pvalue_heartdevelopment)
heartdevelopment = protein_transcripts_logfc_pvalue_heartdevelopment[,-2]

row.names(heartdevelopment) = heartdevelopment$GeneName
heartdevelopment = heartdevelopment[,-1]
dim(heartdevelopment)

rownames(heartdevelopment)

heartdevelopment_proteins_consistent = c("ACAN", "ADRBK1", "BMP10", "BMPR1A", 
                        "ENG", "FGF12", "FLRT3", "INSR","MAP2K1", 
                        "MAPK1", "MB", "MDM2", "MED1", "NOTCH1", "NRG1", "NRP1",   
                        "NTRK3", "PDGFB", "SHC1", "SHH", "SOD2", "SPARC", "TEK",    
                        "TGFBR2", "TGFBR3", "TNNI3", "VCAM1")

heartdevelopment_proteins_consistentdata = heartdevelopment[heartdevelopment_proteins_consistent,]
dim(heartdevelopment_proteins_consistentdata)
write.csv(heartdevelopment_proteins_consistentdata, "heatmap heart development.csv")


heartdevelopment_proteins_consistentdata_matrix = as.matrix(heartdevelopment_proteins_consistentdata)

o1 = seriate(dist(heartdevelopment_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)]), method = "GW")
o2 = seriate(dist(t(heartdevelopment_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)])), method = "GW")


small_mat = heartdevelopment_proteins_consistentdata_matrix[,c(1:5,16,18,20,22)]

small_mat_pv = heartdevelopment_proteins_consistentdata_matrix[,c(6:10,17,19,21,23)]

#         working code for heatmap 

library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp1 = Heatmap(small_mat, name = "log2FC", col = col_fun,
                show_heatmap_legend = F,
                #row_km = 2, column_km = 2,
                column_names_rot = 45,
                row_km = 2,
                rect_gp = gpar(col = "white", lwd = 2),
                
                layer_fun = function(j, i, x, y, w, h, fill) {
                  # restore_matrix() is explained after this chunk of code
                  ind_mat = restore_matrix(j, i, x, y)
                  for(ir in seq_len(nrow(ind_mat))) {
                    # start from the second column
                    for(ic in seq_len(ncol(ind_mat))) {
                      ind = ind_mat[ir, ic] # previous column
                      v = small_mat_pv[i[ind], j[ind]]
                      
                      grid.points(x[ind], y[ind], 
                                  pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(3, "mm"))
                      
                    }
                  }
                }
                
)

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2 (fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4, "cm"))

lgd2 = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "p-value", type = "points", title_position =  "topleft", 
              grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
              pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")
pd_ldgs = packLegend(lgd1, lgd2)

draw(htmp1, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")




##############################################################################
#
#----- intersecting mir-dmr-metabolites with unmatched list of proteins -----
#
setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/10282020_Somalogic_WoundHealing_toScienceWriter/04052022_Manuscript_for_MolecularCell_Submission/05102022_SummaryHeatmaps")

proteinlist = read.csv("heatmap_unmatched_list.csv", row.names = 1)
dim(proteinlist)
dmrmirmetabolite = read.csv("DMR-miR-protein-metabolite-protein_GOOD_V2.csv")
dim(dmrmirmetabolite)
proteindmrmirmetabolite = merge(proteinlist, dmrmirmetabolite, by.x = 0, by.y = "GeneName")
View(proteindmrmirmetabolite)
write.csv(proteindmrmirmetabolite, "umatched prteins dmr mir metabolites.csv", row.names = F)

################################################################################
#
#------ intersect gene names from pathways with the different modalities -----
#
pathwaygenenames = read.csv("pathways_list of GeneNames.csv")
dim(pathwaygenenames)

gwaspqtl = read.csv("table1_GWAS_and_pQTL_simplified.csv")
dim(gwaspqtl)

gwasqptl_pathway = merge(pathwaygenenames, gwaspqtl, by = "GeneName", all.x = T)
dim(gwasqptl_pathway)

write.csv(gwasqptl_pathway, "GWAS-pQTL pathways.csv")

ptsd = read.csv("ptsd.csv")
mdd = read.csv("mdd.csv")
ukbb = read.csv("ukbb.csv")

ptsd_mdd = merge(ptsd, mdd, by="GeneName", all = T)
write.csv(ptsd_mdd, "ptsd-mdd.csv", row.names = F)
ptsd_mdd_ukbb = merge(ptsd_mdd, ukbb, by = "GeneName", all = T)
write.csv(ptsd_mdd_ukbb, "ptsd mdd ukbb.csv", row.names = F)

###########################################################################
#
#---------- multimodal itersection with pathway protein list --------------
#

pathwaygenenames = read.csv("pathways_list of GeneNames.csv")
View(pathwaygenenames)

#------- intersecting dmr and pathways along with protein list ------

dmr = read.csv("DMR inflammation angiogenesis DMR-protein good_v3.csv")
protein_dmr = merge(pathwaygenenames, dmr, by = "GeneName", all.x = T)
dim(protein_dmr)
write.csv(protein_dmr, "DMR inflammation angiogenesis DMR-protein_v4_good.csv", row.names = F)

# ---- intersecting miRs and pathways along with protein list ---------

mir = read.csv("miR inflammation angiogenesis mir-protein good_v2.csv")
protein_mir = merge(pathwaygenenames, mir, by = "GeneName", all.x = T)
dim(protein_mir)
write.csv(protein_mir, "miR inflammation angiogenesis mir-protein_v4_good.csv", row.names = F)

# ---- intersecting metabolites and pathways along with protein list ---------

metabolites = read.csv("metabolites inflammation angiogensis protein-metabolite good_v4.csv")
protein_metabolites = merge(pathwaygenenames, metabolites, by = "GeneName", all.x = T)
dim(protein_metabolites)
write.csv(protein_metabolites, "metabolites inflammation angiogenesis proteins-metabolites_v4_good.csv", row.names = F)

##########################################################################################
#
# -- adding back the marine and SBC-transcript data to the protein data of 87 proteins
#
#======================================================================================

# read the marine and SBC transcript fold changes and p-values
marinsbctranscript = read.csv("hPTSD_marin_sbc-transcripts_data_good.csv")
dim(marinsbctranscript)

# read the protein data across the male military subgroups (of SBC and FCC cohorts)
combinedpathways = read.csv("heatmap combined oxidative stress inflammation angiogenesis heart development.csv", row.names = 1)
dim(combinedpathways)

# - merge the the marine and sbc transcript data with the protein data and pathways
combinedpathwayswithmarinsbctranscripts = merge(combinedpathways, marinsbctranscript, by.x = 0, by.y = "GeneName")
dim(combinedpathwayswithmarinsbctranscripts)

write.csv(combinedpathwayswithmarinsbctranscripts, "heatmap combined pathways proteins postmortem sbc marin.csv")
###############################################################################
#
#--------------  annotation of the agilent platform ------------
#

agilentannotation = read.csv("Agilent_GE_annotation_039494_D_AA_20210927.csv")
dim(agilentannotation)
# [1] 50683    14

ge_toptable = read.csv("hPTSD_transcriptome_DEGs.csv")
dim(ge_toptable)
# [1] 50599    31

get_toptable_annotated = merge(agilentannotation, ge_toptable, by = "ProbeID")

#write.csv(get_toptable_annotated, "hPTSD_transcriptome_DEGs_annotated.csv")

get_toptable_annotated2 = read.csv("hPTSD_transcriptome_DEGs_annotated_v2.csv")

#---- merging protein toptables with protein annotation with geneIDs---

sbc_training_proteins_id = read.csv("SBC_Training_SomalogicFCq-values_geneID.csv")

length(intersect(sbc_training_proteins_id$EntrezGeneID,get_toptable_annotated2$EntrezGeneID))



sbc_training_proteins_transcriptome_common_annotated = merge(get_toptable_annotated2, 
                                                             sbc_training_proteins_id, by = "EntrezGeneID")

dim(sbc_training_proteins_transcriptome_common_annotated)


write.csv(sbc_training_proteins_transcriptome_common_annotated, "sbc_training_proteins_transcriptome_common_annotated.csv")

##################################################################################
##################################################################################
#================================================================================
#
#----------------- heatmap for all combined pathway: 
# inflammation, oxidative stress, angiogenesis, heart development ------------
#
#==============================================================================

########################################################################################################
#======================================================================================================
#
# ----------------- Annotations for GWAS miR pQTL heatmaps --------------------------------
#
#=====================================================================================================

library(dplyr)

gwaspqtlmirpathways = read.csv("GWAS-pQTL-miR pathways_v4_good.csv", row.names = 1)
#mirprotein = mirprotein[-c(35,37,69),]
# subsetting the data based on the common proteins between the two sets of heatmaps
#mirprotein =  mirprotein[!(row.names(mirprotein) %in% common16proteins_toberemoved),]

View(gwaspqtlmirpathways)
gwaspqtlmiroxoinflammation = gwaspqtlmirpathways[gwaspqtlmirpathways$pathway=='inflammatory response'|gwaspqtlmirpathways$pathway=='oxidative stress',]

pvaluemir_oxinfl = gwaspqtlmiroxoinflammation$P.Value_mir
pvalueptsd_oxinfl = gwaspqtlmiroxoinflammation$P_GWAS_pgc_ptsd
pvaluemdd_oxinfl = gwaspqtlmiroxoinflammation$P_GWAS_pgc_mdd
pvalueukbd_oxinfl = gwaspqtlmiroxoinflammation$P_GWAS_ukbb_dep
pvaluenull_oxinfl = gwaspqtlmiroxoinflammation$nopvalue
pvaluepqtlsbc = gwaspqtlmiroxoinflammation$P_pQQTL_SBC
pvaluepqtlfcc = gwaspqtlmiroxoinflammation$P_pQQTL_FCC


length(pvaluepqtlfcc)
pvalue = mirprotein$P.Value_mir 

fold_col_fun = colorRamp2(c(-2, 0, 2), c("cyan", "white", "orange")) 
ha_miR_oxoinflm = rowAnnotation(.= anno_simple(gwaspqtlmiroxoinflammation$logFC_mir,  col = fold_col_fun, na_col = "white", pch = ifelse((pvaluemir_oxinfl<0.05 & !is.na(pvaluemir_oxinfl)),"*","")),
                       mir = anno_text(gwaspqtlmiroxoinflammation$miR_ID,gp = gpar(fontsize = 10, fontface = "bold"), just = "left", location = unit(0.05, "npc")))

fold_col_fun = colorRamp2(c(-2, 0, 2), c("cyan", "white", "orange")) 
ha_gwas_ptsd_oxinfl = rowAnnotation(ptsd = anno_text(gwaspqtloxoinflammation$SNP_ptsd,gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                           .= anno_simple(gwaspqtloxoinflammation$GWAS_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvalueptsd_oxinfl<0.005,"**","*")))
                           #mdd = anno_text(gwaspqtloxoinflammation$SNP_mdd,gp = gpar(fontsize = 8, fontface = "bold"), just = "right", location = unit(1, "npc")),
ha_gwas_ukbd_oxinfl = rowAnnotation(ukbd = anno_text(gwaspqtloxoinflammation$SNP_ukbb,gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                        #...= anno_simple(gwaspqtloxoinflammation$GWAS_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvaluemdd_oxinfl<0.0005,"***","**")),
                        .= anno_simple(gwaspqtloxoinflammation$GWAS_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvalueukbd_oxinfl<0.0005,"***","**")))
                       # ....= anno_simple(gwaspqtloxoinflammation$GWAS_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvaluenull_oxinfl<0.0005,"***","**")),
ha_pqtl_sbc_oxinfl = rowAnnotation(sbc = anno_text(gwaspqtloxoinflammation$SNP_pQTL_SBC,gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                        .= anno_simple(gwaspqtloxoinflammation$pQTL_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvaluepqtlsbc<0.0001,"***","*")))
ha_pqtl_fcc_oxinfl = rowAnnotation(fcc = anno_text(gwaspqtloxoinflammation$SNP_pQTL_FCC,gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                        .= anno_simple(gwaspqtloxoinflammation$pQTL_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvaluepqtlfcc<0.0001,"***","*")))



#draw(ha_gwas_pqtl_oxinfl)

#===================================================================================================================
#
# --------------- heatmap for dmr-proteins -------------------------
#
dmrprotein = read.csv("DMR inflammation angiogenesis DMR-protein_v4_good.csv", row.names = 1)

# subset the DMR file for the inflammatory and oxidative stress pathways
dmroxoinflammation = dmrprotein[dmrprotein$pathway=="inflammatory response"|dmrprotein$pathway=="oxidative stress",]
dim(dmroxoinflammation)

dmrprotein_matrix = as.matrix(dmroxoinflammation[,2:4])

dmrprotein_matrix[is.na(dmrprotein_matrix)] <- 0
#colnames(woundhealing) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_dmrprotein = rowAnnotation(dmrprotein_anno = anno_text(dmroxoinflammation$chr_TSS1500_TSS200, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun_dmr = colorRamp2(c(-0.015, 0, 0.01), c("green", "white", "magenta"))
htmp_dmr_oxoinflam = Heatmap(dmrprotein_matrix, 
                   show_row_names = F, 
                   column_title ="",#cis-regulatory DMRs",
                   column_title_side = c("top"),
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   column_title_rot = 0,
                   show_heatmap_legend = F,
                   rect_gp = gpar(col = "white", lwd = 2),
                   cluster_rows = F,
                   cluster_columns = F, 
                   show_column_names = F,
                   right_annotation = ha_dmrprotein, 
                   col = col_fun_dmr, column_names_rot = 45,
                   width = unit(10, "mm"))

#draw(htmp_dmr_oxoinflam)


#==========================================================================================================================
#
#------------ protein and metabolites heatmap ------------------
#
#

# read metabolites data for the combined pathways
proteinmetabolite = read.csv("metabolites inflammation angiogenesis proteins-metabolites_v4_good.csv", row.names = 1)

View(proteinmetabolite)

# subset the metabolite data for oxidative stress and inflammatory response
metaboliteinflammationoxidative = proteinmetabolite[proteinmetabolite$pathway=="inflammatory response"|proteinmetabolite$pathway=="oxidative stress",]
View(metaboliteinflammationoxidative)


# data frame to matrix
protmetabol_matrix = as.matrix(metaboliteinflammationoxidative[,2:6])

# change NA to zeros
protmetabol_matrix[is.na(protmetabol_matrix)] <- 0

# text annotation for the metabolite heatmap
ha_protmetabol = rowAnnotation(dmrprotein_anno = anno_text(metaboliteinflammationoxidative$metabolites, gp = gpar(fontsize = 8, fontface = "bold"), just = "right", location = unit(1, "npc")))

col_fun_metabol = colorRamp2(c(-0.5, 0, 0.4), c("dodgerblue1", "white", "deeppink1"))
htmp_metaboxoinflam = Heatmap(protmetabol_matrix, 
                          show_row_names = F, 
                          column_title = "",#significant metabolites",
                          column_title_side = c("top"),
                          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          column_title_rot = 0,
                          show_heatmap_legend = F,
                          rect_gp = gpar(col = "white", lwd = 2),
                          cluster_rows = F,
                          cluster_columns = F, 
                          show_column_names = F,
                          left_annotation = ha_protmetabol, 
                          col = col_fun_metabol, column_names_rot = 45,
                          width = unit(15, "mm"))

#################################################################################################################
#===============================================================================================================
#
# ---- heatmap for the oxidative stress and inflammatory response ---------------
#
#=============================================================================================================

combinedpathways = read.csv("heatmap combined pathways proteins postmortem sbc marin_good.csv", row.names = 1)

dim(combinedpathways)
View(combinedpathways)

names(combinedpathways)

oxidoinflammation = combinedpathways[combinedpathways$pathway=="inflammatory response"|combinedpathways$pathway=="oxidative stress",]

oxidoinflammation_matrix = as.matrix(oxidoinflammation[,-30])


#o1 = seriate(dist(oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,20,24,27)]), method = "GW")
#o2 = seriate(dist(t(oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,21,25,28)])), method = "GW")

o1 = seriate(dist(oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,20)]), method = "GW")
o2 = seriate(dist(t(oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,21)])), method = "GW")

#small_mat = oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,20,24,27)]
#small_mat_pv = oxidoinflammation_matrix[,c(6,10,8,9,7,23,17,19,21,25,28)]

small_mat = oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,20)]
small_mat_pv = oxidoinflammation_matrix[,c(6,10,8,9,7,23,17,19,21)]

#         working code for heatmap 

library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp1 = Heatmap(small_mat, name = "log2FC", col = col_fun,
                show_heatmap_legend = F,
                column_km = 2, column_gap = unit(2, "mm"),
                show_column_dend = F,
                show_row_dend = F,
                cluster_columns = F,
                row_title = NULL,
                column_title = NULL,#c("protein", "mRNA"),
                #cluster_rows = T,
                #left_annotation = htmp_metaboxoinflam,
                #left_annotation = ha_miR,
                #right_annotation = ha_gwas_pqtl_oxinfl,
                column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                row_names_side = c("left"),
                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_names_rot = 45,
                #border_gp = gpar(col = "yellow", lty = 2, lwd = 2),
                row_km = 2,
                rect_gp = gpar(col = "white", lwd = 2),
                #top_annotation = HeatmapAnnotation(modality = c(rep("protein",5),rep("mRNA",4) )),
                layer_fun = function(j, i, x, y, w, h, fill) {
                  # restore_matrix() is explained after this chunk of code
                  ind_mat = restore_matrix(j, i, x, y)
                  for(ir in seq_len(nrow(ind_mat))) {
                    # start from the second column
                    for(ic in seq_len(ncol(ind_mat))) {
                      ind = ind_mat[ir, ic] # previous column
                      v = small_mat_pv[i[ind], j[ind]]
                      
                      grid.points(x[ind], y[ind], 
                                  pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(1, "mm"))
                      
                    }
                  }
                }
                
)

#?Heatmap
col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2 (fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4, "cm"))

lgd2 = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "p-value", type = "points", title_position =  "topleft", 
              grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
              pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")
pd_ldgs = packLegend(lgd1, lgd2)

ht_opt(heatmap_column_names_gp = gpar(fontface = "bold"), 
       heatmap_column_title_gp = gpar(fontsize = 10),
       #legend_border = "black",
       #border_gp = gpar(col = "black", lty = 2),
       annotation_border = NULL #TRUE
)

draw(htmp_metaboxoinflam + htmp1 + ha_miR_oxoinflm + htmp_dmr_oxoinflam + ha_gwas_ptsd_oxinfl + ha_gwas_ukbd_oxinfl + 
       ha_pqtl_sbc_oxinfl + ha_pqtl_fcc_oxinfl, ht_gap = unit(c(3,3,3,3,3,3,3), "mm"), 
     main_heatmap = "log2FC", auto_adjust = FALSE)#, rect_gp = gpar(col = "white", lwd = 2), border_gp = gpar(col = "gray", lty = 2))#, cluster_rows= T)#, ht_gap = unit(1, "cm"))#, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")
#draw(htmp1 + htmp_dmr_oxoinflam + ha_gwas_pqtl_oxinfl, auto_adjust = FALSE, cluster_rows= T)#, ht_gap = unit(1, "cm"))#, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")

ht_opt(RESET = TRUE)

# draw the heatmaps and annotations
draw(htmp_responsetowounding + htmp_dmr + htmp_woundhealing + htmp_inflam + htmp_hemostasis + htmp_metabolite,  auto_adjust = FALSE)#, heatmap_legend = packaged_legends, heatmap_legend_side = NULL,  adjust_annotation_extension = TRUE, 
#merge_legend = TRUE)

###########################################################################################################

# ht_list_main = htmp_responsetowounding  %v%  htmp_vasculaturedevelopment
# htmp_dmr_list = htmp_dmr %v% htmp_va_dmr
# htmp_woundhealing_angiogenesis_list = htmp_woundhealing %v% htmp_angiogenesis
# htmp_inflam_heartdevelop_list = htmp_inflam %v% htmp_heartdevelop
# htmp_hemostasis
# htmp_metabolites_list = htmp_metabolite %v% htmp_va_metabolite

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2(fold change)", title_position = "topcenter", direction = "horizontal",
              legend_height = unit(4, "cm"))#, at = c(-1, -0.5, 0, 0.5, 1), labels = c("-1", "-0.5", "0", "0.5", "1"))

draw(lgd1)

ggplot_pathway = Legend(title = "pathway\nenrichment\nq-value", col = ggplot_lgd_col, at = c(1, 2, 3, 4), 
                        labels = c("0.01", "1E^-5", "1E^-10", "1E^-15"))

#lgd_deps = Legend(title = "num.of.DEPs", pch = 1, type = "points", size = sort(unit((count_fdr$count)/5,"mm")), labels = sort(count_fdr$count))
lgd_significant = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "protein/mRNA\np-value", type = "points", title_position =  "leftcenter", 
                         grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
                         pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 9), legend_gp = gpar(col = c(1,7)),  background = "white")
draw(lgd_significant)


library(circlize)
fold_col_fun2 = colorRamp2(c(-2, 0, 2), c("cyan", "white", "orange")) 
col2mir = colorRamp2(c(-2, 0, 2), c("cyan", "white", "orange"))
col3 = colorRamp2(breaks = c(-2,0,2),colors=c("cyan", "white", "orange"))
lgd_miR = Legend(col3, title = "miRNA\nlog2(fold change)",  at = c(-2, -1, 0, 1, 2), labels = c("-2", "-1", "0", "1", "2"),
                  title_position = "topcenter",direction = "horizontal",
                 legend_height = unit(2, "cm"))

draw(lgd_miR)

lgd_dmr = Legend(title = "hypo-methylated  hyper-methylated", title_position = "leftcenter-rot",
                 legend_height = unit(6, "cm"), col_fun_dmr, at = c(-0.015,  0,  0.015), 
                 labels = c("","" ,""))

draw(lgd_dmr)

lgd_sig_mir = Legend(pch = "*", type = "points", labels = "< 0.05", title = "miR: p-value", title_position =  "topleft",
                     size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")

draw(lgd_sig_mir)

lgd_gwas = Legend(pch = c("*", "**","***"), type = "points", labels = c("5e-3 to <1e-2", "5e-4 to <5e-3", "<5e-4"), title = "GWAS\np-value", title_position =  "leftcenter",
                     size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,1,1)),  background = "white")

draw(lgd_gwas)


lgd_pqtl = Legend(pch = c("*","***"), type = "points", labels = c("1e-4 to <1e-2", "<1e-4"), title = "pQTL\np-value", title_position =  "leftcenter",
                       size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,1,1)),  background = "white")

draw(lgd_pqtl)


lgd_count = Legend(labels = c(10,20,30,40,50), title = "pathway\n(protein count)", type = "points", title_position =  "topleft", 
                   grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
                   pch = 16, size = unit(2:6, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = 1),  background = "white")


#package the legends (annotation and heatmap legends)
packaged_legends = packLegend(lgd_sig_mir, lgd_miR,  lgd1, lgd_significant, lgd_dmr, ggplot_pathway, lgd_count,  row_gap = unit(1, "cm"), direction = "horizontal")#, lgd_deps)
dev.off()
draw(packaged_legends)
# draw the heatmaps and annotations
draw(htmp_vasculaturedevelopment + htmp_va_dmr + htmp_angiogenesis + htmp_heartdevelop + htmp_vasculogenesis + htmp_va_metabolite,  auto_adjust = FALSE)#, heatmap_legend = packaged_legends, heatmap_legend_side = NULL,  adjust_annotation_extension = TRUE, 
#merge_legend = TRUE)

#==================================


# decorate the bar graph (describe what it is)
decorate_annotation("..", {
  grid.text("number of \nsignificant \nproteins", -0.4, 0.6, default.units = "npc", just = "right", gp = gpar(fontsize = 8, fontface= "bold", fontfamily="Helvetica"))
  #grid.rect(x = 0.6, y = 0.9, width = 0.5, height = 0.5, draw(lgd_protein_direct_of_regulat), default.units = "native")#, gp = gpar(fill = "tomato1"))# gp = gpar(fill = c("tomato1", "slateblue1")), 
  grid.rect(x = 0.85, y = 0.72,  gp = gpar(fill = "tomato1"),width = 0.07, height = 0.1,)
  grid.rect(x = 0.85, y = 0.6,  gp = gpar(fill = "slateblue1"),width = 0.07, height = 0.1,)
  grid.text("up",0.78,0.76,  default.units = "npc", just = "right", gp = gpar(fontsize = 6, fontface="bold"))
  grid.text("down",0.78,0.62,  default.units = "npc", just = "right", gp = gpar(fontsize = 6, fontface="bold") )
  #grid.text("regulation",0.6, 0.92, default.units = "npc", just = "right" )
})


#=======================================
# decorate the pathwayBP for the ggplot_anno bottom annotation

decorate_annotation("pathwayBP", {
  grid.text("vasculature developmment", 1.05, 0.40, default.units = "npc", just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  grid.text("angiogenesis", 1.05, 0.56, default.units = "npc", just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  grid.text("heart development", 1.05, 0.72, default.units = "npc", just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  grid.text("vasculogenesis", 1.05, 0.88, default.units = "npc", just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  
})


#=========================


# lable the main heatmap, since all of the proteins shown in the heatmap are involved in the pathway
# "vasculature development"
decorate_annotation("pathways", {
  grid.text("vasculature development", -0.6, 0.40, default.units = "npc", just = "left", gp = gpar(fontsize = 10, fontface="bold"))
  
})

#=====================================

# label the methylation heatmap
decorate_annotation("va_dmrprotein_anno", {
  grid.text("differentially methylated \nregions (SBC Training)", -0.2, 1.05, default.units = "npc", just = "left" , gp = gpar(fontsize = 11, fontface="bold"))
  
})


#==================================


#decorate the heatmap orannotation of angiogenesis
decorate_annotation("angiogenesis", {
  grid.text("angiogenesis", -1.8, 1.05, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface = "bold"))
  
})

#==============================

#decorate/label the heart development and text annotation
decorate_annotation("heartdevelopment", {
  grid.text("heart \ndevelopment", -1.4, 1.05, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface = "bold"))
  
})

#decorate/label the heart development and text annotation
decorate_annotation("vasculogenesis", {
  grid.text("vasculogenesis", -1.2, 1.05, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface = "bold"))
  
})

#===============================================


# label the metabolites heatmap and text annotation
decorate_annotation("vas_protmetabol_anno", {
  grid.text("significant metabolites\n(SBC Training)", -0.15, 1.05, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface="bold"))
  
})


#vas_protmetabol_anno = anno_text(vas_proteinmetabolite

#====================================================================

# decorate_annotation(".",{
#   grid.text("* < 0.05", -2.18, 0.98, default.units = "npc", just = "left")
#   grid.text("p-value", -3.1, 0.98, default.units = "npc", just = "right")
# })
decorate_annotation(".",{
  grid.text("* < 0.05", -2.16, -0.77, default.units = "npc", just = "left")
  grid.text("p-value", -3.1, -0.77, default.units = "npc", just = "right")
})

decorate_annotation(".",{
  grid.text("* < 0.05", -2.16, -0.25, default.units = "npc", just = "left")
  grid.text("p-value", -3.1, -0.25, default.units = "npc", just = "right")
})


#==============================================================

list_components()



# heatmap lists for vertical stacking
ht_list_main = htmp_responsetowounding  %v%  htmp_vasculaturedevelopment
htmp_dmr_list = htmp_dmr %v% htmp_va_dmr
htmp_woundhealing_angiogenesis_list = htmp_woundhealing %v% htmp_angiogenesis
htmp_inflam_heartdevelop_list = htmp_inflam %v% htmp_heartdevelop
htmp_hemostasis
htmp_metabolites_list = htmp_metabolite %v% htmp_va_metabolite



draw(ht_list_main, ht_gap = unit(6,  "mm"))

# draw the heatmaps and annotations
draw(ht_list_main + htmp_dmr_list + htmp_woundhealing_angiogenesis_list + htmp_inflam_heartdevelop_list + htmp_metabolites_list,  auto_adjust = FALSE, heatmap_legend = packaged_legends, heatmap_legend_side = NULL,  adjust_annotation_extension = TRUE, 
     merge_legend = TRUE)





################################====================================


#decorate the heatmap orannotation of angiogenesis
decorate_annotation("angiogenesis", {
  grid.text("angiogenesis", -1.2, 1.05, default.units = "npc", just = "left")
  
})

#decorate/label the heart development and text annotation
decorate_annotation("heartdevelopment", {
  grid.text("heart \ndevelopment", -1, 1.06, default.units = "npc", just = "left")
  
})


# lable the main heatmap, since all of the proteins shown in the heatmap are involved in the pathway
# "vasculature development"
decorate_annotation("pathways", {
  grid.text("vasculature development", -0.1, 0.4, default.units = "npc", just = "left")
  
})




########################################################################################################
#======================================================================================================
#
# ----------------- Annotations GWAS miR pQTL 
#          heatmap for angiogenesis and heart development  --------------------------------
#
#=====================================================================================================

library(dplyr)

gwaspqtlmirpathways = read.csv("GWAS-pQTL-miR pathways_v4_good.csv", row.names = 1)
#mirprotein = mirprotein[-c(35,37,69),]
# subsetting the data based on the common proteins between the two sets of heatmaps
#mirprotein =  mirprotein[!(row.names(mirprotein) %in% common16proteins_toberemoved),]

View(gwaspqtlmirpathways)
gwaspqtlmircardioangiogenesis = gwaspqtlmirpathways[gwaspqtlmirpathways$pathway=='angiogenesis'|gwaspqtlmirpathways$pathway=='heart development',]

pvaluemir_cardangiog = gwaspqtlmircardioangiogenesis$P.Value_mir
pvalueptsd_cardangiog = gwaspqtlmircardioangiogenesis$P_GWAS_pgc_ptsd
pvaluemdd_cardangiog = gwaspqtlmircardioangiogenesis$P_GWAS_pgc_mdd
pvalueukbd_cardangiog = gwaspqtlmircardioangiogenesis$P_GWAS_ukbb_dep
pvaluenull_cardangiog = gwaspqtlmircardioangiogenesis$nopvalue
pvaluepqtlsbc_cardangiog = gwaspqtlmircardioangiogenesis$P_pQQTL_SBC
pvaluepqtlfcc_cardangiog = gwaspqtlmircardioangiogenesis$P_pQQTL_FCC



fold_col_fun = colorRamp2(c(-2, 0, 2), c("cyan", "white", "orange")) 
ha_miR_cardangiog = rowAnnotation(.= anno_simple(gwaspqtlmircardioangiogenesis$logFC_mir,  col = fold_col_fun, na_col = "white", pch = ifelse((pvaluemir_cardangiog<0.05 & !is.na(pvaluemir_cardangiog)),"*","")),
                                mir = anno_text(gwaspqtlmircardioangiogenesis$miR_ID,gp = gpar(fontsize = 10, fontface = "bold"), just = "left", location = unit(0.05, "npc")))

 
ha_gwas_ptsd_cardangiog = rowAnnotation(ptsd = anno_text(gwaspqtlmircardioangiogenesis$SNP_ptsd,gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                                    .= anno_simple(gwaspqtlmircardioangiogenesis$GWAS_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvalueptsd_cardangiog<0.005,"**","*")))
#mdd = anno_text(gwaspqtloxoinflammation$SNP_mdd,gp = gpar(fontsize = 8, fontface = "bold"), just = "right", location = unit(1, "npc")),
ha_gwas_mdd_cardangiog = rowAnnotation(mdd = anno_text(gwaspqtlmircardioangiogenesis$SNP_mdd, gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                                        .= anno_simple(gwaspqtlmircardioangiogenesis$GWAS_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvaluemdd_cardangiog<0.005,"**","*")))
ha_gwas_ukbd_cardangiog = rowAnnotation(ukbd = anno_text(gwaspqtlmircardioangiogenesis$SNP_ukbb,gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                                    .= anno_simple(gwaspqtlmircardioangiogenesis$GWAS_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvalueukbd_cardangiog<0.0005,"***","**")))
# ....= anno_simple(gwaspqtlmircardioangiogenesis$GWAS_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvaluenull_oxinfl<0.0005,"***","**")),
ha_pqtl_sbc_cardangiog = rowAnnotation(sbc = anno_text(gwaspqtlmircardioangiogenesis$SNP_pQTL_SBC,gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                                   .= anno_simple(gwaspqtlmircardioangiogenesis$pQTL_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvaluepqtlsbc_cardangiog<0.0001,"***","*")))
ha_pqtl_fcc_cardangiog = rowAnnotation(fcc = anno_text(gwaspqtlmircardioangiogenesis$SNP_pQTL_FCC,gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")),
                                   .= anno_simple(gwaspqtlmircardioangiogenesis$pQTL_log2FC,  col = fold_col_fun, na_col = "white", pch = if_else(pvaluepqtlfcc_cardangiog<0.0001,"***","*")))



#draw(ha_gwas_pqtl_cardangiog)

#===================================================================================================================
#
# --------------- heatmap for dmr-proteins -------------------------
#
dmrprotein = read.csv("DMR inflammation angiogenesis DMR-protein_v4_good.csv", row.names = 1)

# subset the DMR file for the inflammatory and oxidative stress pathways
dmrcardioangiogenesis = dmrprotein[dmrprotein$pathway=="angiogenesis"|dmrprotein$pathway=="heart development",]
View(dmrcardioangiogenesis)

dmrprotein_cardangiog_matrix = as.matrix(dmrcardioangiogenesis[,2:3])

dmrprotein_cardangiog_matrix[is.na(dmrprotein_cardangiog_matrix)] <- 0
#colnames(woundhealing) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_dmrprotein_cardangiog = rowAnnotation(dmrprotein_cardangiog = anno_text(dmrcardioangiogenesis$chr_TSS1500_TSS200, gp = gpar(fontsize = 10, fontface = "bold"), just = "left", location = unit(0, "npc")))


col_fun_dmr = colorRamp2(c(-0.015, 0, 0.01), c("green", "white", "magenta"))
htmp_dmr_cardangiog = Heatmap(dmrprotein_cardangiog_matrix, 
                             show_row_names = F, 
                             column_title ="",#cis-regulatory DMRs",
                             column_title_side = c("top"),
                             column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                             column_title_rot = 0,
                             show_heatmap_legend = F,
                             rect_gp = gpar(col = "white", lwd = 2),
                             cluster_rows = F,
                             cluster_columns = F, 
                             show_column_names = F,
                             right_annotation = ha_dmrprotein_cardangiog, 
                             col = col_fun_dmr, column_names_rot = 45,
                             width = unit(8, "mm"))




#==========================================================================================================================
#
#------------ protein and metabolites heatmap ------------------
#
#

# read metabolites data for the combined pathways
proteinmetabolite = read.csv("metabolites inflammation angiogenesis proteins-metabolites_v4_good.csv", row.names = 1)

View(proteinmetabolite)

# subset the metabolite data for oxidative stress and inflammatory response
metabolitecardioangiogenesis = proteinmetabolite[proteinmetabolite$pathway=="angiogenesis"|proteinmetabolite$pathway=="heart development",]
View(metabolitecardioangiogenesis)


# data frame to matrix
protmetabol_cardangiog_matrix = as.matrix(metabolitecardioangiogenesis[,2:4])

# change NA to zeros
protmetabol_cardangiog_matrix[is.na(protmetabol_cardangiog_matrix)] <- 0

# text annotation for the metabolite heatmap
ha_protmetabol_cardangiog = rowAnnotation(dmrprotein_cardangiog = anno_text(metabolitecardioangiogenesis$metabolites, gp = gpar(fontsize = 10, fontface = "bold"), just = "right", location = unit(1, "npc")))

col_fun_metabol = colorRamp2(c(-0.5, 0, 0.4), c("dodgerblue1", "white", "deeppink1"))
htmp_metab_cardangiog = Heatmap(protmetabol_cardangiog_matrix, 
                              show_row_names = F, 
                              column_title = "",#significant metabolites",
                              column_title_side = c("top"),
                              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                              column_title_rot = 0,
                              show_heatmap_legend = F,
                              rect_gp = gpar(col = "white", lwd = 2),
                              cluster_rows = F,
                              cluster_columns = F, 
                              show_column_names = F,
                              left_annotation = ha_protmetabol_cardangiog, 
                              col = col_fun_metabol, column_names_rot = 45,
                              width = unit(10, "mm"))


#################################################################################################################
#===============================================================================================================
#
# ---- heatmap for the angiogenesis and heart development  ---------------
#
#=============================================================================================================

combinedpathways = read.csv("heatmap combined pathways proteins postmortem sbc marin_good.csv", row.names = 1)

dim(combinedpathways)
View(combinedpathways)

names(combinedpathways)

cardioangiogenesis = combinedpathways[combinedpathways$pathway=="angiogenesis"|combinedpathways$pathway=="heart development",]

cardioangiogenesis_matrix = as.matrix(cardioangiogenesis[,-30])


#o1 = seriate(dist(oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,20,24,27)]), method = "GW")
#o2 = seriate(dist(t(oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,21,25,28)])), method = "GW")

o1_ang = seriate(dist(cardioangiogenesis_matrix[,c(1,5,3,4,2,22,16,18,20)]), method = "GW")
o2_ang = seriate(dist(t(cardioangiogenesis_matrix[,c(1,5,3,4,2,22,16,18,21)])), method = "GW")

#small_mat = oxidoinflammation_matrix[,c(1,5,3,4,2,22,16,18,20,24,27)]
#small_mat_pv = oxidoinflammation_matrix[,c(6,10,8,9,7,23,17,19,21,25,28)]

small_mat_ang = cardioangiogenesis_matrix[,c(1,5,3,4,2,22,16,18,20)]
small_mat_pv_ang = cardioangiogenesis_matrix[,c(6,10,8,9,7,23,17,19,21)]

#         working code for heatmap 

library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp1_cardangiog = Heatmap(small_mat_ang, name = "log2FC", col = col_fun,
                show_heatmap_legend = F,
                column_km = 2, column_gap = unit(2, "mm"),
                show_column_dend = F,
                show_row_dend = F,
                cluster_columns = F,
                row_title = NULL,
                column_title = NULL,#c("protein", "mRNA"),
                #cluster_rows = T,
                #left_annotation = htmp_metaboxoinflam,
                #left_annotation = ha_miR,
                #right_annotation = ha_gwas_pqtl_oxinfl,
                column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                row_names_side = c("left"),
                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_names_rot = 45,
                #border_gp = gpar(col = "yellow", lty = 2, lwd = 2),
                row_km = 2,
                rect_gp = gpar(col = "white", lwd = 2),
                #top_annotation = HeatmapAnnotation(modality = c(rep("protein",5),rep("mRNA",4) )),
                layer_fun = function(j, i, x, y, w, h, fill) {
                  # restore_matrix() is explained after this chunk of code
                  ind_mat = restore_matrix(j, i, x, y)
                  for(ir in seq_len(nrow(ind_mat))) {
                    # start from the second column
                    for(ic in seq_len(ncol(ind_mat))) {
                      ind = ind_mat[ir, ic] # previous column
                      v = small_mat_pv_ang[i[ind], j[ind]]
                      
                      grid.points(x[ind], y[ind], 
                                  pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(1, "mm"))
                      
                    }
                  }
                }
                
)

#?Heatmap
col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "log2 (fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4, "cm"))

lgd2 = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "p-value", type = "points", title_position =  "topleft", 
              grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
              pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")
pd_ldgs = packLegend(lgd1, lgd2)

ht_opt(heatmap_column_names_gp = gpar(fontface = "bold"), 
       heatmap_column_title_gp = gpar(fontsize = 10),
       #legend_border = "black",
       border_gp = gpar(col = "black", lty = 2),
       annotation_border = NULL #TRUE
)

draw(htmp_metab_cardangiog + htmp1_cardangiog + ha_miR_cardangiog + htmp_dmr_cardangiog + ha_gwas_ptsd_cardangiog + ha_gwas_mdd_cardangiog + ha_gwas_ukbd_cardangiog + 
       ha_pqtl_sbc_cardangiog + ha_pqtl_fcc_cardangiog, ht_gap = unit(c(3,3,3,3,3,3,3,3), "mm"), 
     main_heatmap = "log2FC", auto_adjust = FALSE)#, rect_gp = gpar(col = "white", lwd = 2), border_gp = gpar(col = "gray", lty = 2))#, cluster_rows= T)#, ht_gap = unit(1, "cm"))#, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")
#draw(htmp1 + htmp_dmr_oxoinflam + ha_gwas_pqtl_oxinfl, auto_adjust = FALSE, cluster_rows= T)#, ht_gap = unit(1, "cm"))#, heatmap_legend = pd_ldgs, heatmap_legend_side = "right")

ht_opt(RESET = TRUE)

# draw the heatmaps and annotations
draw(htmp_responsetowounding + htmp_dmr + htmp_woundhealing + htmp_inflam + htmp_hemostasis + htmp_metabolite,  auto_adjust = FALSE)#, heatmap_legend = packaged_legends, heatmap_legend_side = NULL,  adjust_annotation_extension = TRUE, 
#merge_legend = TRUE)



# ht_list_main = htmp_responsetowounding  %v%  htmp_vasculaturedevelopment
# htmp_dmr_list = htmp_dmr %v% htmp_va_dmr
# htmp_woundhealing_angiogenesis_list = htmp_woundhealing %v% htmp_angiogenesis
# htmp_inflam_heartdevelop_list = htmp_inflam %v% htmp_heartdevelop
# htmp_hemostasis
# htmp_metabolites_list = htmp_metabolite %v% htmp_va_metabolite

# decorate the pathwayBP for the ggplot_anno bottom annotation
decorate_annotation("pathwayBP", {
  grid.text("response to wounding", 1.05, 0.37, default.units = "npc", just = "left")
  grid.text("wound healing", 1.05, 0.53, default.units = "npc", just = "left")
  grid.text("inflammatory response", 1.05, 0.69, default.units = "npc", just = "left")
  grid.text("hemostasis (coagulation)", 1.05, 0.85, default.units = "npc", just = "left")
  
})


# lable the main heatmap, since all of the proteins shown in the heatmap are involved in the pathway
# "response to wounding"
decorate_annotation("pathways", { 
  grid.text("response to wounding", -0.5, 0.4, default.units = "npc", just = "left", gp = gpar(fontsize = 10, fontface="bold"))
  
})

# decorate the bar graph (describe what it is)
decorate_annotation("..", {
  grid.text("number of \nsignificant \nproteins", -0.4, 0.63, default.units = "npc", just = "right", gp = gpar(fontsize = 8, fontface= "bold", fontfamily="Helvetica"))
  #grid.rect(x = 0.6, y = 0.9, width = 0.5, height = 0.5, draw(lgd_protein_direct_of_regulat), default.units = "native")#, gp = gpar(fill = "tomato1"))# gp = gpar(fill = c("tomato1", "slateblue1")), 
  grid.rect(x = 0.85, y = 0.72,  gp = gpar(fill = "tomato1"),width = 0.07, height = 0.1,)
  grid.rect(x = 0.85, y = 0.6,  gp = gpar(fill = "slateblue1"),width = 0.07, height = 0.1,)
  grid.text("up",0.78,0.76,  default.units = "npc", just = "right", gp = gpar(fontsize = 6, fontface="bold"))
  grid.text("down",0.78,0.62,  default.units = "npc", just = "right", gp = gpar(fontsize = 6, fontface="bold") )
  #grid.text("regulation",0.6, 0.92, default.units = "npc", just = "right" )
})

# label the methylation heatmap
decorate_annotation("dmrprotein_anno", {
  grid.text("differentially methylated \nregions (SBC Training)", -1.9, 1.08, default.units = "npc", just = "left" , gp = gpar(fontsize = 11, fontface="bold"))
  
})

#decorate the heatmap orannotation of wound healing
decorate_annotation("coagulation", {
  grid.text("wound healing", -1.5, 1.08, default.units = "npc", just = "left" , gp = gpar(fontsize = 11, fontface="bold"))
  
})

#decorate/label the inflammatory heatmap and text annotation
decorate_annotation("coagulation", {
  grid.text("inflammatory \n response", 0.7, 1.08, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface="bold"))
  
})

#label the hemostasis heatmap or annoation
decorate_annotation("coagulation", {
  grid.text("hemostasis\n(coagulation)", 2.6, 1.08, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface="bold"))
  
})

# label the metabolites heatmap and text annotation
decorate_annotation("coagulation", {
  grid.text("significant metabolites (SBC Training)", 4.5, 1.08, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface="bold"))
  
})





decorate_annotation(".",{
  grid.text("* < 0.05", -2.16, -0.77, default.units = "npc", just = "left")
  grid.text("p-value", -3.1, -0.77, default.units = "npc", just = "right")
})



#==============================================================================================
################################################################################################

# pQTL box plot for select proteins

# read the pQTL datasets for four SNPs and proteins
pqtldata = read.csv("SBC_SNP_protein_data_v5.csv")
View(pqtldata)

pqtldata2 = pqtldata[pqtldata$GeneName!="MAPKAPK3" & pqtldata$GeneName!="IL1RL1",]
levels(as.factor(pqtldata2$GeneName))
pqtldata3 = na.omit(pqtldata2)
library(ggplot2)
# Basic box plot
# p <- ggplot(pqtldata3, aes(x=as.factor(SNP), y=protein)) + 
#   geom_boxplot()
# p
# # Box plot with dot plot
# p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# # Box plot with jittered points
# # 0.2 : degree of jitter in x direction
# p + geom_jitter(shape=16, position=position_jitter(0.2))

# Change box plot colors by groups
p<-ggplot(pqtldata3, aes(x=as.factor(SNP), y=protein, fill = as.factor(SNP))) +
  scale_fill_brewer(palette="RdBu") + geom_boxplot() 
p2 = p + facet_grid(SNP_GeneName ~ Cohort, scales="free") + theme_minimal()
p2 + xlab("minor allels") + ylab("protein levels") + labs(fill = "allel") +
    theme(legend.position="bottom") + theme(legend.key.width  = unit(0.3, 'cm')) +
   theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text( size=14, face="bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.text=element_text(size=10, face = "bold"),
    legend.title =element_text(size=10, face = "bold"))
  #+ theme(legend.title = element_blank()) + 
  #theme(legend.position='none')
 
# to resize the ggplot legend 
theme(legend.key.height= unit(4, 'cm'),
      legend.key.width= unit(1, 'cm'))
#or both at the same time
theme(legend.key.width  = unit(0.2, 'cm'))





########################################################################################################################
#=============================================================================
#
# --- number of down and up regulated proteins for each cohort ad top annotaiton
#

deps_long = read.csv("01_001_hptsd_somalogic_DEPs_long.csv", row.names = 1)
names(deps_long)


# exclude the civilian and the female cohorts


deps_long_significant = deps_long[deps_long$pvalue<0.05, ]

deps_long_significant_up = deps_long_significant[deps_long_significant$logFC>0,]
deps_long_significant_down = deps_long_significant[deps_long_significant$logFC<0,]


cohort_logfc_down = xtabs(~ cohort, data = deps_long_significant_down)
#View(cohort_logfc_down)

cohort_logfc_up = xtabs(~ cohort, data = deps_long_significant_up)
#View(cohort_logfc_up)

deps_down_up_significant = cbind(cohort_logfc_down, cohort_logfc_up)
#View(deps_down_up_significant)

# exclude civilians and females
deps_nocivilianfemales = deps_down_up_significant[-c(2,3,7),]

rownames(deps_nocivilianfemales) = c( "SBC Training", "FCC Longitudinal", "FCC Validation", "FCC Subthreshold", "SBC Testing")

# file for bar plot
depdata = deps_nocivilianfemales[c(1,5,3,4,2),]
colnames(depdata) = c("down","up")

depdata = as.data.frame(depdata)
View(depdata)

depdata[,"cohort"] = c("SBC Training", "SBC Testing", "FCC Validation", "FCC Subthreshold", "FCC Longitudinal")

# reshape depdata from wide to long format
library(data.table)
depdata_long <- melt(setDT(depdata), id.vars = "cohort", variable.name = "regulation")
View(depdata_long)

# library(reshape2)
# depdata_long2 <- melt(depdata, id.vars = "cohort")
# View(depdata_long2)

# draw the barplots for the depdata_long
library(ggplot2)

depdata_long_sorted <- depdata_long[rev(order(depdata_long$regulation)),]
View(depdata_long_sorted)

p<-ggplot(depdata_long_sorted, aes(x=cohort, y=value, fill = regulation)) +
  geom_bar(stat="identity",  position=position_dodge2()) + theme_minimal() + 
  scale_x_discrete(limits=rev(levels(as.factor(depdata_long_sorted$cohort)))) +
  scale_fill_manual(values=c("gainsboro", "grey" ))+
  scale_color_manual(values=c("black", "black" )) +
  theme_classic()
p

p2 = p + xlab("") + ylab("number of significant proteins") + labs(fill = "") +
  theme(legend.position=c(0.7,0.9)) + theme(legend.key.size =  unit(0.5, 'cm')) #+ #theme(legend.title = element_blank()) + 
 #scale_fill_grey() + theme_classic()
  #scale_fill_brewer(palette="Reds") + theme_minimal()
  #scale_fill_manual(values=c("#9999CC", "#CC6666" ))+
  
p2 + theme(axis.text.x = element_text(face = "bold", color="black", angle = 45, vjust = 1, hjust=1),
           axis.text.y = element_text(face="bold"),
           legend.text = element_text(face = "bold"),
           axis.title.y = element_text(face="bold")) +
     geom_col(width = 0.1, position = position_dodge(0.1)) 

   
"gainsboro"
"lightblue1"
"palevioletred2"
"salmon"
"snow3"


##############################################################################
#
#-- intersect between pathway_protein list and gwas-pqtl supplementary table
#
#===========================================================================


setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/10282020_Somalogic_WoundHealing_toScienceWriter/04052022_Manuscript_for_MolecularCell_Submission/05102022_SummaryHeatmaps")
protein_pathway_list = read.csv("miR inflammation angiogenesis mir-protein_v4_good.csv")

proteins_pathway = protein_pathway_list[,1:2]
View(proteins_pathway)

setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/10282020_Somalogic_WoundHealing_toScienceWriter/04052022_Manuscript_for_MolecularCell_Submission/05222022_Files_for_Submission")
gwas_pqtl = read.csv("GWAS_and_pQTL_supplementary table.csv")

View(gwas_pqtl)


gwas_pqtl_pathway = merge(proteins_pathway, gwas_pqtl, by.x = "X", by.y = "GeneName")

View(gwas_pqtl_pathway)

write.csv(gwas_pqtl_pathway, "GWAS pQTL inflammation oxidative sress angiogenesis.csv")




#############################################################################################################################
###########################################################################################################################
#======================================================================================================================
######################################################################################################################
# Response to wounding: wound healing, inflammatory response, hemostasis for all cohorts except civilians and females
####################################################################################################################
#=================================================================================================================
#
# build the complex graph piece by piece for wound healing response to wounding inflammation and hemostasis
rm(list = ls())
#
# for all cohorts excluding the femals and civilians
#
#===============================================================================

# packages and libraries

library(seriation)
library(ComplexHeatmap)
library(circlize)
library(gridtext)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(GetoptLong)
library(methods)
library(reshape)


#========================================================
#
# ----  protein - metabolites interactome ----------
#
# -- reading the metabolite data along with the proein ids

rm(list =ls())
setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/10282020_Somalogic_WoundHealing_toScienceWriter/04052022_Manuscript_for_MolecularCell_Submission/05102022_SummaryHeatmaps")


proteinmetabolite = read.csv("metabolites inflammation angiogenesis proteins-metabolites_v4_good.csv", row.names = 1)

#proteinmetabolite = proteinmetabolite[-c(35,37,69),]
#View(proteinmetabolite)

# subsetting the data based on the common proteins between the two sets of heatmaps
#proteinmetabolite =  proteinmetabolite[!(row.names(proteinmetabolite) %in% common16proteins_toberemoved),]

View(proteinmetabolite)
metaboliteangiogenesisheartevelopment = proteinmetabolite[proteinmetabolite$pathway=="angiogenesis"|proteinmetabolite$pathway=="heart development",]
View(metaboliteangiogenesisheartevelopment)

metaboliteinflammationoxidative = proteinmetabolite[proteinmetabolite$pathway=="inflammatory response"|proteinmetabolite$pathway=="oxidative stress",]
View(metaboliteinflammationoxidative)


# data frame to matrix
protmetabol_matrix = as.matrix(metaboliteangiogenesisheartevelopment[,2:4])
protmetabol_matrix = as.matrix(metaboliteinflammationoxidative[,2:7])
# change NA to zeros
protmetabol_matrix[is.na(protmetabol_matrix)] <- 0

# text annotation for the metabolite heatmap
ha_protmetabol = rowAnnotation(dmrprotein_anno = anno_text(metaboliteinflammationoxidative$metabolites, gp = gpar(fontsize = 9, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun_metabol = colorRamp2(c(-0.5, 0, 0.4), c("dodgerblue1", "white", "deeppink1"))
htmp_metabolite = Heatmap(protmetabol_matrix, 
                          show_row_names = F, 
                          column_title = "",#"significant\nmetabolites",
                          column_title_side = c("top"),
                          column_title_gp = gpar(fontsize = 6, fontface = "bold"),
                          column_title_rot = 0,
                          show_heatmap_legend = F,
                          rect_gp = gpar(col = "white", lwd = 2),
                          cluster_rows = F,
                          cluster_columns = F, 
                          show_column_names = F,
                          right_annotation = ha_protmetabol, 
                          col = col_fun_metabol, column_names_rot = 45,
                          width = unit(10, "mm"))


draw(htmp_metabolite)

#===================================================================================================================
#
# --------------- dmr-proteins
#
dmrprotein = read.csv("DMR inflammation angiogenesis DMR-protein_v4_good.csv", row.names = 1)

dmroxoinflammation = dmrprotein[dmrprotein$pathway=="inflammatory response"|dmrprotein$pathway=="oxidative stress",]
View(dmroxoinflammation)
dim(dmroxoinflammation)

dmrprotein_matrix = as.matrix(dmroxoinflammation[,2:4])


dmrangiocardiac = dmrprotein[dmrprotein$pathway=="angiogenesis"|dmrprotein$pathway=="heart development",]
View(dmrangiocardiac)
dim(dmrangiocardiac)

dmrprotein_matrix = as.matrix(dmrangiocardiac[,2:4])

#woundhealing_pv = wound_inflam_coag_proteinspathways[,22:24]
dmrprotein_matrix[is.na(dmrprotein_matrix)] <- 0
#colnames(woundhealing) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_dmrprotein = rowAnnotation(dmrprotein_anno = anno_text(dmrangiocardiac$chr_TSS1500_TSS200, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun_dmr = colorRamp2(c(-0.015, 0, 0.01), c("green", "white", "magenta"))
htmp_dmr = Heatmap(dmrprotein_matrix, 
                   show_row_names = F, 
                   column_title ="", #"differentially methylated promoter regions",
                   column_title_side = c("top"),
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   column_title_rot = 0,
                   show_heatmap_legend = F,
                   rect_gp = gpar(col = "white", lwd = 2),
                   cluster_rows = F,
                   cluster_columns = F, 
                   show_column_names = F,
                   right_annotation = ha_dmrprotein, 
                   col = col_fun_dmr, column_names_rot = 45,
                   width = unit(10, "mm"))

draw(htmp_dmr)

#=====================================================================================================
#
# -- protein -- pathways 
#

wound_inflam_coag_proteinspathways = read.csv("wounding_inflammation_wounding inflammation coagulation proteins-pathways good_v2.csv", row.names = 1)
names(wound_inflam_coag_proteinspathways)
#wound_inflam_coag_proteinspathways = wound_inflam_coag_proteinspathways[-c(35,37,69),]

# subsetting the data based on the common proteins between the two sets of heatmaps
#wound_inflam_coag_proteinspathways =  wound_inflam_coag_proteinspathways[!(row.names(wound_inflam_coag_proteinspathways) %in% common16proteins_toberemoved),]

View(wound_inflam_coag_proteinspathways)

responsetowounding = as.matrix(wound_inflam_coag_proteinspathways[,c(1,4:6,8)])
responsetowounding_pv = wound_inflam_coag_proteinspathways[,c(9,12:14,16)]
colnames(responsetowounding)
colnames(responsetowounding_pv)

o1 = seriate(dist(responsetowounding), method = "GW")
o2 = seriate(dist(t(responsetowounding)), method = "GW")

woundhealing = as.matrix(wound_inflam_coag_proteinspathways[,19:21])
colnames(woundhealing)
#woundhealing_pv = wound_inflam_coag_proteinspathways[,22:24]
woundhealing[is.na(woundhealing)] <- 0
colnames(woundhealing) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_woundhealing = rowAnnotation(coagulation = anno_text(wound_inflam_coag_proteinspathways$Wound_healing, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp_woundhealing = Heatmap(woundhealing, 
                            show_row_names = F, 
                            column_title ="" ,#"wound \nhealing",
                            column_title_side = c("top"),
                            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                            column_title_rot = 0,
                            show_heatmap_legend = F,
                            rect_gp = gpar(col = "white", lwd = 2),
                            cluster_rows = F,
                            cluster_columns = F, 
                            show_column_names = T,
                            na_col = "white",
                            column_names_gp = gpar(fontsize = 5, fontface = "bold"),
                            right_annotation = ha_woundhealing, 
                            col = col_fun, column_names_rot = 45,
                            width = unit(12, "mm"))


#inflammation = inflammation[-67,]
inflammation = as.matrix(wound_inflam_coag_proteinspathways[,27:29])
View(inflammation)
colnames(inflammation)
#inflammation_pv = wound_inflam_coag_proteinspathways[,30:32]
inflammation[is.na(inflammation)] <- 0
colnames(inflammation) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_inflam = rowAnnotation(coagulation = anno_text(wound_inflam_coag_proteinspathways$Inflammatory_response, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp_inflam = Heatmap(inflammation, 
                      show_row_names = F, 
                      column_title ="", #"inflammatory\nresponse",
                      column_title_side = c("top"),
                      #column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      column_title_rot = 0,
                      show_heatmap_legend = F,
                      rect_gp = gpar(col = "white", lwd = 2),
                      cluster_rows = F,
                      cluster_columns = F, 
                      show_column_names = T,
                      na_col = "white",
                      column_names_gp = gpar(fontsize = 5, fontface = "bold"),
                      right_annotation = ha_inflam, 
                      col = col_fun, column_names_rot = 45,
                      width = unit(12, "mm"))



hemostasis = as.matrix(wound_inflam_coag_proteinspathways[,35:37])
#hemostasis = hemostasis[-c(35,37),]
#hemostasis = hemostasis[-67,]
View(hemostasis)
colnames(hemostasis)
hemostasis[is.na(hemostasis)] <- 0
colnames(hemostasis) = c("SBC Training" , "SBC Testing",  "FCC Validation")


ha_hemostat = rowAnnotation(coagulation = anno_text(wound_inflam_coag_proteinspathways$hemostasis, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp_hemostasis = Heatmap(hemostasis, 
                          show_row_names = F, 
                          column_title = "",#"hemostasis\n(coagulation)",
                          column_title_side = c("top"),
                          #column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          column_title_rot = 0,
                          show_heatmap_legend = F,
                          rect_gp = gpar(col = "white", lwd = 2),
                          cluster_rows = F,
                          cluster_columns = F, 
                          show_column_names = T,
                          na_col = "white",
                          column_names_gp = gpar(fontsize = 5, fontface = "bold"),
                          right_annotation = ha_hemostat, 
                          col = col_fun, column_names_rot = 45,
                          width = unit(12, "mm"))


colnames(responsetowounding) = c("SBC Training", "FCC Longitudinal", "FCC Validation",
                                 "FCC Subthreshold",  "SBC Testing")
colnames(responsetowounding)




#         working code for heatmap 
#=============================================================================
#
# --- number of down and up regulated proteins for each cohort ad top annotaiton
#

deps_long = read.csv("01_001_hptsd_somalogic_DEPs_long.csv", row.names = 1)
names(deps_long)


# exclude the civilian and the female cohorts


deps_long_significant = deps_long[deps_long$pvalue<0.05, ]

deps_long_significant_up = deps_long_significant[deps_long_significant$logFC>0,]
deps_long_significant_down = deps_long_significant[deps_long_significant$logFC<0,]


cohort_logfc_down = xtabs(~ cohort, data = deps_long_significant_down)
#View(cohort_logfc_down)

cohort_logfc_up = xtabs(~ cohort, data = deps_long_significant_up)
#View(cohort_logfc_up)

deps_down_up_significant = cbind(cohort_logfc_down, cohort_logfc_up)
#View(deps_down_up_significant)

deps_down_up_matrix = as.matrix(deps_down_up_significant)
View(deps_down_up_matrix)
deps_nocivilianfemales_matrix = deps_down_up_matrix[-c(2,3,7),]
View(deps_nocivilianfemales_matrix)
rownames(deps_nocivilianfemales_matrix) = c( "SBC Training", "FCC Longitudinal", "FCC Validation", "FCC Subthreshold", "SBC Testing")

deps_nocivilianfemales_matrix_orders = deps_nocivilianfemales_matrix[c(1,5,3,2,4),]

ha_deps = HeatmapAnnotation(.. = anno_barplot(deps_nocivilianfemales_matrix_orders, name = "total\nDEPs",
                                              show_annotation_name = F, annotation_name_side = "left", gp = gpar(fill = c("slateblue1", "tomato1")), gap = unit(8, "mm"), height = unit(1, "cm")),
                            pathways = anno_empty(border = F, height = unit(3, "mm")))


# lgd_protein_direct_of_regulat = Legend(at = 1:2, labels = c("up", "down"), title = "protein\nregulation \ndirection", legend_gp = gpar(fill = c("tomato1", "slateblue1")),
#                                        grid_height = unit(0.5, "cm"), grid_width = unit(6, "mm"))
#=========================================================================================
#
#------------ ggplot2 plot for consistent pathways as bottom annotation ----------------
#
setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/05252020_summaryAndConclusionGraphs")

# pathways_allcohorts = read.csv("01_hPTSD_SomaScan_meta_pathways_processes_and_mf_cc_across_allCohorts.csv")
# dim(pathways_allcohorts)
# bp_pathways = pathways_allcohorts[pathways_allcohorts$Pathway_Category=="bp",]
# levels(factor(bp_pathways$Pathway))
# 
# wound_inflam_coag_pathways = bp_pathways[(bp_pathways$Pathway == "Wound healing" | bp_pathways$Pathway == "Response to wounding" |
#                                             bp_pathways$Pathway == "Inflammatory response" | bp_pathways$Pathway == "Hemostasis"),] 
# 
# 
# levels(factor(wound_inflam_coag_pathways$Pathway))
# names(wound_inflam_coag_pathways)
# 
# #exclude the civilians and females
# levels(factor(wound_inflam_coag_pathways$Cohort))
# wound_inflam_coag_pathways_noCiviliansFemales = wound_inflam_coag_pathways[!(wound_inflam_coag_pathways$Cohort == "Civilian" | wound_inflam_coag_pathways$Cohort == "Female_AD" | wound_inflam_coag_pathways$Cohort == "Female_VT"),]
# 
# write.csv(wound_inflam_coag_pathways_noCiviliansFemales, "wound_inflam_coag_pathways_noCiviliansFemales.csv", row.names = F)

# wound healing related consistent pathways
wound_inflam_coag_pathways_noCiviliansFemales = read.csv("wound_inflam_coag_pathways_noCiviliansFemales_v2.csv")
View(wound_inflam_coag_pathways_noCiviliansFemales)
# vasculature development related consistent pathways
#vascul_angiogenesis_heartdevelop_pathways  = read.csv("vascul_angiogenesis_heartdevelop_pathways.csv")


levels(factor(wound_inflam_coag_pathways_noCiviliansFemales$Pathway))
names(wound_inflam_coag_pathways_noCiviliansFemales)
levels(factor(wound_inflam_coag_pathways_noCiviliansFemales$Cohort))

# excluding civilian and female cohorts
#vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales = vascul_angiogenesis_heartdevelop_pathways[!(vascul_angiogenesis_heartdevelop_pathways$Cohort=="Civilian"|vascul_angiogenesis_heartdevelop_pathways$Cohort=="Female_AD"|vascul_angiogenesis_heartdevelop_pathways$Cohort=="Female_VT"),]

# combining wound healing and vasculature development related pathways
#woundhealing_vasculaturedevelopment_nofemalenocivilian = rbind(wound_inflam_coag_pathways_noCiviliansFemales, vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales)
#View(woundhealing_vasculaturedevelopment_nofemalenocivilian)

levels(factor(wound_inflam_coag_pathways_noCiviliansFemales$Cohort))

pal <- wes_palette("Zissou1", 100, type = "continuous")
#pal2 = colorRamp2(c(0, 2, 4),c("yellow", "orange", "red"))
wound_inflam_coag_pathways_noCiviliansFemales$Pathway = factor(wound_inflam_coag_pathways_noCiviliansFemales$Pathway, levels=c("Response to wounding","Wound healing", "Inflammatory response", "Hemostasis"))


wound_inflam_coag_pathways_noCiviliansFemales$Cohort <- factor(wound_inflam_coag_pathways_noCiviliansFemales$Cohort, levels=c("SBC Training","SBC Testing", "FCC Validation", "FCC Longitudinal",  "FCC Subthreshold"))
ggplot() + geom_point(data = wound_inflam_coag_pathways_noCiviliansFemales, aes(x = Cohort, y = Pathway, size = numbProteins, color = log(-log(FDR,10),2))) +
  scale_colour_gradient2(low = "goldenrod1", mid = "gold1", high = "darkorange1", midpoint = 1, space = "Lab",  guide = "colourbar") +
  theme(text=element_text(family="mono", size = 16)) + theme_bw() +
  labs(size = "number of\n significant proteins\n(cases vs controls)\n enriched to the \ncorresponding pathway" ,y= "pathway or localization ", color = "FDRcorr. enrichment\n log2(-log10(q-value))") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

#subgroup = sample(letters[1], 12, replace = TRUE, prob = 1)
subgroup = length(colnames(responsetowounding))
panel_fun_ggplot2 = function(index, nm) {
  wound_inflam_coag_pathways_noCiviliansFemales$Pathway = factor(wound_inflam_coag_pathways_noCiviliansFemales$Pathway, levels=c("Response to wounding","Wound healing", "Inflammatory response", "Hemostasis"))
  wound_inflam_coag_pathways_noCiviliansFemales$Cohort <- factor(wound_inflam_coag_pathways_noCiviliansFemales$Cohort, levels=c("SBC Training","SBC Testing", "FCC Validation", "FCC Longitudinal",  "FCC Subthreshold"))
  
  g = ggplot(data = wound_inflam_coag_pathways_noCiviliansFemales, aes(x = Cohort, y = Pathway, size = numbProteins, color = log(-log(FDR,10),2))) + 
    geom_point() + theme_bw() + 
    scale_colour_gradient2(low = "goldenrod1", mid = "gold1", high = "darkorange1", midpoint = 1, space = "Lab",  guide = "colourbar") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    theme(panel.border = element_blank(),
          panel.background = element_blank()) 
  
  g = grid.grabExpr(print(g))
  pushViewport(viewport())
  grid.rect()
  grid.draw(g)
  popViewport()
}

ggplot_lgd_col = colorRamp2(c(0, 2, 4), c("goldenrod1", "gold1", "darkorange1"))
#ggplot_lgd = list(title = "log2(-log10(FDR))", col = ggplot_lgd_col, at = c(0, 1, 2, 3, 4), 
#                    labels = c("0", "1", "2", "3", "4"))

ggplot_anno = anno_zoom(align_to = subgroup, which = "column", panel_fun = panel_fun_ggplot2, height = unit(3, "cm"), 
                        extend = unit(0, "mm"), gap = unit(4, "mm"))





##=================================================================
#
#--------- miRs as left annotation
#

mirprotein = read.csv("wounding inflammation coagulation mir-protein good_v2.csv", row.names = 1)
#mirprotein = mirprotein[-c(35,37,69),]
# subsetting the data based on the common proteins between the two sets of heatmaps
#mirprotein =  mirprotein[!(row.names(mirprotein) %in% common16proteins_toberemoved),]

View(mirprotein)


pvalue = mirprotein$P.Value_mir 

fold_col_fun = colorRamp2(c(-2, 0, 2), c("cyan", "white", "orange")) 
ha_miR = rowAnnotation(microRNA = anno_text(mirprotein$miR_ID,gp = gpar(fontsize = 8, fontface = "bold"), just = "right", location = unit(1, "npc")),
                       .= anno_simple(mirprotein$logFC_mir,  col = fold_col_fun, na_col = "white", pch = ifelse((pvalue<0.05 & !is.na(pvalue)),"*","")))


View(responsetowounding)
colnames(responsetowounding)

responsetowounding_ordered = responsetowounding[,c(1,5,3,2,4)]
View(responsetowounding_ordered)

col_fun = colorRamp2(c(-0.4, 0, 0.35), c("blue", "white", "red"))
htmp_responsetowounding = Heatmap(responsetowounding_ordered, col = col_fun,
                                  show_heatmap_legend = F,
                                  top_annotation = ha_deps,
                                  column_title = "",#"response to wounding",
                                  column_title_side = c("top"),
                                  #column_title_gp = gpar(fontsize = 7, fontface ="bold"),
                                  column_title_rot = 0,
                                  bottom_annotation = HeatmapAnnotation(pathwayBP = ggplot_anno),
                                  left_annotation = ha_miR,
                                  show_row_names = T,
                                  row_names_gp = gpar(fontsize = 8, fontface ="bold"),
                                  #right_annotation = ha_hemostasis,
                                  #right_annotation = htmp_dmr,
                                  #right_annotation = htmp_metabolite,
                                  #cluster_rows = F,
                                  cluster_columns = F,
                                  #right_annotation = ha3,
                                  show_row_dend = F,
                                  show_column_dend = F,
                                  show_column_names = T,
                                  #row_km = 2, column_km = 2,
                                  column_names_rot = 40,
                                  column_names_gp = gpar(fontsize = 8, fontface ="bold"),
                                  row_km = 2,
                                  row_title = NULL,
                                  rect_gp = gpar(col = "white", lwd = 2),
                                  
                                  layer_fun = function(j, i, x, y, w, h, fill) {
                                    # restore_matrix() is explained after this chunk of code
                                    ind_mat = restore_matrix(j, i, x, y)
                                    for(ir in seq_len(nrow(ind_mat))) {
                                      # start from the second column
                                      for(ic in seq_len(ncol(ind_mat))) {
                                        ind = ind_mat[ir, ic] # previous column
                                        v = responsetowounding_pv[i[ind], j[ind]]
                                        
                                        grid.points(x[ind], y[ind], 
                                                    pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(1, "mm"))
                                        
                                      }
                                    }
                                  }
)     



# draw the heatmaps and annotations
draw(htmp_responsetowounding + htmp_dmr + htmp_woundhealing + htmp_inflam + htmp_hemostasis + htmp_metabolite,  auto_adjust = FALSE)#, heatmap_legend = packaged_legends, heatmap_legend_side = NULL,  adjust_annotation_extension = TRUE, 
#merge_legend = TRUE)



# ht_list_main = htmp_responsetowounding  %v%  htmp_vasculaturedevelopment
# htmp_dmr_list = htmp_dmr %v% htmp_va_dmr
# htmp_woundhealing_angiogenesis_list = htmp_woundhealing %v% htmp_angiogenesis
# htmp_inflam_heartdevelop_list = htmp_inflam %v% htmp_heartdevelop
# htmp_hemostasis
# htmp_metabolites_list = htmp_metabolite %v% htmp_va_metabolite

# decorate the pathwayBP for the ggplot_anno bottom annotation
decorate_annotation("pathwayBP", {
  grid.text("response to wounding", 1.05, 0.37, default.units = "npc", just = "left")
  grid.text("wound healing", 1.05, 0.53, default.units = "npc", just = "left")
  grid.text("inflammatory response", 1.05, 0.69, default.units = "npc", just = "left")
  grid.text("hemostasis (coagulation)", 1.05, 0.85, default.units = "npc", just = "left")
  
})


# lable the main heatmap, since all of the proteins shown in the heatmap are involved in the pathway
# "response to wounding"
decorate_annotation("pathways", { 
  grid.text("response to wounding", -0.5, 0.4, default.units = "npc", just = "left", gp = gpar(fontsize = 10, fontface="bold"))
  
})

# decorate the bar graph (describe what it is)
decorate_annotation("..", {
  grid.text("number of \nsignificant \nproteins", -0.4, 0.63, default.units = "npc", just = "right", gp = gpar(fontsize = 8, fontface= "bold", fontfamily="Helvetica"))
  #grid.rect(x = 0.6, y = 0.9, width = 0.5, height = 0.5, draw(lgd_protein_direct_of_regulat), default.units = "native")#, gp = gpar(fill = "tomato1"))# gp = gpar(fill = c("tomato1", "slateblue1")), 
  grid.rect(x = 0.85, y = 0.72,  gp = gpar(fill = "tomato1"),width = 0.07, height = 0.1,)
  grid.rect(x = 0.85, y = 0.6,  gp = gpar(fill = "slateblue1"),width = 0.07, height = 0.1,)
  grid.text("up",0.78,0.76,  default.units = "npc", just = "right", gp = gpar(fontsize = 6, fontface="bold"))
  grid.text("down",0.78,0.62,  default.units = "npc", just = "right", gp = gpar(fontsize = 6, fontface="bold") )
  #grid.text("regulation",0.6, 0.92, default.units = "npc", just = "right" )
})

# label the methylation heatmap
decorate_annotation("dmrprotein_anno", {
  grid.text("differentially methylated \nregions (SBC Training)", -1.9, 1.08, default.units = "npc", just = "left" , gp = gpar(fontsize = 11, fontface="bold"))
  
})

#decorate the heatmap orannotation of wound healing
decorate_annotation("coagulation", {
  grid.text("wound healing", -1.5, 1.08, default.units = "npc", just = "left" , gp = gpar(fontsize = 11, fontface="bold"))
  
})

#decorate/label the inflammatory heatmap and text annotation
decorate_annotation("coagulation", {
  grid.text("inflammatory \n response", 0.7, 1.08, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface="bold"))
  
})

#label the hemostasis heatmap or annoation
decorate_annotation("coagulation", {
  grid.text("hemostasis\n(coagulation)", 2.6, 1.08, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface="bold"))
  
})

# label the metabolites heatmap and text annotation
decorate_annotation("coagulation", {
  grid.text("significant metabolites (SBC Training)", 4.5, 1.08, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface="bold"))
  
})





decorate_annotation(".",{
  grid.text("* < 0.05", -2.16, -0.77, default.units = "npc", just = "left")
  grid.text("p-value", -3.1, -0.77, default.units = "npc", just = "right")
})





list_components()





#==============================================================================================
################################################################################################
####################################################################################################
#=============================================================================================================
#
# build the complex graph piece by piece for vasculature development angiogenesis and heart development
rm(list = ls())
#
# for all cohorts excluding civilians and females
#
#===============================================================================

# packages and libraries

library(seriation)
library(ComplexHeatmap)
library(circlize)
library(gridtext)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(GetoptLong)
library(methods)
library(reshape)


#========================================================
#
# ----  protein - metabolites interactome ----------
#
# -- reading the metabolite data along with the proein ids

rm(list =ls())
setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/05252020_summaryAndConclusionGraphs")

#==========================================================
#
#---------- metabolite for vasculature, angiogenesis and heart development ------
#

vas_proteinmetabolite = read.csv("vasculatureangiogenesisheartdevelopment_protein_metabolite_good_v4.csv", row.names = 1)
View(vas_proteinmetabolite)
vas_proteinmetabolite = vas_proteinmetabolite[-c(25,26),]
View(vas_proteinmetabolite)

#subset by common betwen the two groups of heatmaps
#vas_proteinmetabolite = vas_proteinmetabolite[vas_proteinmetabolite==common16proteins_toberemoved,]
#View(vas_proteinmetabolite)


# data frame to matrix
vas_protmetabol_matrix = as.matrix(vas_proteinmetabolite[,1:6])

# change NA to zeros
vas_protmetabol_matrix[is.na(vas_protmetabol_matrix)] <- 0

#View(vas_protmetabol_matrix)
# text annotation for the metabolite heatmap
ha_va_protmetabol = rowAnnotation(vas_protmetabol_anno = anno_text(vas_proteinmetabolite$metabolites, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun_metabol = colorRamp2(c(-0.5, 0, 0.4), c("dodgerblue1", "white", "deeppink1"))
htmp_va_metabolite = Heatmap(vas_protmetabol_matrix, 
                             show_row_names = F, 
                             column_title = "",#"significant\nmetabolites",
                             column_title_side = c("top"),
                             column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                             column_title_rot = 0,
                             show_heatmap_legend = F,
                             rect_gp = gpar(col = "white", lwd = 2),
                             cluster_rows = F,
                             cluster_columns = F, 
                             show_column_names = F,
                             right_annotation = ha_va_protmetabol, 
                             col = col_fun_metabol, column_names_rot = 45,
                             width = unit(15, "mm"))


#
# --------------- vasculature_angiogenesis dmr-proteins
#

va_dmrprotein = read.csv("vasculatureangiogenesisheartdevelopment_DMR_protein_good_v2.csv", row.names = 1)

va_dmrprotein = va_dmrprotein[-c(25,26),]

#View(va_dmrprotein)
names(va_dmrprotein)

va_dmrprotein_matrix = as.matrix(va_dmrprotein[,1:4])
#woundhealing_pv = wound_inflam_coag_proteinspathways[,22:24]
va_dmrprotein_matrix[is.na(va_dmrprotein_matrix)] <- 0
#colnames(woundhealing) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_va_dmrprotein = rowAnnotation(va_dmrprotein_anno = anno_text(va_dmrprotein$chr_TSS1500_TSS200, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun_dmr = colorRamp2(c(-0.015, 0, 0.01), c("green", "white", "magenta"))
htmp_va_dmr = Heatmap(va_dmrprotein_matrix, 
                      show_row_names = F, 
                      column_title ="", #"differentially methylated promoter regions",
                      column_title_side = c("top"),
                      column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                      column_title_rot = 0,
                      show_heatmap_legend = F,
                      rect_gp = gpar(col = "white", lwd = 2),
                      cluster_rows = F,
                      cluster_columns = F, 
                      show_column_names = F,
                      right_annotation = ha_va_dmrprotein, 
                      col = col_fun_dmr, column_names_rot = 45,
                      width = unit(10, "mm"))

#htmp_dmr_list = htmp_dmr %v% htmp_va_dmr


#======================================================================
#
# -- vasculature angiogenesis and heart developments protein -- pathways 
#

vascul_angiog_heart_proteinspathways = read.csv("vasculatureangiogenesisheartdevelopment_protein_pathway_good_v2.csv", row.names = 1)
#View(vascul_angiog_heart_proteinspathways)
vascul_angiog_heart_proteinspathways = vascul_angiog_heart_proteinspathways[-c(25:26),]

names(vascul_angiog_heart_proteinspathways)

# excluding civilians and females
vasculaturedevelopment = as.matrix(vascul_angiog_heart_proteinspathways[,c(1,4:6,8)])
vasculaturedevelopment_pv = vascul_angiog_heart_proteinspathways[,c(9,12:14,16)]

o1 = seriate(dist(vasculaturedevelopment), method = "GW")
o2 = seriate(dist(t(vasculaturedevelopment)), method = "GW")

angiogenesis = as.matrix(vascul_angiog_heart_proteinspathways[,c(19,21,20)])
colnames(angiogenesis)
angiogenesis[is.na(angiogenesis)] <- 0
colnames(angiogenesis) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_angiogenesis = rowAnnotation(angiogenesis = anno_text(vascul_angiog_heart_proteinspathways$angiogenesis, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp_angiogenesis = Heatmap(angiogenesis, 
                            show_row_names = F, 
                            column_title ="",#angiogenesis",
                            column_title_side = c("top"),
                            column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                            column_title_rot = 0,
                            show_heatmap_legend = F,
                            rect_gp = gpar(col = "white", lwd = 2),
                            cluster_rows = F,
                            cluster_columns = F, 
                            #show_column_names = F,
                            na_col = "white",
                            column_names_gp = gpar(fontsize = 7, fontface = "bold"),
                            right_annotation = ha_angiogenesis, 
                            col = col_fun, column_names_rot = 45,
                            width = unit(12, "mm"))

#
#----------- vasculogenesis -----------------------
#
vasculogenesis = as.matrix(vascul_angiog_heart_proteinspathways[,c(29,31,30)])
colnames(vasculogenesis)
#View(vasculogenesis)
vasculogenesis[is.na(vasculogenesis)] <- 0
colnames(vasculogenesis) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_vasculogenesis = rowAnnotation(vasculogenesis = anno_text(vascul_angiog_heart_proteinspathways$vasculogenesis, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp_vasculogenesis = Heatmap(vasculogenesis, 
                              show_row_names = F, 
                              column_title ="",#vasculogenesis",
                              column_title_side = c("top"),
                              column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                              column_title_rot = 0,
                              show_heatmap_legend = F,
                              rect_gp = gpar(col = "white", lwd = 2),
                              cluster_rows = F,
                              cluster_columns = F,
                              #show_column_names = F,
                              na_col = "white",
                              column_names_gp = gpar(fontsize = 7, fontface = "bold"),
                              right_annotation = ha_vasculogenesis, 
                              col = col_fun, column_names_rot = 45,
                              width = unit(12, "mm"))



#
# ----------- heart development
#
heartdevelopment = as.matrix(vascul_angiog_heart_proteinspathways[,c(24,26,25)])
colnames(heartdevelopment)
#View(heartdevelopment)
heartdevelopment[is.na(heartdevelopment)] <- 0
colnames(heartdevelopment) = c("SBC Training" , "SBC Testing",  "FCC Validation")

ha_heartdevelop = rowAnnotation(heartdevelopment = anno_text(vascul_angiog_heart_proteinspathways$heart_development, gp = gpar(fontsize = 8, fontface = "bold"), just = "left", location = unit(0, "npc")))

col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp_heartdevelop = Heatmap(heartdevelopment, 
                            show_row_names = F, 
                            column_title ="",#heart\ndevelopment",
                            column_title_side = c("top"),
                            column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                            column_title_rot = 0,
                            show_heatmap_legend = F,
                            rect_gp = gpar(col = "white", lwd = 2),
                            cluster_rows = F,
                            cluster_columns = F,
                            #show_column_names = F,
                            na_col = "white",
                            column_names_gp = gpar(fontsize = 7, fontface = "bold"),
                            right_annotation = ha_heartdevelop, 
                            col = col_fun, column_names_rot = 45,
                            width = unit(12, "mm"))



#         w
#=============================================================================
#
# --- number of down and up regulated proteins for each cohort ad top annotaiton
#

deps_long = read.csv("01_001_hptsd_somalogic_DEPs_long.csv", row.names = 1)
names(deps_long)


# exclude the civilian and the female cohorts


deps_long_significant = deps_long[deps_long$pvalue<0.05, ]

deps_long_significant_up = deps_long_significant[deps_long_significant$logFC>0,]
deps_long_significant_down = deps_long_significant[deps_long_significant$logFC<0,]


cohort_logfc_down = xtabs(~ cohort, data = deps_long_significant_down)
#View(cohort_logfc_down)

cohort_logfc_up = xtabs(~ cohort, data = deps_long_significant_up)
#View(cohort_logfc_up)

deps_down_up_significant = cbind(cohort_logfc_down, cohort_logfc_up)
#View(deps_down_up_significant)

deps_down_up_matrix = as.matrix(deps_down_up_significant)
#View(deps_down_up_matrix)
deps_nocivilianfemales_matrix = deps_down_up_matrix[-c(2,3,7),]
#View(deps_nocivilianfemales_matrix)
rownames(deps_nocivilianfemales_matrix) = c( "SBC Training", "FCC Longitudinal", "FCC Validation", "FCC Subthreshold", "SBC Testing")

deps_nocivilianfemales_matrix_orders = deps_nocivilianfemales_matrix[c(1,5,3,2,4),]

ha_deps = HeatmapAnnotation(.. = anno_barplot(deps_nocivilianfemales_matrix_orders, name = "total\nDEPs",
                                              show_annotation_name = F, annotation_name_side = "left", gp = gpar(fill = c("slateblue1", "tomato1")), gap = unit(8, "mm"), height = unit(1, "cm")),
                            pathways = anno_empty(border = F, height = unit(3, "mm")))


# lgd_protein_direct_of_regulat = Legend(at = 1:2, labels = c("up", "down"), title = "protein\nregulation \ndirection", legend_gp = gpar(fill = c("tomato1", "slateblue1")),
#                                        grid_height = unit(0.5, "cm"), grid_width = unit(6, "mm"))





#=========================================================================================
#
#------------ ggplot2 plot for consistent pathways as bottom annotation ----------------
#
setwd("~/Documents/   000-007March2019_SomaScanDataHPTSD/18Oct2019_SomaScan_ConfounderAssessed/03262020_SomaLogicPaperv2/05252020_summaryAndConclusionGraphs")

# pathways_allcohorts = read.csv("01_hPTSD_SomaScan_meta_pathways_processes_and_mf_cc_across_allCohorts.csv")
# dim(pathways_allcohorts)
# bp_pathways = pathways_allcohorts[pathways_allcohorts$Pathway_Category=="bp",]
# levels(factor(bp_pathways$Pathway))
# 
# wound_inflam_coag_pathways = bp_pathways[(bp_pathways$Pathway == "Wound healing" | bp_pathways$Pathway == "Response to wounding" |
#                                             bp_pathways$Pathway == "Inflammatory response" | bp_pathways$Pathway == "Hemostasis"),] 
# 
# 
# levels(factor(wound_inflam_coag_pathways$Pathway))
# names(wound_inflam_coag_pathways)
# 
# #exclude the civilians and females
# levels(factor(wound_inflam_coag_pathways$Cohort))
# wound_inflam_coag_pathways_noCiviliansFemales = wound_inflam_coag_pathways[!(wound_inflam_coag_pathways$Cohort == "Civilian" | wound_inflam_coag_pathways$Cohort == "Female_AD" | wound_inflam_coag_pathways$Cohort == "Female_VT"),]
# 
# write.csv(wound_inflam_coag_pathways_noCiviliansFemales, "wound_inflam_coag_pathways_noCiviliansFemales.csv", row.names = F)

# wound healing related consistent pathways
#wound_inflam_coag_pathways_noCiviliansFemales = read.csv("wound_inflam_coag_pathways_noCiviliansFemales.csv")
#View(wound_inflam_coag_pathways_noCiviliansFemales)
# vasculature development related consistent pathways
vascul_angiogenesis_heartdevelop_pathways  = read.csv("vascul_angiogenesis_heartdevelop_pathways_v2.csv")


levels(factor(vascul_angiogenesis_heartdevelop_pathways$Pathway))
names(vascul_angiogenesis_heartdevelop_pathways)
levels(factor(vascul_angiogenesis_heartdevelop_pathways$Cohort))

# excluding civilian and female cohorts
vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales = vascul_angiogenesis_heartdevelop_pathways[!(vascul_angiogenesis_heartdevelop_pathways$Cohort=="Civilian"|vascul_angiogenesis_heartdevelop_pathways$Cohort=="Female_AD"|vascul_angiogenesis_heartdevelop_pathways$Cohort=="Female_VT"),]

# combining wound healing and vasculature development related pathways
#woundhealing_vasculaturedevelopment_nofemalenocivilian = rbind(wound_inflam_coag_pathways_noCiviliansFemales, vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales)
#View(woundhealing_vasculaturedevelopment_nofemalenocivilian)

#levels(factor(woundhealing_vasculaturedevelopment_nofemalenocivilian$Cohort))

pal <- wes_palette("Zissou1", 100, type = "continuous")
#pal2 = colorRamp2(c(0, 2, 4),c("yellow", "orange", "red"))
vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales$Pathway = factor(vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales$Pathway, levels=c("Vasculature development", "Angiogenesis", "Heart development", "Vasculogenesis"))


vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales$Cohort <- factor(vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales$Cohort, levels=c("SBC Training","SBC Testing", "FCC Validation", "FCC Longitudinal",  "FCC Subthreshold"))
ggplot() + geom_point(data = vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales, aes(x = Cohort, y = Pathway, size = numbProteins, color = log(-log(FDR,10),2))) +
  scale_colour_gradient2(low = "goldenrod1", mid = "gold1", high = "darkorange1", midpoint = 1, space = "Lab",  guide = "colourbar") +
  theme(text=element_text(family="mono", size = 16)) + theme_bw() +
  labs(size = "number of\n significant proteins\n(cases vs controls)\n enriched to the \ncorresponding pathway" ,y= "pathway or localization ", color = "FDRcorr. enrichment\n log2(-log10(q-value))") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

cohorts_length = c("Training_VT", "Longitudinal_AD", "Testing_AD",
                   "Subthreshold_AD",  "Validating_VT")
#subgroup = sample(letters[1], 12, replace = TRUE, prob = 1)
subgroup = length(cohorts_length)
panel_fun_ggplot2 = function(index, nm) {
  vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales$Pathway = factor(vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales$Pathway, levels=c("Vasculature development", "Angiogenesis", "Heart development"))
  vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales$Cohort <- factor(vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales$Cohort, levels=c("SBC Training","SBC Testing", "FCC Validation", "FCC Longitudinal",  "FCC Subthreshold"))
  
  g = ggplot(data = vascul_angiogenesis_heartdevelop_pathways_noCivilianFemales, aes(x = Cohort, y = Pathway, size = numbProteins, color = log(-log(FDR,10),2))) + 
    geom_point() + theme_bw() + 
    scale_colour_gradient2(low = "goldenrod1", mid = "gold1", high = "darkorange1", midpoint = 1, space = "Lab",  guide = "colourbar") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    theme(panel.border = element_blank(),
          panel.background = element_blank()) 
  
  g = grid.grabExpr(print(g))
  pushViewport(viewport())
  grid.rect()
  grid.draw(g)
  popViewport()
}

ggplot_lgd_col = colorRamp2(c(0, 2, 4), c("goldenrod1", "gold1", "darkorange1"))
#ggplot_lgd = list(title = "log2(-log10(FDR))", col = ggplot_lgd_col, at = c(0, 1, 2, 3, 4), 
#                    labels = c("0", "1", "2", "3", "4"))

ggplot_anno = anno_zoom(align_to = subgroup, which = "column", panel_fun = panel_fun_ggplot2, height = unit(2.4, "cm"), 
                        extend = unit(0, "mm"), gap = unit(3, "mm"))



##==========================================================================================
#
#--------- vasculature development angiogenesis heart development miRs as left annotation
#

va_mirprotein = read.csv("vasculatureangiogenesisheartdevelopment_miR_protein_good_v2.csv", row.names = 1)

va_mirprotein = va_mirprotein[-c(25:26),]
#View(va_mirprotein)

va_pvalue = va_mirprotein$P.Value_mir 

fold_col_fun = colorRamp2(c(-2, 0, 2), c("cyan", "white", "orange")) 
ha_va_miR = rowAnnotation(va_microRNA = anno_text(va_mirprotein$miR_ID,gp = gpar(fontsize = 8, fontface = "bold"), just = "right", location = unit(1, "npc")),
                          .= anno_simple(va_mirprotein$logFC_mir,  col = fold_col_fun, na_col = "white", pch = ifelse((va_pvalue<0.05 & !is.na(va_pvalue)),"*","")))


colnames(vasculaturedevelopment) = c("SBC Training",  "FCC Longitudinal", "FCC Validation",
                                     "FCC Subthreshold",  "SBC Testing")
#View(vasculaturedevelopment)

vasculaturedevelopment_ordered =vasculaturedevelopment[,c(1,5,3,2,4)]
#View(vasculaturedevelopment_ordered)
# function for the main heatmap
col_fun = colorRamp2(c(-0.5, 0, 0.4), c("blue", "white", "red"))
htmp_vasculaturedevelopment = Heatmap(vasculaturedevelopment_ordered, name = "log2FC", col = col_fun,
                                      show_heatmap_legend = F,
                                      top_annotation = ha_deps,
                                      #column_title = ""#"vasculature development",
                                      #column_title_side = c("top"),
                                      #column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                      column_title_rot = 0,
                                      bottom_annotation = HeatmapAnnotation(pathwayBP = ggplot_anno),
                                      left_annotation = ha_va_miR,
                                      show_row_names = T,
                                      row_names_gp = gpar(fontsize = 8, fontface ="bold"),
                                      #right_annotation = ha_hemostasis,
                                      #cluster_rows = F,
                                      cluster_columns = F,
                                      #right_annotation = ha3,
                                      show_row_dend = F,
                                      show_column_dend = F,
                                      #row_km = 2, column_km = 2,
                                      column_names_rot = 40,
                                      column_names_gp = gpar(fontsize = 8, fontface ="bold"),
                                      row_km = 2,
                                      row_title = NULL,
                                      rect_gp = gpar(col = "white", lwd = 2),
                                      
                                      layer_fun = function(j, i, x, y, w, h, fill) {
                                        # restore_matrix() is explained after this chunk of code
                                        ind_mat = restore_matrix(j, i, x, y)
                                        for(ir in seq_len(nrow(ind_mat))) {
                                          # start from the second column
                                          for(ic in seq_len(ncol(ind_mat))) {
                                            ind = ind_mat[ir, ic] # previous column
                                            v = vasculaturedevelopment_pv[i[ind], j[ind]]
                                            
                                            grid.points(x[ind], y[ind], 
                                                        pch = 16, gp = gpar(col = ifelse(v < 0.05, "black", ifelse(v>=0.05 && v<0.1, "yellow", NA))), size = unit(1, "mm"))
                                            
                                          }
                                        }
                                      }
)     


#common15proteins = intersect(rownames(vasculaturedevelopment_ordered), rownames(responsetowounding_ordered))
#common16proteins_toberemoved = c(common15proteins,"A2M")

#View(common16proteins_toberemoved)

#===============================================================================================================================================

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun2, title = "protein: log2(fold change)", title_position = "leftcenter-rot",
              legend_height = unit(4.5, "cm"))

ggplot_pathway = Legend(title = "pathway\nenrichment\nq-value", col = ggplot_lgd_col, at = c(1, 2, 3, 4), 
                        labels = c("0.01", "1E^-5", "1E^-10", "1E^-15"))

#lgd_deps = Legend(title = "num.of.DEPs", pch = 1, type = "points", size = sort(unit((count_fdr$count)/5,"mm")), labels = sort(count_fdr$count))
lgd_significant = Legend(labels = c("p < 0.05", "0.05 <= p < 0.1 "), title = "protein\np-value", type = "points", title_position =  "topleft", 
                         grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
                         pch = 16, size = unit(5, "mm"), labels_gp = gpar(fontsize = 9), legend_gp = gpar(col = c(1,7)),  background = "white")
lgd_miR = Legend(title = "microRNA: log2(fold change)", col = fold_col_fun, at = c(-2, -1, 0, 1, 2), 
                 labels = c("-2", "-1", "0", "1", "2"), title_position = "leftcenter-rot",
                 legend_height = unit(5, "cm"))
lgd_dmr = Legend(title = "hypo-methylated  hyper-methylated", title_position = "leftcenter-rot",
                 legend_height = unit(6, "cm"), col = col_fun_dmr, at = c(-0.015,  0,  0.015), 
                 labels = c("","" ,""))
lgd_sig_mir = Legend(pch = "*", type = "points", labels = "< 0.05", title = "miR\np-value", title_position =  "topleft",
                     size = unit(5, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = c(1,7)),  background = "white")


lgd_count = Legend(labels = c(10,20,30,40,50), title = "pathway\n(protein count)", type = "points", title_position =  "topleft", 
                   grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
                   pch = 16, size = unit(2:6, "mm"), labels_gp = gpar(fontsize = 11), legend_gp = gpar(col = 1),  background = "white")

#package the legends (annotation and heatmap legends)
packaged_legends = packLegend(lgd_sig_mir, lgd_miR,  lgd1, lgd_significant, lgd_dmr, ggplot_pathway, lgd_count,  row_gap = unit(1, "cm"), direction = "horizontal")#, lgd_deps)
dev.off()
draw(packaged_legends)
# draw the heatmaps and annotations
draw(htmp_vasculaturedevelopment + htmp_va_dmr + htmp_angiogenesis + htmp_heartdevelop + htmp_vasculogenesis + htmp_va_metabolite,  auto_adjust = FALSE)#, heatmap_legend = packaged_legends, heatmap_legend_side = NULL,  adjust_annotation_extension = TRUE, 
#merge_legend = TRUE)

#==================================


# decorate the bar graph (describe what it is)
decorate_annotation("..", {
  grid.text("number of \nsignificant \nproteins", -0.4, 0.6, default.units = "npc", just = "right", gp = gpar(fontsize = 8, fontface= "bold", fontfamily="Helvetica"))
  #grid.rect(x = 0.6, y = 0.9, width = 0.5, height = 0.5, draw(lgd_protein_direct_of_regulat), default.units = "native")#, gp = gpar(fill = "tomato1"))# gp = gpar(fill = c("tomato1", "slateblue1")), 
  grid.rect(x = 0.85, y = 0.72,  gp = gpar(fill = "tomato1"),width = 0.07, height = 0.1,)
  grid.rect(x = 0.85, y = 0.6,  gp = gpar(fill = "slateblue1"),width = 0.07, height = 0.1,)
  grid.text("up",0.78,0.76,  default.units = "npc", just = "right", gp = gpar(fontsize = 6, fontface="bold"))
  grid.text("down",0.78,0.62,  default.units = "npc", just = "right", gp = gpar(fontsize = 6, fontface="bold") )
  #grid.text("regulation",0.6, 0.92, default.units = "npc", just = "right" )
})


#=======================================
# decorate the pathwayBP for the ggplot_anno bottom annotation

decorate_annotation("pathwayBP", {
  grid.text("vasculature developmment", 1.05, 0.40, default.units = "npc", just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  grid.text("angiogenesis", 1.05, 0.56, default.units = "npc", just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  grid.text("heart development", 1.05, 0.72, default.units = "npc", just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  grid.text("vasculogenesis", 1.05, 0.88, default.units = "npc", just = "left", gp = gpar(fontsize = 9, fontface = "bold"))
  
})


#=========================


# lable the main heatmap, since all of the proteins shown in the heatmap are involved in the pathway
# "vasculature development"
decorate_annotation("pathways", {
  grid.text("vasculature development", -0.6, 0.40, default.units = "npc", just = "left", gp = gpar(fontsize = 10, fontface="bold"))
  
})

#=====================================

# label the methylation heatmap
decorate_annotation("va_dmrprotein_anno", {
  grid.text("differentially methylated \nregions (SBC Training)", -0.2, 1.05, default.units = "npc", just = "left" , gp = gpar(fontsize = 11, fontface="bold"))
  
})


#==================================


#decorate the heatmap orannotation of angiogenesis
decorate_annotation("angiogenesis", {
  grid.text("angiogenesis", -1.8, 1.05, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface = "bold"))
  
})

#==============================

#decorate/label the heart development and text annotation
decorate_annotation("heartdevelopment", {
  grid.text("heart \ndevelopment", -1.4, 1.05, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface = "bold"))
  
})

#decorate/label the heart development and text annotation
decorate_annotation("vasculogenesis", {
  grid.text("vasculogenesis", -1.2, 1.05, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface = "bold"))
  
})

#===============================================


# label the metabolites heatmap and text annotation
decorate_annotation("vas_protmetabol_anno", {
  grid.text("significant metabolites\n(SBC Training)", -0.15, 1.05, default.units = "npc", just = "left", gp = gpar(fontsize = 11, fontface="bold"))
  
})


#vas_protmetabol_anno = anno_text(vas_proteinmetabolite

#====================================================================

# decorate_annotation(".",{
#   grid.text("* < 0.05", -2.18, 0.98, default.units = "npc", just = "left")
#   grid.text("p-value", -3.1, 0.98, default.units = "npc", just = "right")
# })
decorate_annotation(".",{
  grid.text("* < 0.05", -2.16, -0.77, default.units = "npc", just = "left")
  grid.text("p-value", -3.1, -0.77, default.units = "npc", just = "right")
})

decorate_annotation(".",{
  grid.text("* < 0.05", -2.16, -0.25, default.units = "npc", just = "left")
  grid.text("p-value", -3.1, -0.25, default.units = "npc", just = "right")
})


#==============================================================

list_components()



#package the legends (annotation and heatmap legends)
packaged_legends = packLegend(lgd1, lgd_significant, lgd_dmr, lgd_miR, lgd_count, ggplot_pathway, row_gap = unit(0.5, "cm"))#, direction = "horizontal")#, lgd_deps)


# heatmap lists for vertical stacking
ht_list_main = htmp_responsetowounding  %v%  htmp_vasculaturedevelopment
htmp_dmr_list = htmp_dmr %v% htmp_va_dmr
htmp_woundhealing_angiogenesis_list = htmp_woundhealing %v% htmp_angiogenesis
htmp_inflam_heartdevelop_list = htmp_inflam %v% htmp_heartdevelop
htmp_hemostasis
htmp_metabolites_list = htmp_metabolite %v% htmp_va_metabolite



draw(ht_list_main, ht_gap = unit(6,  "mm"))

# draw the heatmaps and annotations
draw(ht_list_main + htmp_dmr_list + htmp_woundhealing_angiogenesis_list + htmp_inflam_heartdevelop_list + htmp_metabolites_list,  auto_adjust = FALSE, heatmap_legend = packaged_legends, heatmap_legend_side = NULL,  adjust_annotation_extension = TRUE, 
     merge_legend = TRUE)









################################====================================


#decorate the heatmap orannotation of angiogenesis
decorate_annotation("angiogenesis", {
  grid.text("angiogenesis", -1.2, 1.05, default.units = "npc", just = "left")
  
})

#decorate/label the heart development and text annotation
decorate_annotation("heartdevelopment", {
  grid.text("heart \ndevelopment", -1, 1.06, default.units = "npc", just = "left")
  
})


# lable the main heatmap, since all of the proteins shown in the heatmap are involved in the pathway
# "vasculature development"
decorate_annotation("pathways", {
  grid.text("vasculature development", -0.1, 0.4, default.units = "npc", just = "left")
  
})













