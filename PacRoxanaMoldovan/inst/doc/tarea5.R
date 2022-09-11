## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
doAll = TRUE
dir_datos = "../data/"
centrado="center"

## ---- Paquetes, eval=doAll, warning=FALSE-------------------------------------
pacman::p_load(edgeR,SummarizedExperiment, ggplot2, pbapply, limma, 
               DESeq2, ReportingTools)

## ---- Datos1, eval=doAll, warning=FALSE---------------------------------------
datos="PRJNA821482.rda"
load(paste(dir_datos, datos, sep=""))

## ---- pData, eval=doAll, warning=FALSE----------------------------------------
# Datos fenotípicos de la primera muestra:
colData(PRJNA821482) 

## ---- edgeR_DGEList, eval=doAll, warning=FALSE--------------------------------
dge <- edgeR::DGEList(counts = assay(PRJNA821482), 
               group = colData(PRJNA821482)$Tratamiento)

## ---- edgeR_summary, eval=doAll, warning=FALSE--------------------------------
summary(dge$samples[,"lib.size"])

## ---- edgeR_porcent, eval=doAll, warning=FALSE--------------------------------
total = nrow(rowData(PRJNA821482))
# funcion que calcula el procentaje de genes no anotados:
porcent <- function(x, t){
  s = sum(is.na(x))
  paste0( s/t *100, "%")
}

paste0("Encontramos ", total, " genes, de estos no presentan anotación:")
s <- sapply(rowData(PRJNA821482), porcent, total)  
s

## ---- edgeR_cont, eval=doAll, warning=FALSE-----------------------------------
# Conteos por millon:  
c <- cpm(PRJNA821482)
mantener <- rowSums(c > 0,5) > 2
table(mantener)

## ---- edgeR_descart, eval=doAll, warning=FALSE--------------------------------

dge <- dge[mantener, keep.lib.sizes=FALSE]

## ---- edgeR_summary2, eval=doAll, warning=FALSE-------------------------------
summary(dge$samples[,"lib.size"])

## ---- edgeR_plot, eval=doAll, warning=FALSE, out.width = "80%", fig.asp = 0.8, fig.width = 8, fig.cap= "Fig1. Representación del tamaño de las librerias por muestra."----
ggplot2::ggplot(dge$samples, aes(y= lib.size, x =row.names(dge$samples),
                                 fill=lib.size)) + 
                geom_col() + coord_flip() + theme_minimal() +
                xlab("Muestra") +  ylab("Tamaño libreria")


## ---- edgeR_dispC, eval=doAll, warning=FALSE----------------------------------

dge.c = edgeR::estimateCommonDisp(dge)
et.c = exactTest(dge.c)

paste0("La dispersión obtenida es ", dge.c$common.dispersion)

## ---- edgeR_sig, eval=doAll, warning=FALSE------------------------------------
topTags(et.c, 2)

## ---- edgeR_sig2, eval=doAll, warning=FALSE-----------------------------------
sum(et.c$table[,"PValue"]< 0.05)

## ---- edgeR_Ajuste, eval=doAll, warning=FALSE---------------------------------
p_valores <- et.c$table[,"PValue"]
pv_ajust <- p.adjust(p_valores, "BH")

## ---- edgeR_Ajuste_sig, eval=doAll, warning=FALSE-----------------------------

sum(pv_ajust < 0.05)

## ---- edgeR_dispT, eval=doAll, warning=FALSE----------------------------------

dge.t = edgeR::estimateTagwiseDisp(dge.c)
et.t = exactTest(dge.t)

## ---- edgeR_sigT, eval=doAll, warning=FALSE-----------------------------------
topTags(et.t, 2)

## ---- edgeR_sig2T, eval=doAll, warning=FALSE----------------------------------
sum(et.t$table[,"PValue"]< 0.05)

## ---- edgeR_AjusteT, eval=doAll, warning=FALSE--------------------------------
p_valores_t <- et.t$table[,"PValue"]
pv_ajust_t <- p.adjust(p_valores_t, "BH")

## ---- edgeR_Ajuste_sigT, eval=doAll, warning=FALSE----------------------------

sum(pv_ajust_t < 0.05)

## ---- edgeR_Result, eval=doAll, warning=FALSE---------------------------------
# Combinamos los datos 
result <- cbind(rowData(PRJNA821482)[mantener, c("ENTREZID", "ENSEMBL")],
                Pval_Common = et.c$table[,"PValue"],
                adj_Pval_Common = pv_ajust, 
                Pval_Tagwise = et.t$table[,"PValue"],
                adj_Pval_Tagwise = pv_ajust_t)  

# Nos quedamos con los significativos:
sign <- (pv_ajust_t < 0.05 | pv_ajust < 0.05)  
Dat_inf <-result[sign,]

## ---- Parametros_Informe1, eval=doAll, warning=FALSE--------------------------
final_data = Dat_inf

entrezid = final_data$ENTREZID
ENTREZID= ifelse(is.na(entrezid), NA,
                 paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=",
                  entrezid,"'>",entrezid,"</a>"))


ensembl = final_data$ENSEMBL
ENSAMBL = ifelse(is.na(ensembl), NA,
                      paste0("<a href='http://www.ensembl.org/Macaca_mulatta",
                             "/Gene/Summary?db=core;g=",
                             ensembl ,"'>",ensembl,"</a>"))


## ---- Informe1, eval=doAll, warning=FALSE-------------------------------------
info = cbind(EntrezID = ENTREZID, Ensembl = ENSAMBL, 
             final_data[,3:6],  stringsAsFactors = F )

f_out = "T5_edgeR_clasico.DE"
report_dir = "Reports/"
titulo = paste0("PRJNA821482"," : Analisis de expresion ",
                "diferencial (edgeR clasico).")

report1 = ReportingTools::HTMLReport(shortName = f_out, title = titulo,
                                     reportDirectory = report_dir)
publish(info,report1)
finish(report1)

## ---- edgeR2, eval=doAll, warning=FALSE---------------------------------------
which(is.na(colData(PRJNA821482)$TipoCelular))
which(is.na(colData(PRJNA821482)$Tratamiento))


## ---- edgeR2_DGEList, eval=doAll, warning=FALSE-------------------------------
dge2 <- edgeR::DGEList(counts = assay(PRJNA821482))

## ---- edgeR2_descart, eval=doAll, warning=FALSE-------------------------------
# Nos quedamo con la que presnetaban conteos altos del apartado anterior:
dge2 <- dge2[mantener, keep.lib.sizes=FALSE]

## ---- edgeR2_Variables, eval=doAll, warning=FALSE-----------------------------
# Tomamos las variables a estudiar.
tratamiento <- factor(colData(PRJNA821482)$Tratamiento)
tipo_celular <- factor(colData(PRJNA821482)$TipoCelular)

table(tratamiento, tipo_celular)

## ---- edgeR2_comb_fac, eval=doAll, warning=FALSE------------------------------
# Creamos una función que concatene y asigne los nombres:
concat <- function(x, t=tratamiento, tc=tipo_celular){ #Numero de fila
  a = paste(t[x], tc[x],sep="_")
  a
}

# Creamos la variable factor:
fac = lapply(list(1:ncol(PRJNA821482)), FUN=concat)
tt = factor(unlist(fac))
tt

## ---- edgeR2_factores, eval=doAll, warning=FALSE------------------------------
# Determinamos la matriz de diseño:  
matriz <- model.matrix(~ 0 + tt)
colnames(matriz) <- levels(tt)
matriz


## ---- edgeR2_disp, eval=doAll, warning=FALSE----------------------------------
dge2 = edgeR::estimateDisp(dge2, matriz)

## ---- edgeR2_disp_info1, eval=doAll, warning=FALSE----------------------------
dge2$common.dispersion

## ---- edgeR2_disp_info2, eval=doAll, warning=FALSE----------------------------
summary(dge2$tagwise.dispersion)

## ---- edgeR2_disp_info3, eval=doAll, warning=FALSE----------------------------
summary(dge2$trended.dispersion)

## ---- edgeR2_BCV, eval=doAll, warning=FALSE, out.align=centrado,out.width = "80%", fig.asp = 0.8, fig.width = 5, fig.cap= "Fig2. Representación del coeficiente de variación biológico (BCV) frente a la abundancia de los genes."----
# Mostramos los resultados
edgeR::plotBCV(dge2)

## ---- edgeR2_ajust, eval=doAll, warning=FALSE---------------------------------
# Ajustamos los modelos lineales generalizados:  
fit <- edgeR::glmFit(dge2, design=matriz) 
# Influencia de las variables en los conteos:
lrt1 = edgeR::glmLRT(fit, coef=1:4)
topTags(lrt1)

## ---- edgeR2_contrast, eval=doAll, warning=FALSE------------------------------
contrast <- limma::makeContrasts( C_cel = (Normoxia_HSC + Hypoxia_HSC)
                                           - (Normoxia_K562 +Hypoxia_K562),
                                  C_HSC = Normoxia_HSC - Hypoxia_HSC,
                                  C_K562 = Normoxia_K562 - Hypoxia_K562,
                                  C_total = (Normoxia_HSC - Hypoxia_HSC) 
                                        - (Normoxia_K562 - Hypoxia_K562),
                                  levels = matriz)  


## ---- edgeR2_cel, eval=doAll, warning=FALSE-----------------------------------
lrt_cel = edgeR::glmLRT(fit, contrast = contrast[,"C_cel"]) 
p_ajust_cel = p.adjust(lrt_cel$table[,"PValue"])
sum(p_ajust_cel< 0.05)

## ---- edgeR2_cel_top, eval=doAll, warning=FALSE-------------------------------
# Los más significativos son:
topTags(lrt_cel)

## ---- edgeR2_HSC, eval=doAll, warning=FALSE-----------------------------------
lrt_HSC = edgeR::glmLRT(fit, contrast = contrast[,"C_HSC"]) 
p_ajust_HSC = p.adjust(lrt_HSC$table[,"PValue"])
sum(p_ajust_HSC< 0.05)

## ---- edgeR2_HSC_top, eval=doAll, warning=FALSE-------------------------------
# Los más significativos son:
topTags(lrt_HSC)

## ---- edgeR2_K562, eval=doAll, warning=FALSE----------------------------------
lrt_K562 = edgeR::glmLRT(fit, contrast = contrast[,"C_K562"]) 
p_ajust_K562 = p.adjust(lrt_K562$table[,"PValue"])
sum(p_ajust_K562< 0.05)

## ---- edgeR2_K562top, eval=doAll, warning=FALSE-------------------------------
# Los más significativos son:
topTags(lrt_K562)

## ---- edgeR2_total, eval=doAll, warning=FALSE---------------------------------
lrt_total = edgeR::glmLRT(fit, contrast = contrast[,"C_total"]) 
p_ajust_total = p.adjust(lrt_total$table[,"PValue"])
sum(p_ajust_total< 0.05)

## ---- edgeR2_total_top, eval=doAll, warning=FALSE-----------------------------
# Los más significativos son:
topTags(lrt_total)

## ---- edgeR2_Result, eval=doAll, warning=FALSE--------------------------------
# Combinamos los datos 
result2 <- cbind(rowData(PRJNA821482)[mantener, c("ENTREZID", "ENSEMBL")],
                lrt_total$table, adj_Pval= p_ajust_total)  

# Nos quedamos con los significativos:
sign <- (p_ajust_total < 0.05)  
Dat_inf2 <-result2[sign,]

## ---- Parametros_Informe2, eval=doAll, warning=FALSE--------------------------
final_data = Dat_inf2

entrezid = final_data$ENTREZID
ENTREZID= ifelse(is.na(entrezid), NA,
                 paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=",
                  entrezid,"'>",entrezid,"</a>"))


ensembl = final_data$ENSEMBL
ENSAMBL = ifelse(is.na(ensembl), NA,
                      paste0("<a href='http://www.ensembl.org/Macaca_mulatta",
                             "/Gene/Summary?db=core;g=",
                             ensembl ,"'>",ensembl,"</a>"))


## ---- Informe2, eval=doAll, warning=FALSE-------------------------------------
info = cbind(EntrezID = ENTREZID, Ensembl = ENSAMBL, 
             final_data[,3:7],  stringsAsFactors = F )

f_out = "T5_edgeR_contraste.DE"
report_dir = "Reports/"
titulo = paste0("PRJNA821482"," : Analisis de expresion diferencial",
                " en hipoxia en ambas lineas celulares.")

report2 = ReportingTools::HTMLReport(shortName = f_out, title = titulo,
                                     reportDirectory = report_dir)
publish(info,report2)
finish(report2)

## ---- edgeR3_ajust, eval=doAll, warning=FALSE---------------------------------
# Ajustamos los modelos lineales por cuasiverosimilitud.:  
qlfit <- edgeR::glmQLFit(dge2, design=matriz)  

# Influencia de las variables en los conteos:
Ftest1 <- edgeR::glmLRT(qlfit, coef=1:4)
topTags(Ftest1)

## ---- edgeR3_cel, eval=doAll, warning=FALSE-----------------------------------
ql_lrt_cel = edgeR::glmLRT(qlfit, contrast = contrast[,"C_cel"]) 
p_ajust_cel = p.adjust(ql_lrt_cel$table[,"PValue"])
sum(p_ajust_cel< 0.05)

## ---- edgeR3_cel_top, eval=doAll, warning=FALSE-------------------------------
# Los más significativos son:
topTags(ql_lrt_cel)

## ---- edgeR3_HSC, eval=doAll, warning=FALSE-----------------------------------
ql_lrt_HSC = edgeR::glmLRT(qlfit, contrast = contrast[,"C_HSC"]) 
p_ajust_HSC = p.adjust(ql_lrt_HSC$table[,"PValue"])
sum(p_ajust_HSC< 0.05)

## ---- edgeR3_HSC_top, eval=doAll, warning=FALSE-------------------------------
# Los más significativos son:
topTags(ql_lrt_HSC)

## ---- edgeR3_K562, eval=doAll, warning=FALSE----------------------------------
ql_lrt_K562 = edgeR::glmLRT(qlfit, contrast = contrast[,"C_K562"]) 
p_ajust_K562 = p.adjust(ql_lrt_K562$table[,"PValue"])
sum(p_ajust_K562< 0.05)

## ---- edgeR3_K562top, eval=doAll, warning=FALSE-------------------------------
# Los más significativos son:
topTags(ql_lrt_K562)

## ---- edgeR3_total, eval=doAll, warning=FALSE---------------------------------
ql_lrt_total = edgeR::glmLRT(qlfit, contrast = contrast[,"C_total"]) 
p_ajust_total = p.adjust(ql_lrt_total$table[,"PValue"])
sum(p_ajust_total< 0.05)

## ---- edgeR3_total_top, eval=doAll, warning=FALSE-----------------------------
# Los más significativos son:
topTags(ql_lrt_total)

## ---- edgeR3_Result, eval=doAll, warning=FALSE--------------------------------
# Combinamos los datos 
result3 <- cbind(rowData(PRJNA821482)[mantener, c("ENTREZID", "ENSEMBL")],
                lrt_total$table, adj_Pval= p_ajust_total)  

# Nos quedamos con los significativos:
sign <- (p_ajust_total < 0.05)  
Dat_inf3 <-result2[sign,]

## ---- Parametros_Informe3, eval=doAll, warning=FALSE--------------------------
final_data = Dat_inf3

entrezid = final_data$ENTREZID
ENTREZID= ifelse(is.na(entrezid), NA,
                 paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=",
                  entrezid,"'>",entrezid,"</a>"))


ensembl = final_data$ENSEMBL
ENSAMBL = ifelse(is.na(ensembl), NA,
                      paste0("<a href='http://www.ensembl.org/Macaca_mulatta",
                             "/Gene/Summary?db=core;g=",
                             ensembl ,"'>",ensembl,"</a>"))


## ---- Informe3, eval=doAll, warning=FALSE-------------------------------------
info = cbind(EntrezID = ENTREZID, Ensembl = ENSAMBL, 
             final_data[,3:7],  stringsAsFactors = F )

f_out = "T5_edgeR_cuasiverosimilitud.DE"
report_dir = "Reports/"
titulo = paste0("PRJNA821482"," : Analisis de expresion diferencial en hipoxia",
                " en ambas lineas celulares mediante cuasiverosimilitud.")

report3 = ReportingTools::HTMLReport(shortName = f_out, title = titulo,
                                     reportDirectory = report_dir)
publish(info,report3)
finish(report3)

## ---- DESeq2_mat, eval=doAll, warning=FALSE-----------------------------------
matriz

## ---- DESeq2_1, eval=doAll, warning=FALSE-------------------------------------
dds <- DESeq2::DESeqDataSet(PRJNA821482, design = (matriz))
dds

## ---- DESeq2_2, eval=doAll, warning=FALSE-------------------------------------
mantener <- rowSums(counts(dds) > 0,5) > 2
dds = dds[mantener,]

## ---- DESeq2_dif, eval=doAll, warning=FALSE-----------------------------------
dds <- DESeq(dds)

## ---- DESeq2_resname, eval=doAll, warning=FALSE-------------------------------
resultsNames(dds)

## ---- DESeq2_res, eval=doAll, warning=FALSE-----------------------------------
# Contraste:
C_total <- limma::makeContrasts( ((Normoxia_HSC - Hypoxia_HSC) 
                                  - (Normoxia_K562 - Hypoxia_K562)),
                                  levels = matriz)  
res_dds = results(dds, contrast = C_total,
                  alpha = 0.05,  pAdjustMethod = "BH")

res_dds

table(res_dds$padj<0.05)

## ---- DESeq2_Result, eval=doAll, warning=FALSE--------------------------------
# Combinamos los datos 
result4 <- cbind(rowData(PRJNA821482)[mantener, c("ENTREZID", "ENSEMBL")],
                res_dds)  

# Nos quedamos con los significativos:
sign <- (p_ajust_total < 0.05)  
Dat_inf4 <-result2[sign,]

## ---- Parametros_Informe4, eval=doAll, warning=FALSE--------------------------
final_data = Dat_inf4

entrezid = final_data$ENTREZID
ENTREZID= ifelse(is.na(entrezid), NA,
                 paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=",
                  entrezid,"'>",entrezid,"</a>"))


ensembl = final_data$ENSEMBL
ENSAMBL = ifelse(is.na(ensembl), NA,
                      paste0("<a href='http://www.ensembl.org/Macaca_mulatta",
                             "/Gene/Summary?db=core;g=",
                             ensembl ,"'>",ensembl,"</a>"))


## ---- Informe4, eval=doAll, warning=FALSE-------------------------------------
info = cbind(EntrezID = ENTREZID, Ensembl = ENSAMBL, 
             final_data[,3:7],  stringsAsFactors = F )

f_out = "T5_DESeq2.DE"
report_dir = "Reports/"
titulo = paste0("PRJNA821482"," : Analisis de expresion diferencial en hipoxia",
                " en ambas lineas celulares mediante cuasiverosimilitud.")

report4 = ReportingTools::HTMLReport(shortName = f_out, title = titulo,
                                     reportDirectory = report_dir)
publish(info,report4)
finish(report4)

