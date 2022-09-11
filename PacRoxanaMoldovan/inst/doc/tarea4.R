## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
doAll = TRUE
done = FALSE
dir_datos = "../data/"

## ---- Paquetes, eval=doAll, warning=FALSE-------------------------------------
pacman::p_load(pbapply)
pacman::p_load(SummarizedExperiment,EnrichmentBrowser,GSEABase, tami, GSA)
pacman::p_load(ReportingTools)
pacman::p_load(pd.cyngene.1.0.st,org.Mmu.eg.db, GO.db)

## ---- Datos1, eval=doAll, warning=FALSE---------------------------------------
datos="geod75759_rma.rda"
load(paste(dir_datos, datos, sep=""))

## ----generacion_GeneSetCollection_Go, eval=done, warning=FALSE----------------
#  
#  # Tomamos la información de Go:
#  frame <- toTable(org.Mmu.egGO2EG)
#  # Creamos un data frame con el identificador Go, la evidencia y el id del gen:
#  goframeData <- data.frame(frame$go_id, frame$Evidence, frame$gene_id)
#  goFrame = AnnotationDbi::GOFrame(goframeData, organism = "Macaca mulatta")
#  goAllFrame = AnnotationDbi::GOAllFrame(goFrame)
#  mccGOgsc_all = GSEABase::GeneSetCollection(goAllFrame, setType = GOCollection())
#  
#  # Guardamos la colección Go con todos los genes de Macaca mulatta:
#  save(mccGOgsc_all, file = paste0(dir_datos,"mccGOgsc_all.rda"))

## ----Carga_Go1, eval=doAll, warning=FALSE-------------------------------------
# Cargamos la colección de Go:
load(paste0(dir_datos,"mccGOgsc_all.rda"))

## ----Eliminar_genes_Go, eval=doAll, warning=FALSE-----------------------------
# Tomamos los que tienen correspondencia con nuestros identificadores:

#Función subsettingGeneSet:
SelectSubGeneSet <- function(GS, G_ID){ 
  # Introducimos el GeneSet y los ENTREZID
  geneIds(GS) = geneIds(GS)[is.element( geneIds(GS), G_ID)]
  GS
}

# Grupos de genes Go con nuestros identificadores:
mccGOgsc_reducido <- pblapply(mccGOgsc_all, SelectSubGeneSet,
                       G_ID=fData(geod75759_rma)$ENTREZID)
mccGOgsc <- GeneSetCollection(mccGOgsc_reducido)


## ----Tamaño_coleccion_Go, eval=doAll, warning=FALSE---------------------------
paste0("Tenemos ", length(mccGOgsc_all), " grupos segun la clsificación GO.")
paste0("El primer grupo contiene ", length(geneIds(mccGOgsc_all[[1]])), " genes.")

# Una vez eliminados:
paste0("Tenemos ", length(mccGOgsc), " grupos segun la clsificación GO.")
paste0("El primer grupo contiene ", length(geneIds(mccGOgsc[[1]])), " genes.")


## ----No_genes_coleccion_Go, eval=doAll, warning=FALSE-------------------------

count_Go <- pbsapply(geneIds(mccGOgsc), length)
summary(count_Go)

a = mccGOgsc[(lapply(geneIds(mccGOgsc), length) <= 1)]
paste0(length(a), " grupos GO presenta 0 o un único gen.")


## ----ELiminar_grupos_pequenios_Go, eval=doAll, warning=FALSE------------------
# Tomamos los grupos con más de un gen
mccGOgsc = mccGOgsc[(lapply(geneIds(mccGOgsc), length) > 1)]
paste0("Al final nos quedamos con ", length(mccGOgsc), "grupos Go")
count_Go <- pbsapply(geneIds(mccGOgsc), length)
summary(count_Go)
length(mccGOgsc)

## ----Guardado_final_Go, eval=doAll, warning=FALSE-----------------------------
# Guardamos el conjunto:
save(mccGOgsc, file = paste0(dir_datos,"mccGOgsc.rda"))

## ----generacion_GeneSetCollection_KEGG, eval=done, warning=FALSE--------------
#  # Obtenemos la coleccion de genes de KEGG:
#  mccKEGGgsc_all <- EnrichmentBrowser::getGenesets(org = "mcc", db ="kegg",
#                                                  gene.id.type = "ENTREZID",
#                                              return.type = "GeneSetCollection")
#  mccKEGGgsc_all
#  
#  
#  # Guardamos los datos:
#  save(mccKEGGgsc_all, file = (paste0(dir_datos, "mccKEGGgsc_all.rda")))

## ----Carga_KEGG1, eval=doAll, warning=FALSE-----------------------------------
# Cargamos los datos:
load(paste0(dir_datos, "mccKEGGgsc_all.rda"))

## ----Eliminar_genes_KEGG, eval=doAll, warning=FALSE---------------------------
#Función subsettingGeneSet
SelectSubGeneSet <- function(GS, G_ID){ 
  # Introducimos el GeneSet y los ENTREZID
  geneIds(GS) = geneIds(GS)[is.element( geneIds(GS), G_ID)]
  GS
}

# Grupos de genes KEGG con nuestros identificadores:
mccKEGGgsc_reducido <- pblapply(mccKEGGgsc_all, SelectSubGeneSet,
                       G_ID=fData(geod75759_rma)$ENTREZID)
mccKEGGgsc <- GeneSetCollection(mccKEGGgsc_reducido)
# Guardamos los datos:
save(mccKEGGgsc, file = (paste0(dir_datos, "mccKEGGgsc.rda")))

## ----Tamaño_coleccion_KEGG, eval=doAll, warning=FALSE-------------------------
paste0("Tenemos ", length(mccKEGGgsc_all), " grupos segun la clsificación KEGG.")
paste0("El primer grupo contiene ", length(geneIds(mccKEGGgsc_all[[1]])), " genes.")

# Una vez eliminados:
paste0("Tenemos ", length(mccKEGGgsc), " grupos segun la clsificación KEGG.")
paste0("El primer grupo contiene ", length(geneIds(mccKEGGgsc[[1]])), " genes.")


## ----Ej_grupos_KEGG, eval=doAll, warning=FALSE--------------------------------
head(names(mccKEGGgsc))

## ----No_genes_coleccion_KEGG, eval=doAll, warning=FALSE-----------------------

count_Go <- pbsapply(geneIds(mccKEGGgsc), length)
summary(count_Go)

a = mccKEGGgsc[(lapply(geneIds(mccKEGGgsc), length) <= 1)]
paste0(length(a), " grupos KEGG presenta menos de 3 genes.")


## ----Guardado_final_KEGG, eval=doAll, warning=FALSE---------------------------
# Guardamos los datos:
save(mccKEGGgsc, file = (paste0(dir_datos, "mccKEGGgsc.rda")))
load(paste0(dir_datos, "mccKEGGgsc.rda"))

## ---- clase_ExpressionSet, eval=doAll, warning=FALSE--------------------------
class(geod75759_rma)

## ---- Conversion_SummarizedExperiment1, eval=doAll, warning=FALSE-------------
se75759 <- makeSummarizedExperimentFromExpressionSet(geod75759_rma)
se75759

## ---- Conversion_SummarizedExperiment2, eval=doAll, warning=FALSE-------------
head(rowData(se75759), 2)

## ---- ExpressionSet1, eval=doAll, warning=FALSE-------------------------------
se75759 = EnrichmentBrowser::probe2gene(se75759, from = "PROBEID",
                                        to = "ENTREZID")
se75759
rowData(se75759)

## ---- Group, eval=doAll, warning=FALSE----------------------------------------
tratamiento <- colData(se75759)$Sample_treatment
tratamiento
grupos <- ifelse(tratamiento == "Formaldehyde", 0, 1)
colData(se75759) <- cbind(colData(se75759), GROUP=grupos)
head(colData(se75759),1)

## ----Limma, eval=doAll, warning=FALSE-----------------------------------------
se75759 = EnrichmentBrowser::deAna(expr = se75759, 
                                   padj.method = "BH", 
                                   de.method = "limma")
head(rowData(se75759),2 )

## ----Limma_sig, eval=doAll, warning=FALSE-------------------------------------
FDR = 0.1
sum(rowData(se75759)$ADJ.PVAL < FDR )

## ----Fisher_KEGG, eval=doAll, warning=FALSE-----------------------------------
se75759.oraKEGG <- EnrichmentBrowser::sbea(method = "ora", se = se75759,
                                           padj.method = "BH",
                                           gs=mccKEGGgsc, perm = 0,
                                           alpha = 0.05)
gsRanking(se75759.oraKEGG, signif.only = F)

## ----Fisher_GO, eval=doAll, warning=FALSE-------------------------------------
se75759.oraGO <- EnrichmentBrowser::sbea(method = "ora", se = se75759, 
                                         padj.method = "BH",
                                         gs=mccGOgsc, perm = 0,
                                         alpha = 0.05)
gsRanking(se75759.oraGO, signif.only = F)

## ----GSA_KEGG, results="hide", eval=doAll, warning=FALSE, message=FALSE-------
se75759.gsaKEGG <- EnrichmentBrowser::sbea(method = "gsa", se = se75759,
                                           padj.method = "BH",
                                           gs=mccKEGGgsc,alpha = 0.05)   
gsRanking(se75759.gsaKEGG, signif.only = F)

## ----GSA_KEGG_sig, eval=doAll, warning=FALSE----------------------------------
se75759.gsaKEGG$nr.sigs

## ----GSA_KEGG_informe, eval=doAll, warning=FALSE------------------------------
# Tomamos los resultados:
tabla = se75759.gsaKEGG$res.tbl

# Creamos la url para cada identificador:
ID = tabla$GENE.SET
KEGGID = ifelse(is.na(ID), NA, paste0("<a href='https://www.genome.jp/",
                                      "dbget-bin/www_bget?", ID, "'>",
                                      ID, "</a>"))

final_info_KEGG = data.frame(Grupo_KEGG = KEGGID, Score = tabla$SCORE,
                             P_val = tabla$PVAL, p_ajust=tabla$ADJ.PVAL  )  


## ----GSA_KEGG_informe2, eval=doAll, warning=FALSE-----------------------------
info = final_info_KEGG

f_out = "T4_GSA_KEGG.DE"
report_dir = "Reports/"
titulo = "geod75759 : Resultados GSA de grupos KEGG."

report1 = ReportingTools::HTMLReport(shortName = f_out, 
                                     title = titulo,
                                     reportDirectory = report_dir)
publish(info,report1)
finish(report1)


## ----GSA_tami_competitive_GO, eval=doAll, warning=FALSE-----------------------
geod75759_GO_comp <- GeneSetTest(x = geod75759_rma,
                              y=factor(pData(geod75759_rma)$Sample_treatment),
                              test = rowtmod, association="statistic",
                              correction="BH", GeneNullDistr = "randomization",
                              GeneSetNullDistr = "competitive",
                              alternative="two-sided", nmax = 100,
                              id = "ENTREZID", gsc=geneIds(mccGOgsc),
                              descriptive=mean, foutput = "geod75759_rma_comp")
geod75759_GO_comp_df = tidy(geod75759_GO_comp)

sig_GO_comp = geod75759_GO_comp_df$adjp<0.05

## ----GSA_tami_competitive_GO_sig, eval=doAll, warning=FALSE-------------------
sum(sig_GO_comp)

head(geod75759_GO_comp_df[sig_GO_comp,])

## ----GSA_comp_GO_informe, eval=doAll, warning=FALSE---------------------------
# Tomamos los resultados:
tabla2 = geod75759_GO_comp_df[sig_GO_comp,]

# Creamos la url para cada identificador:
ID = tabla2$GO
GOID = ifelse(is.na(ID), NA, tami::go2url(ID))

final_comp_GO = data.frame(Grupo_GO = GOID, satistic = tabla2$statistic,
                             p_value = tabla2$rawp, p_ajust=tabla2$adjp)  


## ----GSA_comp_GO_informe2, eval=doAll, warning=FALSE--------------------------
info = final_comp_GO

f_out = "T4_GSA_comp_GO.DE"
report_dir = "Reports/"
titulo = "geod75759 : Resultados GSA de grupos Go (hipotesis competitiva). "

reportGO1 = ReportingTools::HTMLReport(shortName = f_out, 
                                       title = titulo,
                                       reportDirectory = report_dir)
publish(info,reportGO1)
finish(reportGO1)


## ----GSA_tami_self-contained_GO, eval=doAll, warning=FALSE--------------------
geod75759_GO_self <- GeneSetTest(x = geod75759_rma,
                              y=factor(pData(geod75759_rma)$Sample_treatment),
                              test = rowtmod, association="statistic",
                              correction="BH", GeneNullDistr = "randomization",
                              GeneSetNullDistr = "self-contained",
                              alternative="two-sided", nmax = 100,
                              id = "ENTREZID", gsc=geneIds(mccGOgsc),
                              descriptive=mean, foutput = "geod75759_rma_self")
geod75759_GO_self_df = tidy(geod75759_GO_self)

sig_GO_self = geod75759_GO_self_df$adjp<0.05

## ----GSA_tami_contained_GO_sig, eval=doAll, warning=FALSE---------------------
sum(sig_GO_self)

head(geod75759_GO_self_df[sig_GO_self,])

## ----GSA_self_GO_informe, eval=doAll, warning=FALSE---------------------------
# Tomamos los resultados:
tabla3 = geod75759_GO_self_df[sig_GO_self,]

# Creamos la url para cada identificador:
ID = tabla3$GO
GOID = ifelse(is.na(ID), NA, tami::go2url(ID))

final_self_GO = data.frame(Grupo_GO = GOID, satistic = tabla3$statistic,
                            p_value = tabla3$rawp, p_ajust=tabla3$adjp)  


## ----GSA_self_GO_informe2, eval=doAll, warning=FALSE--------------------------
info = final_self_GO

f_out = "T4_GSA_self_GO.DE"
report_dir = "Reports/"
titulo = "geod75759 : Resultados GSA de grupos Go (hipotesis autocomparativa). "

reportGO2 = ReportingTools::HTMLReport(shortName = f_out, 
                                      title = titulo,
                                      reportDirectory = report_dir)
publish(info,reportGO2)
finish(reportGO2)


