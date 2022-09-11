## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
doAll = TRUE
dir_datos = "../data/"

## ---- Paquetes, eval=doAll, warning=FALSE-------------------------------------
pacman::p_load(ArrayExpress, genefilter, limma)
pacman::p_load(ReportingTools)

## ---- Datos, eval=doAll, warning=FALSE----------------------------------------
datos="geod75759_rma.rda"
load(paste(dir_datos, datos, sep=""))

## ---- pData, eval=doAll, warning=FALSE----------------------------------------
# Datos fenotípicos de la primera muestra:
(pData(geod75759_rma)[1,])

## ---- Def_factor, eval=doAll, warning=FALSE-----------------------------------
# Tomamos el factor a estudiar
tratamiento <- factor(pData(geod75759_rma)$Sample_treatment)

## ---- T_test, eval=doAll, warning=FALSE---------------------------------------
# Test de la t:
tt <- genefilter::rowttests(geod75759_rma, tratamiento)
head(tt)


## ---- P_valor, eval=doAll, warning=FALSE--------------------------------------
# P.value:
alfa = 0.05  # Fijamos un alfa de 0,05
pv_original = tt$p.value
table( pv_original <= alfa)

## ---- Significativos, eval=doAll, warning=FALSE-------------------------------

# Mostramos los genes más significativos:
sig <- tt[order(tt$p.value),]
head(sig)


## ---- Bonferroni, eval=doAll, warning=FALSE-----------------------------------

# Ajustamos por el método de Bonferroni:
pv_Bf <- p.adjust(tt$p.value, method = "bonferroni")
# Mostramos los genes significativos:
table( pv_Bf < alfa)


## ---- Benjamini-Hochberg, eval=doAll, warning=FALSE---------------------------

# Ajustamos por el método de Benjamini-Hochberg:
pv_BH <- p.adjust(pv_original, method = "BH")
# Mostramos los genes significativos:
table( pv_BH < alfa)


## ---- Data_frame, eval=doAll, warning=FALSE-----------------------------------

#Creación del data frame:
df1 =cbind(fData(geod75759_rma), Pvalue = pv_original, 
           adj.Pval.Bf =pv_Bf, adj.Pval.BH=pv_BH)
head(df1)


## ---- Parametros_Informe, eval=doAll, warning=FALSE---------------------------
final_data = df1

ID = ifelse(is.na(final_data$PROBEID), NA, final_data$PROBEID)

entrezid = final_data$ENTREZID
ENTREZID= ifelse(is.na(entrezid), NA,
                 paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=",
                  entrezid,"'>",entrezid,"</a>"))

ensemblT = final_data$ENSEMBLTRANS
ENSEMBLTRANS = ifelse(is.na(ensemblT), NA,
                      paste0("<a href='http://www.ensembl.org/",
                             "Macaca_mulatta/Transcript/Idhistory?t=",
                             ensemblT ,"'>",ensemblT,"</a>"))

ensembl = final_data$ENSAMBL
ENSAMBL = ifelse(is.na(ensembl), NA,
                      paste0("<a href='http://www.ensembl.org/Macaca_mulatta",
                             "/Gene/Summary?db=core;g=",
                             ensembl ,"'>",ensembl,"</a>"))

GENSYMBOL = ifelse(is.na(final_data$GENSYMBOL), NA, final_data$GENSYMBOL)


## ---- Informe, eval=doAll, warning=FALSE--------------------------------------
info = data.frame(ID=as.character(final_data$PROBEID), EntrezID = ENTREZID, 
                  EnsemblT=ENSEMBLTRANS,
                  EnsemblGen = ENSAMBL, GenSymbol = GENSYMBOL,
                  Pvalue = final_data$Pvalue,
                  adj.Pval.bonferroni = as.numeric(final_data$adj.Pval.Bf),
                  adj.Pval.BH = as.numeric(final_data$adj.Pval.BH), 
                  stringsAsFactors = F )

f_out = "T2_Ttest.DE"
report_dir = "Reports/"
titulo = paste0("geod75759"," : Analisis de expresion ",
                "diferencial (t-test con varianzas iguales)")

report1 = ReportingTools::HTMLReport(shortName = f_out, title = titulo,
                                     reportDirectory = report_dir)
publish(info,report1)
finish(report1)

## ---- Matriz, eval=doAll, warning=FALSE---------------------------------------

tratamiento <- pData(geod75759_rma)$Sample_treatment

# Matriz de diseño:
matriz <- model.matrix(~ 0 + tratamiento)
colnames(matriz) <- c( "Formaldehido", "Control")
matriz


## ---- Modelos_lineales, eval=doAll, warning=FALSE-----------------------------
ajusteA <- lmFit(geod75759_rma, matriz)
matriz_contrastes <- makeContrasts( dif = Control - Formaldehido,
                                    levels = matriz)
matriz_contrastes

## ---- Estimacion, eval=doAll, warning=FALSE-----------------------------------

ajusteB <- contrasts.fit(ajusteA, matriz_contrastes)
ajusteC <- eBayes(ajusteB)

## ---- Significativo_limma, eval=doAll, warning=FALSE--------------------------

sig <- topTable(ajusteC, number = 10, coef = 1, adjust.method ="BH")
sig

## ---- Variables, eval=doAll, warning=FALSE------------------------------------
# Tomamos las variables a estudiar.
tratamiento <- factor(pData(geod75759_rma)$Sample_treatment)
tejido <- factor(pData(geod75759_rma)$Sample_type)

table(tratamiento, tejido)

## ---- Combinar_factores, eval=doAll, warning=FALSE----------------------------

# Creamos una función que concatene y asigne los nombres:
concat <- function(x, t=tratamiento, tj=tejido){ #Numero de fila
  a = paste(t[x], tj[x],sep="_")
  a
}

# Creamos la variable factor:
fac = lapply(list(1:ncol(geod75759_rma)), FUN=concat)
trat_tej = factor(unlist(fac))
trat_tej


## ---- Matriz_2factores, eval=doAll, warning=FALSE-----------------------------

# Determinamos la matriz de diseño:
matriz2 <- model.matrix(~ 0 + trat_tej)
colnames(matriz2) <- levels(trat_tej)
matriz2

## ---- Modelos_lineales2, eval=doAll, warning=FALSE----------------------------
# Ajuste:
ajusteA2 <- lmFit(geod75759_rma, matriz2)

## ---- Contrastes2 , eval=doAll, warning=FALSE---------------------------------
#Construcción de los contrastes:
matriz_contrastes2 <- makeContrasts( dif1 = Filtered_Air_NE - Formaldehyde_NE ,
                                     dif2 = Filtered_Air_WBC- Formaldehyde_WBC,
                                     difT = (Filtered_Air_NE - Formaldehyde_NE)
                                          - (Filtered_Air_WBC- Formaldehyde_WBC),
                                     levels = matriz2)
matriz_contrastes2

## ---- Estimacion2, eval=doAll, warning=FALSE----------------------------------

ajusteB2 <- contrasts.fit(ajusteA2, matriz_contrastes2)
ajusteC2 <- eBayes(ajusteB2)

## ---- Sig1_limma2, eval=doAll, warning=FALSE----------------------------------
# Mostramos los genes con la expresión diferencial más significativa:
TT1 <-topTable(ajusteC2, coef = 1, adjust.method ="BH", number = Inf,)
head(TT1, 3)

## ---- Sig2_limma2, eval=doAll, warning=FALSE----------------------------------
# Mostramos los genes con la expresión diferencial más significativa:
TT2 <- topTable(ajusteC2, coef = 2, adjust.method ="BH", number = Inf)
head(TT2, 3)

## ---- Sig3_limma2, eval=doAll, warning=FALSE----------------------------------
# Mostramos los genes con la expresión diferencial más significativa:
TT3 <- topTable(ajusteC2, coef = 3, adjust.method ="BH", number = Inf)
head(TT3, 3)


## ---- Sig1_limma2_relax, eval=doAll, warning=FALSE----------------------------
sum(TT1$adj.P.Val < 0.1)

## ---- Sig1_limma2_num, eval=doAll, warning=FALSE------------------------------
head(TT1, 14)

## ---- Data_frame2, eval=doAll, warning=FALSE----------------------------------
TT1_ord = topTable(ajusteC2, coef = 1, adjust.method ="BH",
                   sort.by = "none", number = Inf)
TT2_ord = topTable(ajusteC2, coef = 2, adjust.method ="BH",
                   sort.by = "none", number = Inf)
TT3_ord = topTable(ajusteC2, coef = 3, adjust.method ="BH",
                   sort.by = "none", number = Inf)

#Creación del data frame con los p valores y los p valores ajustados.
df2 =cbind(fData(geod75759_rma), Pvalue.dif1 = TT1_ord$P.Value, 
           adj.Pvalue.dif1 = TT1_ord$adj.P.Val, Pvalue.dif2 = TT2_ord$P.Value,
           adj.Pvalue.dif2 = TT2_ord$adj.P.Val, Pvalue.difT = TT3_ord$P.Value,
           adj.Pvalue.difT = TT3_ord$adj.P.Val)
head(df2)


## ---- Parametros_Informe2, eval=doAll, warning=FALSE--------------------------
final_data = df2[df2$adj.Pvalue.dif1 < 0.1,]

ID = ifelse(is.na(final_data$PROBEID), NA, final_data$PROBEID)

entrezid = final_data$ENTREZID
ENTREZID= ifelse(is.na(entrezid), NA,
                 paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=",
                  entrezid,"'>",entrezid,"</a>"))

ensemblT = final_data$ENSEMBLTRANS
ENSEMBLTRANS = ifelse(is.na(ensemblT), NA,
                      paste0("<a href='http://www.ensembl.org/",
                             "Macaca_mulatta/Transcript/Idhistory?t=",
                             ensemblT ,"'>",ensemblT,"</a>"))

ensembl = final_data$ENSAMBL
ENSAMBL = ifelse(is.na(ensembl), NA,
                      paste0("<a href='http://www.ensembl.org/Macaca_mulatta",
                             "/Gene/Summary?db=core;g=",
                             ensembl ,"'>",ensembl,"</a>"))

GENSYMBOL = ifelse(is.na(final_data$GENSYMBOL), NA, final_data$GENSYMBOL)


## ---- Informe2, eval=doAll, warning=FALSE-------------------------------------
info2 = data.frame(ID=as.character(final_data$PROBEID), EntrezID = ENTREZID, 
                  EnsemblT=ENSEMBLTRANS,
                  EnsemblGen = ENSAMBL, GenSymbol = GENSYMBOL,
                  Pvalue.dif1 = final_data$Pvalue.dif1,
                  adj.Pvalue.dif1 = final_data$adj.Pvalue.dif1,
                  Pvalue.dif2 = final_data$Pvalue.dif2,
                  adj.Pvalue.dif2 = final_data$adj.Pvalue.dif2,
                  Pvalue.difT = final_data$Pvalue.difT,
                  adj.Pvalue.difT = final_data$adj.Pvalue.difT,
                  stringsAsFactors = FALSE )

f_out2 = "T2_limma.DE"
report_dir2 = "Reports/"
titulo2 = paste0("geod75759"," : Analisis limma para dos covariables.")

report2 = ReportingTools::HTMLReport(shortName = f_out2, title = titulo2,
                                     reportDirectory = report_dir2)
publish(info2,report2)
finish(report2)

