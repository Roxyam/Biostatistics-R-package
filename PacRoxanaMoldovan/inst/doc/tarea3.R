## ----Parametros ,include=FALSE------------------------------------------------
doAll = TRUE #Lo que debemos hacer cada vez que se ejecuta
done = FALSE #Marcamos así los pasos ya terminado
doafter = TRUE #Pasos finales una vez tenemos los datos
dir_datos = "../data/"

## ---- Librerias, eval=doafter, warning=FALSE----------------------------------
# Carga de paquetes:
pacman::p_load(Rsamtools,GenomicFeatures, GenomicAlignments,S4Vectors)

## ----Conteo, eval=done, warning=FALSE-----------------------------------------
#  
#  # Listamos los fichero sort.bam con su ruta:
#  indir = getwd()
#  fich_bam <- list.files(path = paste0(indir, "/PRJNA821482"),
#                         pattern = "*.sort.bam", full.names = TRUE)
#  
#  #Creamos una referencia a estos:
#  bamLst = Rsamtools::BamFileList(fich_bam, index=character(),
#                                  yieldSize=100000, obeyQname = TRUE)
#  
#  
#  # Fijamos la ruta del fichero genes.gtf
#  ruta_gtf= paste0(indir,"/PRJNA821482/referencia/Homo_sapiens/NCBI/GRCh38",
#                   "/Annotation/Archives/archive-2015-08-11-09-31-31",
#                   "/Genes/genes.gtf")
#  
#  # Generamos un objeto de tipo TxDb a partir de estos:
#  txdb = makeTxDbFromGFF(ruta_gtf, format="gtf", organism = "Homo sapiens")
#  
#  # Obtenemos los datos de anotación:
#  genes = exonsBy(txdb, by="gene")
#  
#  # Relacionamos las lecturas con las características genómicas:
#  PRJNA821482 <- summarizeOverlaps(features = genes,
#                                   read = bamLst,
#                                   mode= "Union",
#                                   singleEnd=FALSE, # Usa pair-ends
#                                   ignore.strand=TRUE,
#                                   fragments=FALSE)
#  
#  save(PRJNA821482, file=(paste0(dir_datos,"PRJNA821482_brutos.rda")))

## ----Carga_datos, eval=doafter, warning=FALSE, echo=F-------------------------
# Cargamos los datos:
load(paste0(dir_datos,"PRJNA821482_brutos.rda"))


## ----Dim, eval=doafter, warning=FALSE-----------------------------------------
pacman::p_load(SummarizedExperiment, GenomicRanges)
class(PRJNA821482)
dim(PRJNA821482)

## ----Head, eval=doafter, warning=FALSE----------------------------------------
head(assay(PRJNA821482),2)

## ----Variables_fenotipicas, eval=doafter, warning=FALSE-----------------------

# Obtenemos los metadatos:
metadat <- read.table( paste0(dir_datos, "SummarizedExperiment/SraRunTable.txt"),
                      sep = ",", header = T)

# Creamos las variables de tipo factor
TipoCelular <- factor(c(0,0,0,0,0,0,1,1,1,1,1,1), levels = 0:1,
                     labels = c("HSC", "K562"))
Tratamiento <- factor(c(0,0,0,1,1,1,0,0,0,1,1,1), levels = 0:1,
                      labels = c("Hypoxia", "Normoxia"))
# Añadimos la información a colData
d <- DataFrame(Run=metadat$Run, Muestra = metadat$Library.Name,
               TipoCelular=TipoCelular, Tratamiento=Tratamiento,  
               BioProject = metadat$BioProject )
colData(PRJNA821482, use.names=T) = d

# Asignamos estos nombres a las columnas
colnames(PRJNA821482) <- colData(PRJNA821482)$Run

head(colData(PRJNA821482), 1)


## ----Expermient_data, eval=doafter--------------------------------------------
info_ex <- list(name = "PRJNA821482", organism = "Homo sapiens",
                type = "Expression profiling by high throughput sequencing",
                lab = "Beijing Institute of Radiation Medicine",
                contact ="quanc1989@gmail.com",
                title = paste0("Hypoxia promotes erythroid differentiation",
                " throughautophagy induced by the HIF/REDD1/mTORC1 signaling",
                " axis"),
                abstract = paste0("For hematopoietic stem and progenitor cells",
                " (HSPCs), hypoxia is a specific microenvironment known as the",
                " hypoxic niche. The role of hypoxia in maintaining quiescence,",
                " self-renewal and proliferation of HSPCs has been more",
                " elaborated. What remains unclear is   the regulation of",
                " hypoxia on HSPC differentiation, especially on the",
                " differentiation of HSPCs into erythrocytes. In this study,",
                " we show that hypoxia evidently accelerates erythroid",
                " differentiation, and in this process effectively enhances",
                " autophagy which rapidly supplies cellular components through",
                " degradation of organelles or proteins to meet the needs of",
                " cell stress or development. In order to ascertain whether",
                " autophagy mediates the erythroid differentiation promoted",
                " by hypoxia, we analyze the perturbation of erythroid",
                " differentiation after pharmacological and genetic",
                " interference with autophagy, and address that autophagy is",
                " required for hypoxia-accelerated erythroid differentiation.",
                " Transcriptomics reveals HIF-1 and mTORC1 pathways may",
                " function as the upstream of autophagy to regulate erythroid",
                " differentiation under hypoxia. We further determine that",
                " mTORC1 signaling is suppressed by hypoxia to relieve its",
                " inhibition of autophagy. In addition, we also discover a",
                " new regulatory pattern that with the process of erythroid",
                " differentiation, mTORC1 activity gradually decreases and",
                " autophagy activity increases concomitantly. Moreover, we",
                " provide evidence that HIF-1 target gene REDD1, which",
                " inhibits mTORC1 activity, is dramatically upregulated under",
                " hypoxia. Silencing REDD1 expression activates mTORC1",
                " signaling and impairs the enhanced autophagy and the",
                " accelerated erythroid differentiation under hypoxia.",
                " Together, our study reveals that hypoxia is conductive to",
                " accelerating erythroid differentiation, in which autophagy",
                " induced by the HIF-1/REDD1/mTORC1 pathway plays a pivotal",
                " role."),
                url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199778")
metadata(PRJNA821482) <- info_ex


## ----Anotacion, eval=doafter,warning=FALSE------------------------------------
pacman::p_load(AnnotationDbi,org.Hs.eg.db)

# Tomamos otros identificadores
rowData(PRJNA821482)
anot <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(PRJNA821482),
                             columns = c("ENTREZID","ENSEMBL"),
                              keytype = "ALIAS")

# Eliminamos las sondas y genes repetidos:
a <- match(unique(anot$ENTREZID),anot$ENTREZID)
datos_anot1 <- anot[a,]

b <- match(unique(datos_anot1$ALIAS),datos_anot1$ALIAS)
datos_anot2 <- datos_anot1[b,]

# Hacemos la correspondencia entre los identificadores y las columnas
indice = match(rownames(PRJNA821482), datos_anot2$ALIAS)
rowData(PRJNA821482) = datos_anot2[indice,]

# Asignamos los datos y eliminamos las sondas que no se corresponden con genes:
PRJNA821482 <- PRJNA821482[which(!is.na(rowData(PRJNA821482)$ENTREZID)),]

## ----Info, eval=doafter, warning=FALSE----------------------------------------
(PRJNA821482)

## ----Guardado_final, eval=doafter, warning=FALSE------------------------------
save(PRJNA821482, file=(paste0(dir_datos,"PRJNA821482.rda")))

