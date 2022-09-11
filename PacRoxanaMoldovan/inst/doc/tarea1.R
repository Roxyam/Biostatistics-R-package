## ----setup,include=FALSE------------------------------------------------------
# Parámetros iniciales:
doAll = TRUE
done = FALSE
dir_datos = "../data/"
dir_fig = "../Figuras/"

## ---- Paquetes, eval=doAll, warning=FALSE-------------------------------------
library(PacRoxanaMoldovan)
# Paquetes directamente relacionados con los datos y su procesamiento:
pacman::p_load(ArrayExpress, Biobase, oligo, pd.cyngene.1.0.st)
# Otros paquetes empleados:
pacman::p_load(magick, stringr, jsonlite)

## ---- Descarga, eval=done, message=FALSE, warning=FALSE-----------------------
#  # Descarga de datos desde ArrayExpress:
#  geod75759_raw = ArrayExpress("E-GEOD-75759")
#  # Guardado de datos:
#  f_out = "geod75759_raw.rda"
#  save(geod75759_raw, file = paste0(dir_datos,f_out))

## ---- Carga_datos, eval=doAll, message=FALSE, warning=FALSE-------------------
# Carga de datos del fichero
load(paste0(dir_datos,"geod75759_raw.rda"))

## ----Clase, eval=doAll--------------------------------------------------------
  class(geod75759_raw)

## ----Plataforma Anotacion, eval=doAll-----------------------------------------
  cat("La plataforma empleada para la anotación es",
      annotation(geod75759_raw), "\n")

## ----Nº Sondas y muestras, eval=doAll-----------------------------------------
n_sondas <- (dim(exprs(geod75759_raw)))[1]
n_muestras <- dim(exprs(geod75759_raw))[2]
  cat("El experimento presenta ", n_sondas, " sondas que se corresponden con ",
      n_muestras, " muestras.")


## ----Image, eval=doAll, message=FALSE, warning=FALSE--------------------------
n_muestra <- sample(1:n_muestras, 1) # Escogemos una muestra al azar
png(paste0(dir_fig,"microarray", n_muestra, "_image.png"))
image(geod75759_raw[,n_muestra])
dev.off()

## ----Fig_1, eval=doAll, message=FALSE, warning=FALSE, fig.margin = TRUE, fig.width=3.5, fig.height=3.5, fig.cap = "Fig1. Imagen digital de los niveles de expresión de una muestra perteneciente a geod75759."----
fig1 <- knitr:: include_graphics(paste0(dir_fig,"microarray",
                                        n_muestra, "_image.png"))
fig1

## ---- MA-plot, eval=doAll, message=FALSE, warning=FALSE-----------------------
# Representación y guardado de los gráficos:
for (i in 1:n_muestras){
  png(paste0(dir_fig,"MAplot_geod75759_raw_", i, ".png"))
  oligo :: MAplot(geod75759_raw, which=i,plot.method="smoothScatter")
  dev.off()
}

## ----Fig_2, eval=doAll, message=FALSE, warning=FALSE, echo=FALSE, fig.margin = TRUE, fig.cap="Fig2. MA-plot de las 18 muestras del estudio geod75759 antes del procesamiento."----

figuras <- list.files(path = dir_fig, pattern = "MAplot_geod75759_raw_")
# Abrimos todas las filas:
m_plots = c()
for (fig in figuras){
  f <- magick::image_read(paste0(dir_fig,fig))
  m_plots = c(m_plots, image_scale(f, 225))

}

# Agrupamos las figuras por filas:
f1 <- image_append(c(m_plots[[1]], m_plots[[2]], m_plots[[3]]))
f2 <- image_append(c(m_plots[[4]], m_plots[[5]], m_plots[[6]]))
f3 <- image_append(c(m_plots[[7]], m_plots[[8]], m_plots[[9]]))
f4 <- image_append(c(m_plots[[10]], m_plots[[11]], m_plots[[12]]))
f5 <- image_append(c(m_plots[[13]], m_plots[[14]], m_plots[[15]]))
f6 <- image_append(c(m_plots[[16]], m_plots[[17]], m_plots[[18]]))

image_append(c(f1,f2,f3,f4,f5,f6), stack = TRUE)


## ---- Densidad, eval=doAll, message=FALSE, warning=FALSE----------------------

# Representación y guardado de los estimadores de densidad:
png(paste0(dir_fig,"Densidad_geod75759_raw.png"))
oligo::hist(geod75759_raw, which="all",
            main= "A.Estimadores kernel de geod75759 sin normalizar.",
            ylab="Densidad", xlab="Log2 Intensidad", xaxt="n")
dev.off()

## ---- BoxPlot, eval=doAll, message=FALSE, warning=FALSE-----------------------
# Representación y guardado del box-plot:
png(paste0(dir_fig,"BoxPlot_geod75759_raw.png"))
oligo::boxplot(geod75759_raw, which = "all", 
               main= "B.Boxplot geod75759 sin normalizar",
               ylab="Expresión", xlab="Muestras", xaxt="n")
               axis(1, at=c(1:18),labels = paste0("M",c(1:18)))
dev.off()

## ----Fig_3, eval=doAll, message=FALSE, warning=FALSE, echo=FALSE, fig.margin = TRUE, fig.cap=" Fig3. Evaluación de la variabilidad de los niveles de expresión entre arrays en el experimento geod75759.  a.Estimador kernel de densidad de los niveles de intensidad entre microarrays. b. Diagrama de cajas con bigotes de los niveles de expresión a nivel de sonda."----
fig3a <- magick::image_read(paste0(dir_fig,"Densidad_geod75759_raw.png"))
fig3b <- magick::image_read(paste0(dir_fig,"BoxPlot_geod75759_raw.png"))
image_append(c(image_scale(fig3a,350), image_scale(fig3b,350)))

## ---- RMA, eval=doAll, message=FALSE, warning=FALSE---------------------------
# Normalización:
geod75759_rma <- rma(geod75759_raw)

## ----Nº Sondas y muestras 2, eval=doAll---------------------------------------
# Numero de sondas y de muestras:
n_sondas2 <- (dim(exprs(geod75759_rma)))[1]
n_muestras2 <- dim(exprs(geod75759_rma))[2]
  cat("Trans la normalización tenemos ", n_sondas2, 
      " sondas, que se corresponden con ", n_muestras2, " muestras.")

## ---- Densidad_rma, eval=doAll, message=FALSE, warning=FALSE------------------

# Representación y guardado de los estimadores de densidad:
png(paste0(dir_fig,"Densidad_geod75759_rma.png"))
oligo::hist(geod75759_rma,
            main= "A.Estimadores kernel de geod75759 tras normalizar.",
            ylab="Densidad", xlab="Log2 Intensidad", xaxt="n")
dev.off()

## ---- BoxPlot_rma, eval=doAll, message=FALSE, warning=FALSE-------------------
# Representación y guardado del box-plot:
png(paste0(dir_fig,"BoxPlot_geod75759_rma.png"))
oligo::boxplot(geod75759_rma, 
               main= "B.Boxplot geod75759 tras normalizar",
               ylab="Expresión", xlab="Muestras", xaxt="n")
               axis(1, at=c(1:18),labels = paste0("M",c(1:18)))
dev.off()

## ----Fig_4, eval=doAll, message=FALSE, warning=FALSE, echo=FALSE, fig.margin = TRUE, fig.cap=" Fig4. Evaluación de la variabilidad de los niveles de expresión entre arrays en el experimento geod75759.  a.Estimador kernel de densidad de los niveles de intensidad entre microarrays. b. Diagrama de cajas con bigotes de los niveles de expresión a nivel de sonda."----
fig3a <- magick::image_read(paste0(dir_fig,"Densidad_geod75759_rma.png"))
fig3b <- magick::image_read(paste0(dir_fig,"BoxPlot_geod75759_rma.png"))
image_append(c(image_scale(fig3a,350), image_scale(fig3b,350)))

## ----Expermient_data, eval=doAll----------------------------------------------
dat_ex <- new("MIAME", name = "E-GEOD-75759",
              lab = "Environmental Sciences and Engineering, UNC-Chapel Hill",
              contact ="rfry@unc.edu",
              title = paste0("Expression data from male cynomolgus macaques",
              " exposed to 6 ppm formaldehyde."),
              abstract = paste0("Exogenous formaldehyde disrupts genomic/",
              " epigenomic profiles in the rodent nose and white blood cells",
              " (WBCs) related to inflammation and immune signaling, although",
              " it does not reach the circulating blood. We aimed to compare",
              " and contrast alterations in genomic signaling in the nose and",
              " circulating blood of non-human primates exposed to",
              " formaldehyde. We used microarrays to identify transcripts",
              " that were diffentially expressed in response to formaldehyde",
              " inhalation exposure. A total of 14 primates received two",
              " consecutive days of 6-hour whole body inhalation exposures",
              " consisting of either filtered air (n = 6) or a target of 6",
              " ppm formaldehyde (n = 8). To assess formaldehyde-induced",
              " changes in genome-wide gene expression profiles, RNA samples",
              " extracted from the nasal epithelium and circulating WBCs were",
              " labeled and hybridized to the Affymetrix Cynomolgus Macaque",
              " Gene 1.0 ST Array."),
              url = paste0("https://www.ebi.ac.uk/arrayexpress/",
                           "experiments/E-GEOD-75759/"),
              other = list(notes = paste0("Released on 30 June 2016, ",
                                          "last updated on 24 July 2016")))
experimentData(geod75759_rma) = dat_ex
experimentData(geod75759_rma)


## ----Datos_fenotipicos, eval=doAll--------------------------------------------
# Mostramos la variable fenotipicas asociadas a nuestra muestra:
datos_fenotipicos <- pData(geod75759_rma)
head(datos_fenotipicos, n=1)

## ----Definicion_factores, eval=doAll------------------------------------------

# Factor que define el tipo de muestra:
tipo_muestra <- factor(c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1),
                       levels = 0:1, labels = c("WBC", "NE"))

# Factor que describe el tratamiento de la muestra:
tratamiento <- factor(c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1),levels = 0:1,
                      labels = c("Formaldehyde", "Filtered_Air"))

# Completamos los datos con los nuevos factores en el orden conveniente
datos_fenotipicos2 <- data.frame("Source.Name" = datos_fenotipicos$Source.Name,
                                 "Sample_type" = tipo_muestra, 
                                 "Sample_treatment" = tratamiento, 
                                 datos_fenotipicos[,2:37])

pData(geod75759_rma) <- datos_fenotipicos2

## ----Mostrar_pData, eval=doAll------------------------------------------------
# Mostramos los datos fenotípicos de la primera y última muestra:
(pData(geod75759_rma)[c(1,14),])

## ----Mostrar_factores, eval=doAll---------------------------------------------
# Descripción de la covariante Sample_type:
table(geod75759_rma$Sample_type)

# Descripción de la covariante Sample_type:
table(geod75759_rma$Sample_treatment)

## ---- Lectura_datos_anotacion, eval=doAll-------------------------------------
dir_expression_set = paste0(dir_datos, "ExpressionSet/")
platform = read.table(paste0(dir_expression_set,  
"TFS-Assets_LSG_Support-Files_CynGene-1_0-st-v1-na34-cm1-transcript-csv", "/",
                             "CynGene-1_0-st-v1.na34.cm1.transcript.csv"),
                      sep=",", header = TRUE)

## ---- Datos_anotacion_disponibles, eval=doAll---------------------------------

# Mostramos los identificadores de las columnas:
colnames(platform)

## ---- Itentificador_sonda, eval=doAll-----------------------------------------

# Mostramos los primeros identificadores.
head(rownames(fData(geod75759_rma)))
head(platform$probeset_id)

## ---- Correspondencia_anotacion, eval=doAll-----------------------------------
# Ahora que que tenemos la información de las sondas buscamos correspondencias:
ID = data.frame("id"=as.numeric(featureNames(geod75759_rma)))
anot <- merge(ID, platform, by.x="id", by.y="probeset_id", all.x=TRUE)

head(anot, 2)

## ----Anotacion1, eval=doAll, message=FALSE, warning=FALSE---------------------

# Eliminamos las sondas sin gene_assignment (no se corresponden con genes):
anot2 <- anot[!(anot$gene_assignment == "---"),]

# Nos centramos en el gene_assigment que tiene la info que nos interesa:
# Múltiples anotaciones par aun mismo gen estas separadas por ' /// ':
opciones_anot <- (str_split(anot2$gene_assignment, " /// "))


# Creamos el data frame que contendrá la información:
datos_anot <- data.frame("accession" = character(), "gene_symbol" = character(),
                         "gene_title" = character(), "cytoband"== character(), 
                         "entrez_gene_id"=numeric())

# Nos quedamos con las anotaciones de Ensembl, que comienzan por 'ENS':
for (x in (1:length(opciones_anot))){ # Recorremos la lista de todas las sondas.
  opcion = opciones_anot[[x]]

  for (i in opcion){ # Recorremos las posibles anotaciones.
    
    # Comprobamos si los primeros caracteres coinciden con 'ENS':
    j = substr(i, start = 1, stop = 3)
    if (j == "ENS"){
      
      # Separamos la información de la anotación por ' // '
      inf = (str_split(i, " // "))[[1]]
      
      #Descartamos las sondas que tengan '---' en el elemento 5.
      if (inf[5]!= "---"){
        
        # Tomamos el identificador de dicha posición:
        id = anot2[x,1]
        nuevo <- data.frame (PROBEID = id, ENTREZID = inf[5],
                             ENSEMBLTRANS = inf[1] , GENSYMBOL = inf[2])
        # Lo añadimos al data frame las columnas deseadas:
        datos_anot <- rbind(datos_anot,nuevo)
      }
    }
  }
}

## ----Visualizazión_datos_anot, eval=doAll-------------------------------------
head(datos_anot)


## ----Visualizazión_ensembl, eval=doAll----------------------------------------
head(datos_anot$ENSEMBLTRANS)


## ----Descarga_VGNC, eval=doAll------------------------------------------------
db <- fromJSON(paste0(dir_expression_set, "macaque_vgnc_gene_set_All.json"))
colnames(db)

## ----Correspondencia_VGNC, eval=doAll-----------------------------------------
a <- (merge(datos_anot, db, by.x="ENTREZID", by.y="ncbi_id", all.x=TRUE))
datos_anot <- data.frame(a[,1:4], ENSAMBL=a$ensembl_gene_id)

## ----Visualizazión_anotacion, eval=doAll--------------------------------------
head(datos_anot)

## ----Comprobacion_sondas_unicas, eval=doAll-----------------------------------
p <- table(datos_anot$PROBEID)>1
table(p)

p2 <- which(p)

## ----Comprobacion_sondas_unicas2, eval=doAll----------------------------------
  cat("Si, encontramos", length(p2),
      "sondas que se corresponden a más de un gen.")

## ----Comprobacion_sondas_unicas3, eval=doAll----------------------------------
head (p2)

## ----Comprobacion_genes_unicos, eval=doAll------------------------------------
g <- table(datos_anot$ENTREZID)>1
table(g)
g2 <- which(g)

## ----Comprobacion_genes_unicos2, eval=doAll-----------------------------------
  cat("Si, encontramos", length(g2),
      "genes que se corresponden a más de una sonda.")

## ----Comprobacion_genes_unicos3, eval=doAll-----------------------------------
head (g2)

## ----Eliminar_repeticiones, eval=doAll----------------------------------------
# Eliminamos sondas y genes repetidos:
a = match(unique(datos_anot$PROBEID),datos_anot$PROBEID)
datos_anot2 <- datos_anot[a,]
b = match(unique(datos_anot2$ENTREZID),datos_anot2$ENTREZID)
datos_anot3 <- datos_anot2[b,]

## ---- eval=doAll--------------------------------------------------------------
length(unique(datos_anot3$PROBEID)) == length(unique(datos_anot3$ENTREZID))

## ---- eval=doAll--------------------------------------------------------------
# Finalmente vemos si según los identificadores ENSEMBL hay repeticiones:
e <- (table(datos_anot3$ENSEMBL)>1)
table(e)


## ----Asignar_fData, eval=doAll------------------------------------------------
# Asignamos los datos a fData:
c <- match(featureNames(geod75759_rma),datos_anot3$PROBEID)
fData(geod75759_rma) <- datos_anot3[c,]
# Eliminamos las sondas cuya información no está en los datos de anotación:
geod75759_rma <- geod75759_rma[which(!is.na(fData(geod75759_rma)$PROBEID),)]


## ---- Guardado_rma, eval=doAll, message=FALSE, warning=FALSE------------------
# Guardado de datos:
f_out = "geod75759_rma.rda"
save(geod75759_rma, file = paste0(dir_datos,f_out))

# Carga de datos del fichero
load(paste0(dir_datos,f_out))

