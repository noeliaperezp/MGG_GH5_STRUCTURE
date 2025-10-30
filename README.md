
# MGG GH5 - ANÁLISIS DE ESTRUCTURA POBLACIONAL 

El estudio de la estructura de una población ofrece información valiosa sobre su historia y origen. En esta sesión práctica analizaremos la estructura y el flujo génico de poblaciones humanas utilizando herramientas estadísticas muy comunes en genética de poblaciones, como los análisis de componentes principales (PCA) y los modelos de mezcla (ej., ADMIXTURE o STRUCTURE).

El método de componentes principales, o PCA, busca identificar los principales ejes de variación (i.e., componentes principales) dentro de un conjunto de datos, asumiendo independencia entre ellos. El primer componente principal (PC1) captura la mayor diferencia genética entre individuos o poblaciones, el PC2 la siguiente mayor diferencia, y así sucesivamente hasta explicar la mayor parte de la variación genética total. A menudo se visualiza mediante un diagrama de dispersión bidimensional, donde cada punto representa un individuo y los puntos que son genéticamente similares se agrupan más próximos entre sí. Los análisis de PCA permiten detectar aislamiento por distancia, reflejando en ocasiones la distribución geográfica. Estos análisis son efectivos cuando una gran proporción de la variación es capturada por un número pequeño de componentes principales.

ADMIXTURE y STRUCTURE utilizan el mismo modelo estadístico, pero como ADMIXTURE permite analizar un mayor número de muestras y marcadores (SNPs) en menos tiempo, nos centraremos en él. ADMIXTURE utiliza el método de máxima verosimilitud para asignar individuos a K subpoblaciones ancestrales, donde el número K se especifica de antemano. En cada subpoblación (demo) se asume apareamiento aleatorio. Cada subpoblación se caracteriza por un conjunto concreto de frecuencias alélicas para cada locus, y se asume que no hay desequilibrio de ligamiento (LD) entre marcadores. Por este motivo se recomienda evitar los SNPs estrechamente ligados para garantizar que el LD que intenta minimizar el programa se deba a la estructura poblacional y no a la proximidad física de los marcadores. El modelo asigna probabilidades a cada individuo de pertenecer a cada subpoblación asumida. Un mismo individuo puede tener probabilidades mayores de cero de pertenecer a varias subpoblaciones, indicando posible mezcla o flujo genético.

** FST

## Datos

Trabajaremos con datos de la fase 3 del Proyecto 1000 Genomas (https://www.internationalgenome.org/1000-genomes-summary/). Se trata de una iniciativa internacional, desarrollada entre el 2008 y 2015 con el objetivo de crear un catálogo completo de la variación genética humana mediante la secuenciación del genoma completo de un conjunto diverso de individuos de diferentes poblaciones. Este proyecto proporciona un marco de referencia para comprender nuestra historia evolutiva, la diferenciación entre grupos y las implicaciones biomédicas de dicha variación.

En el transcurso de la práctica trabajaremos principalmente con dos conjuntos de datos:

1.	Poblaciones humanas con muestreo geográfico discreto.
2.	Poblaciones humanas con muestreo geográfico ‘continuo’ o con mayor resolución espacial.
   
Encontrarás información completa sobre las poblaciones y muestras que cuentan con datos disponibles en el portal de datos del IGSR (“The International Genome Sample Resource”, https://www.internationalgenome.org/data-portal/population). El IGSR se estableció con el fin de garantizar el acceso y uso de los datos del Proyecto 1000 Genomas, actualizar los recursos con el ensamblaje de referencia actual, y ampliar el conjunto de datos con nuevos datos procedentes de las muestras del Proyecto 1000 Genomas y de nuevas poblaciones.
Debido a que disponemos de un tiempo limitado, partiremos de una muestra de 20 o 50 individuos ya seleccionados previamente de forma aleatoria por población (Tabla 1), y nos limitaremos a las variantes genéticas localizadas en el cromosoma 22.  

Tabla 1. Poblaciones humanas evaluadas por caso de estudio.

| Datos | Código| Población | Superpoblación / Ascendencia | Individuos |
| ----- | ----- | --------- | ---------------------------- | ---------- |
| Muestreo geográfico discreto (chr22_pop_dist)	 | YRI | Yoruba | África | 20 |
|  | LWK | Luhya | África | 20 |
|  | GBR | Británica | Europa | 20 |
|  | GIH | Gujarati | Asia meridional | 20 |
|  | CHB | Han Chinese | Asia Oriental | 20 |				
|  | PEL | Peruana | América | 20 |	
| Muestreo geográfico ‘continuo’ (chr22_pop_cont) | GBR | Británica | Europa | 50 |
|  | IBS | Ibérica | Europa | 50 |
|  | TSI | Toscana | Europa | 50 |			
|  | FIN | Finlandés | Europa | 50 |			
|  | GIH | Gujarati | Asia meridional | 50 |			


Tan sólo a modo informativo, los archivos VCF (Variant Call Format) de cada cromosoma, con información para el total de los 2504 individuos secuenciados y pertenecientes a 26 poblaciones diferentes, pueden descargarse desde el sitio FTP:  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/. Desde la terminal, utilizaríamos directamente el comando wget seguido de la URL del archivo que se quiere descargar (ej. para descargar el archivo VCF completo del cromosoma 22: wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz).


## Descarga de datos

Tras conectarte al CESGA, muévete al directorio de trabajo creado para esta sesión práctica.

```
cd $LUSTRE/MGG_GH5/STRUCTURE
```
Los datos y scripts con los que trabajaremos están almacenados en este mismo repositorio, creado para esta sesión práctica. Para clonar el repositorio desde la terminal:

```
git clone https://github.com/noeliaperezp/MGG_GH5_STRUCTURE.git .
ls
```

Verás que se han descargado una serie de archivos:
- Scripts con extensión “*.sh*” listos para ejecutar y que incluyen cada uno de los comandos que iremos utilizando paso a paso a lo largo de la práctica.
- Un directorio “data” que contiene los archivos de variantes (VCF) para los dos escenarios que analizaremos (“chr22_pop_dist.vcf.gz” y “chr22_pop_finer.vcf.gz”). En este directorio también encontrarás dos archivos adicionales con los identificadores de las muestras y el código de la población de origen (“poplist_human_pop_dist.txt” y “poplist_human_pop_finer.txt”).
- Un script de R (“structure_figures.R”) con el código que utilizaremos para representar los resultados.
- Un archivo README con el guion de prácticas para esta sesión.

Ahora que tenemos todo lo necesario, empezaremos analizando la estructura poblacional en el conjunto de datos muestreados de forma discreta (**chr22_pop_dist**).


## Filtrado

Carga los módulos necesarios:
```
module load cesga/2020 gcc/system plink/2.00a2.3
```

Crea un nuevo directorio y muévete a él:
```
mkdir data_pruned
cd data_pruned/
```

Retén variantes bialélicas y genera archivos plink (.bed,.fam,.bim):
```
plink2 --vcf ../data/chr22_pop_dist.vcf.gz --make-bed --max-alleles 2 --snps-only --out chr22_pop_dist.SNPs
```

Asigna nuevos identificadores a variantes sin nombre:
```
plink2 --bfile chr22_pop_dist.SNPs --set-missing-var-ids @:# --make-bed --out chr22_pop_dist.SNPid
rm chr22_pop_dist.SNPs*
```

Filtra variantes y muestras con >10% missing calls y maf < 0.05:
```
plink2 --bfile chr22_pop_dist.SNPid --geno 0.1 --mind 0.1 --maf 0.05 --make-bed --out chr22_pop_dist.filtered
rm chr22_pop_dist.SNPid*
```

Filtra variantes con fuerte desequilibrio de ligamiento (LD):
```
plink2 --bfile chr22_pop_dist.filtered --indep-pairwise 50 5 0.9 --out chr22_pop_dist.plink2
plink2 --bfile chr22_pop_dist.filtered --extract chr22_pop_dist.plink2.prune.in --make-bed --out chr22_pop_dist_pruned
rm chr22_pop_dist.filtered*
rm *.in
rm *.out
rm *.log
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/STRUCTURE”):
```
cd ..
pwd
```

De aquí en adelante trabajaremos con los archivos filtrados (“*_pruned*”), que encontrarás en el directorio "*data_pruned/*". 

¿Cuántas muestras y variantes se han eliminado y cuántos retenemos? ¿son suficientes?


## Análisis de componentes principales (PCA)

Crea un nuevo directorio y muévete a él:
```
if [ ! -d results_pca/chr22_pop_dist_pruned ] 
then 
 mkdir -p results_pca/chr22_pop_dist_pruned
fi
cd results_pca/chr22_pop_dist_pruned
```

PCA:
```
plink2 --bfile ../../data_pruned/ chr22_pop_dist_pruned --pca 10 –out chr22_pop_dist_pruned_pca 
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/STRUCTURE”):
```
cd ../..
pwd 
```

En el directorio "*results_pca/*" encontrarás dos outputs, "*.eigenval*" y "*.eigenvec*".

**Transferencia de ficheros del CESGA a tu sistema local**
Representaremos los resultados del PCA en RStudio. Para ello necesitarás transferir los resultados del PCA a tu sistema local. También necesitarás un archivo adicional llamado poplist_human_pop_dist.txt, que contiene una lista de los individuos analizados junto con el código de su población de origen. Para transferir estos archivos puedes usar WinSCP, o alternativamente usando el comando scp desde tu máquina local: 

En tu sistema local, crea y abre la carpeta a la que quieras exportar los archivos. Haz clic con el botón derecho y abre una terminal (ej. Windows PowerShell). 

Ten en cuenta que necesitarás conocer la ruta completa dentro del CESGA que conduce al archivo que quieres transferir. Para ello ejecuta el comando pwd en la terminal del CESGA. Asegúrate de que te encuentras en el directorio “*…/MGG_GH5/STRUCTURE*”, punto base de la estructura de directorios que hemos creado para esta sesión práctica.  

Ahora escribe en la terminal local (ej. Windows PowerShell) los siguientes comandos. Te pedirá confirmar la autenticidad del host, y a continuación tu contraseña de usuario en el CESGA.

Una vez copiado en siguiente código recuerda cambiar:
- *username* por tu nombre de usuario
- */path* por la ruta que te devuelve *pwd* en el CESGA

```
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_pca/chr22_pop_dist_pruned/chr22_pop_dist_pruned_pca.eigenval .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_pca/chr22_pop_dist_pruned/chr22_pop_dist_pruned_pca.eigenvec .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/data/poplist_human_pop_dist.txt .
```

Una vez transferidos los resultados del PCA, abre RStudio y muévete al directorio que has creado con los resultados. Para cambiar de directorio utiliza el comando 

```
setwd(“your/path/”)
```

Presionando tabulador tendrás una guía de los directorios a los que te puedes mover. Una vez en la carpeta copia y ejecuta las líneas de código que encontrarás en el apartado "PCA" en el script structure_figures.R.

> **_PREGUNTA:_** ¿Qué porcentaje de la varianza es explicada por PC1 y PC2 juntos? ¿Qué separa cada PC?


## Análisis de mezcla con ADMIXTURE

Carga los módulos necesarios:
```
module load cesga/2020 admixture/1.3.0
```

Crea un nuevo directorio y muévete a él:
```
if [ ! -d results_admixture/chr22_pop_dist_pruned ] 
then 
 mkdir -p results_admixture/chr22_pop_dist_pruned
fi
cd results_admixture/chr22_pop_dist_pruned
```

ADMIXTURE:
```
NPOP=6
for (( K = 1; K <= $NPOP; K++ ))
do
  admixture --cv ../data_pruned/chr22_pop_dist_pruned.bed $K | tee chr22_pop_dist_pruned_admixture_K${K}.out
done
```

Extrae el cross-validation error(CV) para cada valor de K:
```
awk '/CV/ {print $3,$4}' *.out | cut -c 4,7-20 > chr22_pop_dist_pruned.admixture.cv.error
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/STRUCTURE”):
```
cd ../..
pwd 
```

En el directorio "*results_admixture/*" encontrarás tres outputs, "*.Q*", que contiene las fracciones de ascendencia, "*.P*", con las fecuencias alélicas de las poblaciones ancestrales inferidas, y "*.cv.error*" con los errores para cada valor de K.

**Transferencia de ficheros del CESGA a tu sistema local**
De nuevo, representaremos los resultados del ADMIXTURE en RStudio. Para ello, en la terminal de tu sistema local escribe el siguiente código, cambiando previamente:

- *username* por tu nombre de usuario
- */path* por la ruta que te devuelve *pwd* en el CESGA

```
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_admixture/chr22_pop_dist_pruned/chr22_pop_dist_pruned.admixture.cv.error .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_admixture/chr22_pop_dist_pruned/*P .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_admixture/chr22_pop_dist_pruned/*Q .
```

Una vez transferidos los resultados del ADMIXTURE, abre RStudio y asegúrate de que te encuentras en la carpeta que contiene los resultados. A continuación copia y ejecuta las líneas de código que encontrarás en el apartado "ADMIXTURE" en el script structure_figures.R.

> **_PREGUNTA:_** ¿Qué significa cada K? ¿Qué valor de K deberíamos utilizar?

> **_PREGUNTA:_** ¿Qué conclusiones se sacan del ADMIXTURE?


## Diferenciación poblacional: cálculo de F<sub>ST

Carga los módulos necesarios:
```
module load cesga/2020 plink/1.9b5
```

Crea un nuevo directorio y muévete a él:
```
if [ ! -d results_fst/chr22_pop_dist_pruned ] 
then 
 mkdir -p results_fst/chr22_pop_dist_pruned
fi
cd results_fst/chr22_pop_dist_pruned
```

Calcula el F<sub>ST global y por marcador:
```
plink --bfile ../../data_pruned/chr22_pop_dist_pruned --fst --within ../../data/ poplist_human_pop_dist.txt --out chr22_pop_dist_pruned_fst
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/STRUCTURE”):
```
cd ../..
pwd 
```

En el directorio "*results_fst/*" encontrarás un output principal, "*.fst*", que contiene las estimas por marcador, y un archivo "*.log*" que contiene el valor F<sub>ST promedio.

> **_PREGUNTA:_** ¿Qué nivel de diferenciación genética encontramos entre las poblaciones evaluadas?


## Ejecución de comandos a través de scripts

Para analizar la estructura poblacional del conjunto de datos muestreados de forma 'continua', o con mayor resolución espacial (**chr22_pop_finer**, seguiremos los mismos pasos y lanzaremos los mismos comandos, en esta ocasión desde un script listo para ejecutarse. En el CESGA, el comando *sbatch* se utiliza para enviar trabajos (ej. ejecución de scripts) al sistema de colas. 

Filtrado de datos:
```
sbatch script_data_filtering.sh chr22_pop_finer
```

Puedes comprobar el estado de ejecución de trabajo utilizando el comando
```
squeue
```

Una vez finalizado, lanza los siguientes scripts.

PCA:
```
sbatch script_ run_pca.sh chr22_pop_finer_pruned

```

ADMIXTURE (ten en cuenta que ahora el conjunto de datos incluye muestras de 5 poblaciones diferentes):
```
sbatch script_ run_admixture.sh chr22_pop_finer_pruned 5
```

F<sub>ST:
```
sbatch script_ calculate_fst.sh chr22_pop_finer_pruned data/poplist_human_pop_finer
```

Una vez finalizados, puedes transferir los resultados a tu sistema local para la elaboración de figuras. Recuerda cambiar:
- *username* por tu nombre de usuario
- */path* por la ruta que te devuelve *pwd* en el CESGA

```
# Resultados del PCA
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_pca/chr22_pop_finer_pruned/chr22_pop_finer_pruned_pca.eigenval .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_pca/chr22_pop_finer_pruned/chr22_pop_finer_pruned_pca.eigenvec .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/data/poplist_human_pop_finer.txt .

# Resultados del ADMIXTURE
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_admixture/chr22_pop_finer_pruned/chr22_pop_finer_pruned.admixture.cv.error .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_admixture/chr22_pop_finer_pruned/*P .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_admixture/chr22_pop_finer_pruned/*Q .
```

> **_PREGUNTA:_** ¿Qué porcentaje de la varianza es explicada por PC1 y PC2 juntos? ¿Qué separa cada PC?

> **_PREGUNTA:_** ¿Qué conclusiones se sacan del ADMIXTURE?

> **_PREGUNTA:_** ¿Qué nivel de diferenciación genética encontramos ahora entre las poblaciones evaluadas?

