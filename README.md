
# MGG GH5 - ANÁLISIS DE ESTRUCTURA POBLACIONAL 

El estudio de la estructura de una población ofrece información valiosa sobre su historia y origen. En esta sesión práctica analizaremos la estructura y el flujo génico de poblaciones humanas utilizando herramientas estadísticas muy comunes en genética de poblaciones, como los análisis de componentes principales (PCA) y los modelos de mezcla (ej., ADMIXTURE o STRUCTURE).

El método de componentes principales, o PCA, busca identificar los principales ejes de variación (i.e., componentes principales) dentro de un conjunto de datos, asumiendo independencia entre ellos. El primer componente principal (PC1) captura la mayor diferencia genética entre individuos o poblaciones, el PC2 la siguiente mayor diferencia, y así sucesivamente hasta explicar la mayor parte de la variación genética total. A menudo se visualiza mediante un diagrama de dispersión bidimensional, donde cada punto representa un individuo y los puntos que son genéticamente similares se agrupan más próximos entre sí. Los análisis de PCA permiten detectar aislamiento por distancia, reflejando en ocasiones la distribución geográfica. Estos análisis son efectivos cuando una gran proporción de la variación es capturada por un número pequeño de componentes principales.

ADMIXTURE y STRUCTURE utilizan el mismo modelo estadístico, pero como ADMIXTURE permite analizar un mayor número de muestras y marcadores (SNPs) en menos tiempo, nos centraremos en él. ADMIXTURE utiliza el método de máxima verosimilitud para asignar individuos a K subpoblaciones ancestrales, donde el número K se especifica de antemano. En cada subpoblación (demo) se asume apareamiento aleatorio. Cada subpoblación se caracteriza por un conjunto concreto de frecuencias alélicas para cada locus, y se asume que no hay desequilibrio de ligamiento (LD) entre marcadores. Por este motivo se recomienda evitar los SNPs estrechamente ligados para garantizar que el LD que intenta minimizar el programa se deba a la estructura poblacional y no a la proximidad física de los marcadores. El modelo asigna probabilidades a cada individuo de pertenecer a cada subpoblación asumida. Un mismo individuo puede tener probabilidades mayores de cero de pertenecer a varias subpoblaciones, indicando posible mezcla o flujo genético.

Por otro lado, el F<sub>ST</sub> (Wright 1949), o distancia genética, es una de las principales medidas de diferenciación genética entre y dentro de poblaciones, originalmente definido como "*la correlación entre gametos seleccionados al azar dentro de la misma subpoblación en relación con la población total o, equivalentemente, como la desviación de las frecuencias genotípicas respecto a las expectativas de Hardy-Weinberg en relación con la población total*" (Lohmueller y Nielsen, 2021). Las estimas ofrecen información sobre la fracción de la variación en frecuencias alélicas debida a diferencias entre regiones o subpoblaciones, adscribiéndose la fracción restante a diferencias dentro de ellas. Si no hay subdivisión de la población, el F<sub>ST</sub> es igual a 0. Hay múltiples herramientas que permiten calcular este parámetro. En esta práctica utilizaremos el programa *plink*.


## Datos

Trabajaremos con datos de la fase 3 del Proyecto 1000 Genomas (https://www.internationalgenome.org/1000-genomes-summary/). Se trata de una iniciativa internacional, desarrollada entre el 2008 y 2015 con el objetivo de crear un catálogo completo de la variación genética humana mediante la secuenciación del genoma completo de un conjunto diverso de individuos de diferentes poblaciones. Este proyecto proporciona un marco de referencia para comprender nuestra historia evolutiva, la diferenciación entre grupos y las implicaciones biomédicas de dicha variación.

En el transcurso de la práctica trabajaremos principalmente con dos conjuntos de datos:

1.	Poblaciones humanas con muestreo geográfico discreto.
2.	Poblaciones humanas con muestreo geográfico ‘continuo’ o con mayor resolución espacial.
   
Encontrarás información completa sobre las poblaciones y muestras que cuentan con datos disponibles en el portal de datos del IGSR (“The International Genome Sample Resource”, https://www.internationalgenome.org/data-portal/population). El IGSR se estableció con el fin de garantizar el acceso y uso de los datos del Proyecto 1000 Genomas, actualizar los recursos con el ensamblaje de referencia actual, y ampliar el conjunto de datos con nuevos datos procedentes de las muestras del Proyecto 1000 Genomas y de nuevas poblaciones.
Debido a que disponemos de un tiempo limitado, utilizaremos una muestra de 20 o 50 individuos seleccionados de forma aleatoria por población (Tabla 1), y nos limitaremos a las variantes genéticas localizadas en el cromosoma 22.  

Tabla 1. Poblaciones humanas evaluadas por caso de estudio.

| Datos | Código| Población | Superpoblación / Ascendencia | Individuos |
| ----- | ----- | --------- | ---------------------------- | ---------- |
| Muestreo geográfico discreto (chr22_pop_dist)	 | YRI | Yoruba | África | 20 |
|  | LWK | Luhya | África | 20 |
|  | GBR | Británica | Europa | 20 |
|  | GIH | Gujarati | Asia meridional | 20 |
|  | CHB | Han Chinese | Asia Oriental | 20 |				
|  | PEL | Peruana | América | 20 |	
| Muestreo geográfico ‘continuo’ (chr22_pop_finer) | GBR | Británica | Europa | 50 |
|  | IBS | Ibérica | Europa | 50 |
|  | TSI | Toscana | Europa | 50 |			
|  | FIN | Finlandés | Europa | 50 |			
|  | GIH | Gujarati | Asia meridional | 50 |			


Empezaremos analizando la estructura poblacional en el conjunto de datos muestreados de forma discreta (**chr22_pop_dist**).


## Descarga de datos

Tras conectarte al CESGA, muévete al directorio de trabajo creado para esta sesión práctica.

```
cd $LUSTRE/MGG_GH5/STRUCTURE
```
Los scripts con los que trabajaremos están almacenados en este mismo repositorio. Para clonar el repositorio desde la terminal:

```
git clone https://github.com/noeliaperezp/MGG_GH5_STRUCTURE.git .
ls
```

Verás que se han descargado una serie de archivos:
- Scripts con extensión *.sh* listos para ejecutar y que incluyen cada uno de los comandos que iremos utilizando paso a paso a lo largo de la práctica.
- Un script de R (*structure_figures.R*) con el código que utilizaremos para representar los resultados.
- Un archivo README con el guion de prácticas para esta sesión.

Los archivos VCF (Variant Call Format) de cada cromosoma, con información para el total de los 2504 individuos secuenciados y pertenecientes a 26 poblaciones diferentes, pueden descargarse desde el sitio FTP:  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/. Desde la terminal, utilizaremos directamente el comando *wget* seguido de la URL del archivo que se quiere descargar. Para descargar el archivo VCF completo del cromosoma 22: 
```
mkdir data
cd data/
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```

Para descargar el panel con los identificadores de las poblaciones:
```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
cd ../
```

## Obtención de submuestras

```
mkdir data_subsamples
cd data_subsamples/
```

Empezaremos definiendo las variables que contendrán la información de las poblaciones a retener (Tabla 1) y el número de individuos por población.
```
POP_DIST=("YRI" "LWK" "GBR" "GIH" "CHB" "PEL")   # código de las poblaciones de interés
N_POP_DIST=20   # número de individuos por población
```

Haz una copia de la tabla con los identificadores de las poblaciones y listado de individuos:
```
awk 'NR>1{print $1"\t"$2"\t"$3}' ../data/integrated_call_samples_v3.20130502.ALL.panel > sample_pop.tsv
```

Crea una función que permita seleccionar individuos al azar a partir de la información contenida en la tabla anterior:
```
pick_rnd_ids() {
  local pop=$1   # id de la población
  local n=$2     # número de individuos a seleccionar
  awk -v p="$pop" '$2==p{print $1}' sample_pop.tsv | shuf -n "${n}"
} 
```

Obtén la lista de individuos:
```
for POP in "${POP_DIST[@]}"; do
  pick_rnd_ids "${POP}" "${N_POP_DIST}" | while read IID; do
    echo "${IID}" >> human_pop_dist.ids
    echo -e "0\t${IID}\t${POP}" >> poplist_human_pop_dist.txt
  done
done
```

Extrae de los VCF la submuestra:
```
module load cesga/2020 bcftools/1.19
bcftools view -S human_pop_dist.ids -Oz -o chr22_pop_dist.vcf.gz ../data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
bcftools index chr22_pop_dist.vcf.gz
cd ../
```

## Filtrado de datos

Carga los módulos necesarios:
```
module load cesga/2020 gcc/system plink/2.00a2.3
```

Crea un nuevo directorio y muévete a él:
```
mkdir data_pruned
cd data_pruned/
```

Retén variantes bialélicas (SNPs) y genera los archivos con formato plink (.bed,.fam,.bim):
```
plink2 --vcf ../data_subsamples/chr22_pop_dist.vcf.gz --make-bed --max-alleles 2 --snps-only --out chr22_pop_dist.SNPs
```

Asigna nuevos identificadores a variantes sin nombre:
```
plink2 --bfile chr22_pop_dist.SNPs --set-missing-var-ids @:# --make-bed --out chr22_pop_dist.SNPid
rm chr22_pop_dist.SNPs*
```

Filtra variantes y muestras con >10% missing calls y maf < 0.05:
```
plink2 --bfile chr22_pop_dist.SNPid --geno 0.1 --mind 0.1 --maf 0.05 --make-bed --out chr22_pop_dist_pruned
rm chr22_pop_dist.SNPid*
```

Filtra variantes con fuerte desequilibrio de ligamiento (LD):
```
# NOTA: como algunos análisis se basan precisamente en el LD, crearemos una copia.
plink2 --bfile chr22_pop_dist_pruned --indep-pairwise 50 5 0.9 --out chr22_pop_dist_filtLD.plink2
plink2 --bfile chr22_pop_dist_pruned --extract chr22_pop_dist_filtLD.plink2.prune.in --make-bed --out chr22_pop_dist_pruned_LD
rm chr22_pop_dist_filtLD*
rm *.in
rm *.log
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/STRUCTURE”):
```
cd ..
pwd
```

De aquí en adelante trabajaremos con los archivos filtrados (*_pruned_LD*), que encontrarás en el directorio *data_pruned/*. 

> **_PREGUNTA:_** ¿Cuántas muestras y variantes se han eliminado y cuántos retenemos? ¿son suficientes?


## Análisis de componentes principales (PCA)

Crea un nuevo directorio y muévete a él:
```
if [ ! -d results_pca/chr22_pop_dist_pruned_LD ] 
then 
 mkdir -p results_pca/chr22_pop_dist_pruned_LD
fi
cd results_pca/chr22_pop_dist_pruned_LD
```

PCA:
```
plink2 --bfile ../../data_pruned/chr22_pop_dist_pruned_LD --pca 10 -out chr22_pop_dist_pruned_LD_pca 
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/STRUCTURE”):
```
cd ../..
pwd 
```

En el directorio *results_pca/* encontrarás dos outputs, *.eigenval*, con la cantidad de varianza que explica cada componente principal (o eigenvector) en los datos, y *.eigenvec*, con los vectores que definen la dirección de los componentes principales.

### *Transferencia de ficheros del CESGA a tu sistema local*

Representaremos los resultados del PCA en R, y en concreto RStudio. Para ello necesitarás transferir los resultados del PCA a tu sistema local. También necesitarás un archivo adicional llamado *poplist_human_pop_dist.txt*, que contiene una lista de los individuos analizados junto con el código de su población de origen. Para transferir estos archivos puedes usar el programa WinSCP, o alternativamente usando el comando *scp* desde tu máquina local: 

En tu sistema local, crea y abre la carpeta a la que quieras exportar los archivos. Haz clic con el botón derecho y abre una terminal (ej. Windows PowerShell). 

Ten en cuenta que necesitarás conocer la ruta completa dentro del CESGA que conduce al archivo que quieres transferir. Para ello ejecuta el comando *pwd* en la terminal del CESGA. Asegúrate de que te encuentras en el directorio */MGG_GH5/STRUCTURE*, punto base de la estructura de directorios que hemos creado para esta sesión práctica.  

Ahora escribe en la terminal local (ej. Windows PowerShell) los siguientes comandos. Te pedirá confirmar la autenticidad del host, y a continuación tu contraseña de usuario en el CESGA.

Antes de copiar el siguiente código, recuerda cambiar:
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

Presionando tabulador tendrás una guía de los directorios a los que te puedes mover. Una vez en la carpeta copia y ejecuta las líneas de código que encontrarás en el apartado "PCA" en el script *structure_figures.R*.

> **_PREGUNTA:_** ¿Qué porcentaje de la varianza es explicada por PC1 y PC2 juntos? ¿Qué separa cada PC?


## Análisis de mezcla con ADMIXTURE

Carga los módulos necesarios:
```
module load cesga/2020 admixture/1.3.0
```

Crea un nuevo directorio y muévete a él:
```
if [ ! -d results_admixture/chr22_pop_dist_pruned_LD ] 
then 
 mkdir -p results_admixture/chr22_pop_dist_pruned_LD
fi
cd results_admixture/chr22_pop_dist_pruned_LD
```

Ejecuta el análisis con ADMIXTURE. Por lo general se suelen evaluar múltiples valores de K, desde 1 hasta el nº máximo de poblaciones que tenemos (*NPOP*), y se selecciona aquel valor K que minimiza el error (*cross-validation error (cv)*). 
```
NPOP=6   # cambiar en función de las poblaciones que se estén analizando
for (( K = 1; K <= $NPOP; K++ ))
do
  admixture --cv ../../data_pruned/chr22_pop_dist_pruned_LD.bed $K | tee chr22_pop_dist_pruned_admixture_K${K}.out
done
```

Extrae el error *cv* para cada valor de K:
```
awk '/CV/ {print $3,$4}' *.out | cut -c 4,7-20 > chr22_pop_dist_pruned.admixture.cv.error
```

Retorna al directorio principal (“*$LUSTRE/MGG_GH5/STRUCTURE”):
```
cd ../..
pwd 
```

En el directorio *results_admixture/* encontrarás tres outputs, *.Q*, que contiene las fracciones de ascendencia, *.P*, con las fecuencias alélicas de las poblaciones ancestrales inferidas, y *.cv.error* con los errores para cada valor de K.

### *Transferencia de ficheros del CESGA a tu sistema local*

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
if [ ! -d results_fst/chr22_pop_dist_pruned_LD ] 
then 
 mkdir -p results_fst/chr22_pop_dist_pruned_LD
fi
cd results_fst/chr22_pop_dist_pruned_LD
```

Calcula el F<sub>ST global:
```
plink --bfile ../../data_pruned/chr22_pop_dist_pruned_LD --fst --within ../../data_subsamples/poplist_human_pop_dist.txt --out chr22_pop_dist_pruned_fst
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/STRUCTURE”):
```
cd ../..
pwd 
```

En el directorio *results_fst/* encontrarás un output principal, *.fst*, que contiene las estimas por marcador, y un archivo *.log* que contiene el valor F<sub>ST promedio.

> **_PREGUNTA:_** ¿Qué nivel de diferenciación genética encontramos entre las poblaciones evaluadas?


## (TAREA) Análisis de la estructura poblacional en el conjunto de datos muestreados de forma 'continua' (**chr22_pop_finer**)

Siguiendo los mismos pasos, analiza la estructura poblacional del conjunto de datos muestreados de forma 'continua', o con mayor resolución espacial (**chr22_pop_finer**). Ten en cuenta que, al lanzar los comandos, deberás adaptar los nombres de las variables y de los archivos para trabajar con este conjunto de datos. 

Como alternativa, tras obtener la submuestra, puedes lanzar los scripts 'script_data_filtering.sh','script_run_pca.sh', 'script_run_admixture.sh', y 'script_calculate_fst.sh' con el comando *sbatch*. En este caso, revisa la lista de argumentos requeridos por el script, indicados en el apartado '# Set arguments'. Puedes visualizar el script con el comando *cat* <nombre_del_script>. Una vez enviados los trabajos, puedes revisar su estado mediante el comando *squeue*.

> **_PREGUNTA:_** ¿Qué porcentaje de la varianza es explicada por PC1 y PC2 juntos? ¿Qué separa cada PC?

> **_PREGUNTA:_** ¿Qué conclusiones se sacan del ADMIXTURE?

> **_PREGUNTA:_** ¿Qué nivel de diferenciación genética encontramos ahora entre las poblaciones evaluadas?


# ESTIMACIÓN DEL CENSO EFECTIVO

Estimaremos el censo efectivo de la historia reciente de una población humana de referencia (FIN, filandesa), utilizando el modelo basado en desequilibrio de ligamiento (LD) entre pares de SNPs implementado en el programa GONE2. 

Fundamento del método: cada nueva mutación genera un pequeño LD entre la mutación y otros sitios polimórficos. Este LD aumenta por deriva genética en cada generación, en una cantidad aproximadamente proporcional a 1/*N<sub>e*, mientras que la recombinación lo destruye, más rápido cuanto mayor es la distancia genética entre los loci. Por tanto, el LD entre loci distantes refleja principalmente la deriva reciente y los muy ligados la de un periodo más largo. El programa calcula la contribución de cada generación pasada al LD actual y busca la trayectoria de *N<sub>e* histórico que mejor reproduce el LD observado a todas las distancias genéticas (Santiago et al., 2020. *Mol. Biol. Evol.*).


Clona el repositorio de github de GONE2: 
```
pwd    # asegúrate de que te encuentras en STRUCTURE/
git clone https://github.com/esrud/GONE2
```

Compila el programa:
```
cd GONE2
make gone
cd ../
```

Crea un nuevo directorio y guarda el código de la población en una variable:
```
mkdir effective_size
cd effective_size/
POP="FIN"
```

Extrae los identificadores de las muestras FIN:
```
awk -v p="$POP" '$2==p {print "0", $1}' ../data_subsamples/sample_pop.tsv > "${POP}.ids"
```

Extrae la submuestra de la población FIN de los archivos en formato PLINK ya filtrados (incluido en **chr22_pop_finer**):
```
plink2 --bfile ../data_pruned/chr22_pop_finer_pruned --keep "${POP}.ids" --make-bed --export vcf --out "chr22_${POP}"
```

Lanza GONE2:
```
sbatch ../script_run_GONE2.sh chr22_FIN.vcf
```

