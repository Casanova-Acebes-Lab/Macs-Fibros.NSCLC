# AGENTS.md

# Macs-Fibros.NSCLC

Este repositorio contiene el código reproducible de análisis de single-cell
RNA-seq para el proyecto Macs-Fibros.NSCLC.

El objetivo principal al trabajar en este repositorio es preservar la
reproducibilidad científica y los análisis históricos validados, realizando
cambios pequeños, explícitos y fáciles de revisar.

## Idioma

Comunícate con el usuario en español salvo que pida explícitamente otro idioma.

El código, nombres de variables y comentarios pueden conservar el estilo
existente del repositorio, que mezcla inglés y español.

## Principio fundamental

La corrección científica y la reproducibilidad tienen prioridad sobre la
refactorización o modernización del código.

No cambies decisiones biológicas, parámetros, thresholds, anotaciones,
selección de muestras, identidades de clusters, dimensionalidad, métodos
estadísticos o definiciones de figuras salvo que el usuario lo solicite
explícitamente.

Si un cambio solicitado puede alterar un resultado científico, indícalo
claramente antes de implementarlo.

No sustituyas una decisión histórica por otra simplemente porque parezca más
moderna o convencional.

## Estado actual del repositorio

Este repositorio puede contener trabajo local no versionado o parcialmente
modificado.

En el momento de crear estas instrucciones existían:

- `R/00_Doublets.R` eliminado localmente;
- `R/01_QC.R` modificado;
- varios scripts de análisis y figuras no versionados.

Estos cambios forman parte del trabajo en curso.

Nunca:

- restaures automáticamente archivos eliminados;
- descartes modificaciones locales;
- sobrescribas archivos no versionados;
- ejecutes `git clean`;
- ejecutes `git reset --hard`;
- intentes sincronizar el árbol de trabajo con GitHub por iniciativa propia.

No interpretes un estado Git "sucio" como un problema que deba resolverse.

## Git

Antes de una modificación no trivial, ejecuta:
```bash
git status --short
```

Los comandos de inspección Git son seguros, por ejemplo:

```bash
git status
git status --short
git diff
git diff --stat
git log
git symbolic-ref --short HEAD
```

La versión de Git disponible en alguno de los sistemas utilizados por el
proyecto puede no soportar:

```bash
git branch --show-current
```

En ese caso usa:

```bash
git symbolic-ref --short HEAD
```

No ejecutes sin petición explícita del usuario:

```bash
git add
git commit
git push
git pull
git merge
git rebase
git reset
git checkout
git switch
git clean
```

No crees ni elimines ramas salvo petición explícita.

Nunca uses:

```bash
git add .
git add -A
```

Para cualquier futura operación de staging, añade únicamente rutas
específicamente aprobadas por el usuario.

## Estructura general del pipeline

Los scripts principales se encuentran en `R/`.

El pipeline actual incluye:

- `01_QC.R`: control de calidad, preprocesamiento y clustering global.
- `02_Subcluster.Macros.R`: subclustering de poblaciones mieloides/macrófagos.
- `03_Subcluster.Fibros.R`: subclustering de fibroblastos y mesoteliales.
- `04_CD45neg.R`: análisis/anotación de poblaciones CD45 negativas.
- `05_CancerCells.R`: análisis de células tumorales.
- `06_Tcells.R`: análisis y reclustering de células T.
- `07_Join.Clustering.R`: integración final de anotaciones.
- `08_Figure1.R`: Figure 1.
- `09_Figure2.R`: Figure 2 y análisis de trayectoria/pseudotiempo.
- `09_Figure2_select_root_static.R`: selección manual de raíces de pseudotiempo.
- `10_Figure3.R`: Figure 3.
- `11_Figure4.R`: Figure 4.
- `Functions.R`: funciones auxiliares compartidas.

No asumas que todos estos archivos están ya versionados en Git.

## Directorio de ejecución

Los scripts están diseñados para ejecutarse desde la raíz del repositorio.

Antes de ejecutar un script, comprueba que el directorio actual corresponde a
la raíz del repositorio.

No cambies silenciosamente el sistema de rutas para permitir ejecución desde
otros directorios.

## Configuración de rutas

`config.tsv` define actualmente:

```text
rsdir   = results
datadir = data
outdir  = results/plots
```

Los scripts utilizan rutas relativas al repositorio y `file.path()`.

Prefiere siempre las rutas definidas en `config.tsv`.

No introduzcas rutas absolutas específicas de una máquina salvo que el usuario
lo solicite.

## Entornos de computación

El proyecto se trabaja principalmente desde:

- `cluster1`
- `gurugu`

No asumas que:

- las rutas absolutas son idénticas;
- las versiones de R son idénticas;
- los paquetes instalados son idénticos;
- los entornos Conda son idénticos;
- los recursos computacionales son idénticos.

La información de reproducibilidad y paquetes se conserva en `env/`.

No modifiques archivos de `env/` salvo petición explícita.

## Seguridad en el cluster

`cluster1` es un entorno HPC compartido.

No ejecutes análisis pesados en un nodo login/head sin autorización explícita.

Antes de ejecutar un comando, distingue entre:

1. inspección ligera;
2. prueba pequeña;
3. análisis real o potencialmente costoso.

Puedes ejecutar inspecciones ligeras cuando sean necesarias.

Para análisis que puedan requerir CPU, RAM o tiempo significativo, pide
permiso antes de ejecutarlos.

No envíes trabajos SLURM automáticamente.

No solicites recursos mediante `srun`, `sbatch` u otros sistemas de colas salvo
que el usuario lo pida.

## Entornos y paquetes

No instales, actualices ni elimines sin autorización explícita:

- paquetes de R;
- paquetes Conda;
- paquetes Python;
- módulos;
- dependencias del sistema.

No modifiques:

- `.bashrc`;
- `.profile`;
- configuración global de R;
- configuración global de Conda;
- módulos del sistema.

Si falta una dependencia, informa primero del problema.

## Datos

Trata `data/` como datos de entrada.

No modifiques, regeneres, renombres ni elimines datos de entrada salvo petición
explícita.

Considera especialmente como entradas potencialmente sensibles:

- objetos `.rds`;
- matrices;
- barcodes;
- firmas;
- HVGs;
- coordenadas UMAP históricas;
- referencias;
- datos raw.

No sobrescribas objetos existentes únicamente para probar código.

## Resultados

Los resultados se almacenan principalmente bajo:

```text
results/
├── logs/
├── objects/
├── plots/
├── qc/
└── tables/
```

No sobrescribas resultados importantes existentes salvo que el workflow
solicitado lo requiera claramente.

Cuando hagas una prueba, utiliza cuando sea posible una salida temporal o un
nombre claramente identificable como test.

No añadas grandes resultados generados a Git.

## Reproducibilidad histórica

Este proyecto contiene análisis diseñados para reproducir resultados
previamente validados.

Preserva cuando aparezcan en los scripts:

- número de PCs;
- resoluciones de clustering;
- semillas aleatorias;
- HVGs;
- PCA features;
- coordenadas UMAP históricas;
- selección de células;
- cluster mappings;
- niveles y orden de factores;
- colores;
- thresholds;
- roots de pseudotiempo;
- definiciones de marcadores;
- comparaciones estadísticas.

No cambies estos parámetros para "mejorar" el análisis sin una petición
explícita.

Si encuentras algo que parece inconsistente, descríbelo antes de cambiarlo.

## Seurat

Este proyecto utiliza Seurat.

Antes de modificar código relacionado con un objeto Seurat, determina cuando
sea relevante:

- versión de Seurat;
- assays;
- layers;
- reductions;
- metadata;
- identidades;
- objeto de entrada del workflow.

No asumas semántica de Seurat v4 si el análisis está usando Seurat v5.

No cambies silenciosamente:

- `DefaultAssay`;
- layer;
- normalización;
- variable features;
- variables regresadas;
- PCA;
- número de dimensiones;
- vecinos;
- resolución;
- UMAP;
- identities.

En análisis diferenciales preserva explícitamente:

- población comparada;
- definición de identities;
- dirección de la comparación;
- thresholds;
- corrección por testing múltiple.

## Objetos y anotaciones

Antes de utilizar un objeto de análisis:

1. identifica qué script lo generó;
2. revisa sus metadata relevantes;
3. no asumas nombres de columnas;
4. no reconstruyas una etapa upstream usando información downstream.

Evita cualquier lógica circular.

Por ejemplo, un workflow oficial de clustering no debe utilizar las
anotaciones finales para recrear los clusters que pretende generar.

Las anotaciones celulares son resultados científicos.

No renombres, fusiones, dividas o reasignes poblaciones según tu propia
interpretación salvo petición explícita.

Puedes proponer una interpretación, pero no implementarla sin decisión del
usuario.

## Functions.R

Antes de modificar funciones compartidas, determina qué scripts dependen de
ellas.

Ten en cuenta dos particularidades observadas:

- `read_sample_remove_doublets()` depende de la variable global `datadir`;
- `Functions.R` define una función llamada `scale()`, que puede ocultar
  `base::scale()`.

No "corrijas" automáticamente estas características.

Un cambio en `Functions.R` puede afectar múltiples etapas del pipeline y debe
tratarse como un cambio de alto impacto.

## Figuras

Los scripts de figuras forman parte del workflow de publicación.

No cambies silenciosamente:

- colores;
- nombres de poblaciones;
- orden de factores;
- orden de facetas;
- etiquetas;
- leyendas;
- escalas;
- dimensiones de figuras;
- selección de pathways;
- selección de comparaciones;
- anotaciones estadísticas.

Preserva la identidad visual existente salvo que el usuario pida un cambio.

Las figuras se guardan principalmente bajo:

```text
results/plots/
```

Cuando el análisis de una figura genera datos utilizados directamente para el
panel, exporta también una tabla de source data cuando forme parte del workflow.

Las tablas se almacenan normalmente bajo:

```text
results/tables/
```

Cuando añadas un panel nuevo a una figura existente:

1. no alteres los paneles no relacionados;
2. separa claramente el código por panel;
3. utiliza las mismas definiciones biológicas que el resto de la figura;
4. exporta sus datos fuente cuando corresponda.

## Trayectorias y pseudotiempo

Figure 2 contiene análisis con Monocle y decisiones manuales o biológicamente
informadas.

No selecciones automáticamente roots de pseudotiempo.

No cambies sin petición explícita:

- células utilizadas para construir la trayectoria;
- grafo;
- roots;
- definición de pseudotiempo;
- inclusión/exclusión final de poblaciones.

El script `09_Figure2_select_root_static.R` existe específicamente para apoyar
la selección manual de roots.

## CellChat y comunicación celular

Al trabajar con CellChat u otros análisis de comunicación celular no cambies
sin petición explícita:

- poblaciones;
- bases de datos ligando-receptor;
- sender/receiver;
- thresholds;
- probabilidades;
- filtros;
- pathways seleccionados.

Distingue siempre entre:

1. cambiar únicamente lo que se visualiza;
2. cambiar el análisis o filtrado subyacente.

No confundas ambos tipos de cambio.

## Estilo de R

Conserva el estilo existente del repositorio.

Prefiere código explícito y legible.

Usa cuando corresponda:

```r
suppressPackageStartupMessages({
  library(...)
})
```

Mantén los separadores de sección del estilo:

```r
# ============================================================
# Section
# ============================================================
```

Prefiere:

```r
file.path(...)
```

en lugar de concatenar manualmente rutas.

Crea directorios de forma explícita:

```r
dir.create(
  path,
  recursive = TRUE,
  showWarnings = FALSE
)
```

No añadas nuevas dependencias si el trabajo puede realizarse razonablemente
con las dependencias existentes.

No cambies seeds de manera arbitraria.

Cuando exista aleatoriedad relevante, utiliza una seed explícita.

## Depuración

Cuando aparezca un error:

1. localiza la operación que falla;
2. inspecciona el objeto o variable relevante;
3. realiza primero una prueba pequeña;
4. aplica el cambio mínimo necesario;
5. evita reescribir el script completo si no es necesario.

No ejecutes un pipeline costoso para diagnosticar algo que pueda comprobarse
con objetos ya existentes.

Si el fix cambia el comportamiento científico, indícalo.

## Alcance de los cambios

Si el usuario pide modificar un script, figura, panel o función, limita los
cambios a ese ámbito siempre que sea posible.

No hagas simultáneamente:

- limpieza estética;
- renombrados;
- refactorizaciones;
- actualizaciones de paquetes;
- cambios de estilo;

que no sean necesarios para resolver la petición.

Si necesitas modificar varios archivos, explica por qué.

## Antes de editar

Para una modificación no trivial:

1. ejecuta `git status --short`;
2. lee completamente el script que vas a modificar;
3. revisa los scripts upstream/downstream directamente relacionados si es
   necesario;
4. identifica los objetos de entrada y salida;
5. resume brevemente el cambio que vas a realizar;
6. implementa el cambio mínimo.

No hagas preguntas innecesarias si la petición es clara y de bajo riesgo.

## Después de editar

Después de modificar código:

1. indica qué archivos cambiaste;
2. resume exactamente qué cambió;
3. indica qué comprobaciones realizaste;
4. explica qué no se pudo comprobar;
5. muestra `git status --short`;
6. muestra un resumen de `git diff` cuando resulte útil.

No hagas commit, push ni staging automáticamente.

## Ejecución de scripts

No ejecutes automáticamente scripts completos de análisis sólo porque los hayas
modificado.

Después de una modificación puedes realizar comprobaciones ligeras de sintaxis
o inspección cuando sean seguras.

Antes de ejecutar un análisis completo, pide autorización si puede ser
computacionalmente costoso o modificar resultados existentes.

## Regla final

El usuario toma las decisiones científicas finales.

Tu función es:

- inspeccionar;
- implementar con precisión;
- detectar inconsistencias;
- depurar;
- explicar consecuencias;
- preservar la reproducibilidad.

No cambies silenciosamente el significado científico del análisis.
