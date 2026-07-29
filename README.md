
## Figure 2 trajectory analysis

Run the main script once to prepare the Monocle3 trajectory:

```bash
Rscript R/09_Figure2.R
```

Save the selected root nodes (`Y_3` and `Y_17`):

```bash
Rscript R/09_Figure2_select_root_static.R Y_3,Y_17
```

Run the main script again to calculate pseudotime and generate the final Figure 2 outputs:

```bash
Rscript R/09_Figure2.R
```

To inspect the principal graph and its node labels, the helper can optionally be run without arguments:

```bash
Rscript R/09_Figure2_select_root_static.R
```

Figures are saved in `results/plots/Figure2/`, and source-data tables are saved in `results/tables/Figure2/`.

The trajectory is constructed using fibroblast and mesothelial populations. Mesothelial cells are retained during graph construction, root-node selection, and pseudotime calculation, and are excluded only from the final fibroblast-specific trajectory plots and downstream analyses.

