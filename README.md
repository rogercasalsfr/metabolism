# Expression-informed metabolic modelling

Context-specific genome-scale metabolic models built from bulk RNA-seq, constrained with
measured culture-medium and Seahorse data, to compare metabolic states between treatments
and controls. Developed during a research stay at the Biomedical Center, University of Iceland (2025).

## Workflow

1. **Quantification and counts** (`counts/`, R)
   - `kallisto_to_counts.R`: kallisto transcript estimates to gene-level counts and TPM with tximport.
   - `combat-seq.R`: batch correction of the counts with ComBat-seq (sva).
   - `count_to_tpm.R`: corrected counts to TPM (IOBR).
2. **Model reconstruction** (`fastcormics/`, MATLAB)
   - `fastcormics_in.m`: builds one context-specific model per condition from the TPM matrix with
     [rFASTCORMICS](https://github.com/sysbiolux/rFASTCORMICS), using the COBRA Toolbox and the CPLEX solver.
3. **Transcriptomic analysis** (`anl/rna-seq/`, Python)
   - `rna_anl.ipynb`: differential expression and transcription-factor activity inference with
     [decoupler](https://decoupler-py.readthedocs.io) and the CollecTRI network.
   - `dio_tfs.ipynb`: regulators of selected target genes.
4. **Model analysis** (`anl/metabolism/`, Python)
   - `model_constrains-median.ipynb`: applies medium and Seahorse-derived exchange constraints and runs
     parsimonious FBA with [COBRApy](https://opencobra.github.io/cobrapy/), including glucose carbon allocation.
   - `median_mix.ipynb`, `median_zero.ipynb`: flux comparisons between conditions.
   - `metabolites_fastcormics.ipynb`: extracts the medium definition from the rFASTCORMICS `.mat` files.

## Requirements

- **R:** DESeq2, tximport, biomaRt, sva (ComBat-seq), IOBR
- **MATLAB** (R2019b), [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/), IBM ILOG CPLEX, rFASTCORMICS
- **Python:** cobra, decoupler (≥ 2.0), pandas, scikit-learn, h5py, matplotlib

Input and output paths are set at the top of each script and notebook; adjust them to your data.
The RNA-seq data are not included in this repository.

## License

Code in this repository is released under the [MIT License](LICENSE), except
`fastcormics/fastcormics_in.m`, which is adapted from the rFASTCORMICS example script
(© M. Pires Pacheco, T. Bintener) and remains under rFASTCORMICS'
[LGPL-3.0 licence](https://github.com/sysbiolux/rFASTCORMICS/blob/master/LICENSE).
