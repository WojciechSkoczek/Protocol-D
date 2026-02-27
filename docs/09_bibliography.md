# Protocol D — Bibliography and data-source citations

This repository uses multiple public catalogues and data releases. Please cite the relevant survey/data papers when using Protocol D outputs.

## Data catalogues used as raw inputs (primary citations)

### LoTSS (LOFAR Two-metre Sky Survey)
- **LoTSS DR3** (survey + catalogues): Shimwell et al., *The LOFAR Two-metre Sky Survey: VII. Third Data Release*, arXiv:2602.15949 (2026).  
  Data release portal: LOFAR Surveys “DR3: public data release”.
- **LoTSS DR2** (survey + catalogues): Shimwell et al., *The LOFAR Two-metre Sky Survey — V. Second data release*, A&A **659**, A1 (2022). DOI: 10.1051/0004-6361/202142484. arXiv:2202.11733.

### NVSS
- Condon et al., *The NRAO VLA Sky Survey*, AJ **115**, 1693–1716 (1998). DOI: 10.1086/300337.

### RACS-low (ASKAP)
- Hale et al., *The Rapid ASKAP Continuum Survey Paper II: First Stokes I Source Catalogue Data Release*, PASA **38**, e058 (2021). DOI: 10.1017/pasa.2021.47. arXiv:2109.00956.

### CatWISE2020 (WISE/NEOWISE)
- Marocco et al., *The CatWISE2020 Catalog*, ApJS **253**, 8 (2021). DOI: 10.3847/1538-4365/abd805. arXiv:2012.13084.

### 2MRS (2MASS Redshift Survey)
- Huchra et al., *The 2MASS Redshift Survey—Description and Data Release*, ApJS **199**, 26 (2012). DOI: 10.1088/0067-0049/199/2/26. arXiv:1108.0669.

### Quaia (Gaia–unWISE quasar catalogue)
- Storey-Fisher et al., *Quaia, the Gaia-unWISE Quasar Catalog: An All-sky Spectroscopic Quasar Sample*, ApJ **964**, 69 (2024). DOI: 10.3847/1538-4357/ad1328.  
  Data release: Zenodo “Quaia” record(s), e.g. 10.5281/zenodo.10403370.

## Dipole-tension / number-count dipole references (analysis context)

- Böhme et al., *Overdispersed Radio Source Counts and Excess Radio Dipole Detection*, Phys. Rev. Lett. **135**, 201001 (2025).  
  arXiv:2509.16732. Related DOI: 10.1103/6z32-3zf4.
- Secrest et al., *A Test of the Cosmological Principle with Quasars*, ApJL **908**, L51 (2021). DOI: 10.3847/2041-8213/abdd40. arXiv:2009.14826.
- Dam, Lewis & Brewer, *Testing the cosmological principle with CatWISE quasars: a bayesian analysis of the number-count dipole*, MNRAS **525**, 231–245 (2023). DOI: 10.1093/mnras/stad2322.

## Notes
- **Software citation** is handled separately via `CITATION.cff` (cite Protocol D itself + Zenodo DOI).
- Raw FITS/large catalogue files should remain untracked (`data/raw/*.fits` in `.gitignore`).