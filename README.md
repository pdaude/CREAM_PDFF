# Comparative REview of Algorithm and Method for proton-density fat fraction (PDFF) quantification (CREAM_PDFF)

The  CREAM_PDFF is a multi-language toolbox developed to evaluate state-of-the art open source algorithms for fat-water separation. 

| Reference          | Referred as   | Method       | Code   | 2D/3D | Spectrum choice   | Echo spacing   | Code acceleration  | Year | Code repository    *(If not specified, it is on github.com)*                                                              |
|--------------------|---------------|--------------|--------|-------|-------------------|----------------|--------------------|------|----------------------------------------------------------------------------------|
| Boehm et al. <sup>8</sup>     | VLGCA         | Graph-cut    | Matlab | 2D    | Free              | Uniform        | Parallel computing | 2021 | gitlab.com/christofboehm/fieldmap-graph-cut                                      |
| Bydder et al. <sup>7</sup>     | IDEAL-CE      | IDEAL        | Matlab | 3D    | Model Constrained | Free           | GPU computing      | 2020 | marcsous/pdff                                                                    |
| Andersson et al. <sup>6</sup>  | MSGCA-A       | Graph-cut    | Matlab | 3D    | Free              | Uniform        | No                 | 2018 | Snubben-B/FW-Recon-Spatial-Smoothing                                             |
| Cui et al.<sup>5</sup>         | GOOSE         | Graph-cut    | Matlab | 2D    | Fixed<sup>1</sup>             | Uniform / Free | No                 | 2015 |  [GOOSE code repository](https://research.engineering.uiowa.edu/cbig/content/goose)                        |
| Berglund et al.<sup>4</sup>    | MSGCA-B       | Graph-cut    | Python | 3D    | Free              | Uniform        | No                 | 2017 | bretglun/fwqpbo                                                                  |
| Liu˙ et al. <sup>3</sup>       | B0-NICE       | Region-based | Matlab | 3D    | Fixed<sup>2</sup>            | Uniform        | No                 | 2015 |  [B0-NICE code repository](https://www.mathworks.com/matlabcentral/fileexchange/48313-b0-mapping-b0-nice)    |
| Berglund et al.<sup>2</sup>   | Fatty-Riot-GC | Graph-cut    | Matlab | 3D    | Free              | Uniform        | No                 | 2012 | welcheb/fw_i3cm1i_3pluspoint_berglund_QPBO welcheb/FattyRiot                     |
| Hernando et al. <sup>1</sup>   | Hernando-GC   | Graph-cut    | Matlab | 2D    | Free              | Uniform / Free | No                 | 2012 | [Hernando-GC code repository](https://www.ismrm.org/workshops/FatWater12/data.htm) |

<sup>1</sup> ISMRM challenge spectrum

<sup>2</sup> own spectrum

## Publication 

This open-source is described in the following abstract :

- Daudé P, Kober F, Confort-Gouny S, Bernard M, Rapacchi S. Comparative Review of Algorithm and Method for proton-density fat fraction (PDFF) quantification. Proceedings [Internet]. Virtual Meeting; 2021. p. Poster. Available from: https://mritogether.github.io/files/abstracts/daude.pdf


## Installation guidelines

### VLGCA repository
To convert CUDA file to mex file

```bash
cd fieldmap-graph-cut/code/gandalf/
mexcuda residualcalculation_cuda.cu
```

## References


1. Hernando D, Kellman P, Haldar JP, Liang Z-P. Robust water/fat separation in the presence of large field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2009;NA-NA.

2. Berglund J, Kullberg J. Three-dimensional water/fat separation and T 2* estimation based on whole-image optimization-Application in breathhold liver imaging at 1.5 T. Magn Reson Med. 2012;67:1684–93.

3. Liu J, Drangova M. Method for B0 off-resonance mapping by non-iterative correction of phase-errors (B0-NICE): B0 Mapping with Multiecho Data. Magn Reson Med. 2015;74:1177–88.

4. Berglund J, Skorpil M. Multi-scale graph-cut algorithm for efficient water-fat separation: Multi-Scale Graph-Cut Water/Fat Separation. Magn Reson Med. 2017;78:941–9.

5. Cui C, Wu X, Newell JD, Jacob M. Fat water decomposition using globally optimal surface estimation (GOOSE) algorithm. Magn Reson Med. 2015;73:1289–99.

6. Andersson J, Ahlström H, Kullberg J. Water-fat separation incorporating spatial smoothing is robust to noise. Magn Reson Imaging. 2018;50:78–83.

7. Bydder M, Ghodrati V, Gao Y, Robson MD, Yang Y, Hu P. Constraints in estimating the proton density fat fraction. Magn Reson Imaging. 2020;66:1–8.

8. Boehm C, Diefenbach MN, Makowski MR, Karampinos DC. Improved body quantitative susceptibility mapping by using a variable‐layer single‐min‐cut graph‐cut for field‐mapping. Magn Reson Med. 2021;85:1697–712.











