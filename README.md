# Mobile Element Locator Tool - Long Read Asssembly (MELT-LRA)

MELT-LRA takes as its input a set of [PAV](https://github.com/EichlerLab/pav)-generated variant
calls (i.e., insertions or deletions) and produces as its output the subset of those variants
that appear to correspond with known Mobile Element Insertions (MEIs), annotated with the
evidence for that assertion.

## MELT-LRA output

MELT-LRA produces output in the following formats:

 - VCF
 - Ad-hoc CSV (comma separated value)
 - ASCII art

MELT-LRA ASCII art renditiion of a single MEI call, an Alu on chromosome 1:

```
Location        |ME   |+/-|%ME   |%id   |%cov  | insertion
chr1:710579     |ALU  |+  |100.0%| 97.9%| 86.5%| atactgctat [AAAGAACTGCCCGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCA....+234bp....CGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA] aaagaactgcccaagactgggtaatttata
                                                             <---TSD---->                                                             <------------polyA------------>  <---TSD---->
                                                                         [ALU-----------------------------              ----------ALU>                               
```


## MEI Viewer

MELT-LRA also includes an interactive MEI Viewer implemented using [VueJS](https://vuejs.org/).
The viewer reads and displays MELT-LRA output in both tabular and graphical form. The viewer
allows an analyst to adjust MELT-LRA parameter settings (within the range dictated by the original
pipeline parameters) and see the effect on the resulting MEI callset, providing a way to 
perform manual parameter tuning:

![MEI Viewer](docs/images/mei-viewer-1.png)

Further curation is enabled by the detailed MEI view, which presents a graphical summary of the
evidence used to generate each MEI call:

![Homozygous Alu on chr22](docs/images/mei-viewer-alu-1.png)

![Alu on chr1, reverse strand](docs/images/mei-viewer-alu-2.png)

![MEI signature based on reference diffs](docs/images/mei-viewer-calu-lineu-diffs-1.png)


