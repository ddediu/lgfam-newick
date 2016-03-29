# lgfam-newick: Language family classifications as Newick trees

## Summary

This repository contains the data, R code, outputs and description of a flexible method for generating standardized [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) language family trees with branch lengths from the four most used language classification databases: [Ethnologue](http://www.ethnologue.com/), [WALS](http://wals.info/), [AUTOTYP](http://www.autotyp.uzh.ch/) and [Glottolog](http://glottolog.org/).
The code is released under [GPL v2](http://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html), but the various pieces of input data might be governed by different licenses (specified in the respective folders).

The aims of this project are to:

a) provide several well-known linguistic (genealogical) classifications (currently [WALS](http://wals.info/), [Ethnologue](http://www.ethnologue.com/), [Glottolog](http://glottolog.org/) and [AUTOTYP](http://www.autotyp.uzh.ch/)) in the *de facto* standard [Newick format](https://en.wikipedia.org/wiki/Newick_format), and
b) offer a set of [`R`](http://www.r-project.org/) `S3` classes and functions for reading, converting, writing and working with language family trees.

Also included is code by [Balthasar Bickel](http://www.linguistik.uzh.ch/en/about/mitglieder/bickel.html) that matches tree nodes to datasets and prunes the trees to keep only the nodes that have matching data (the `./code/MatchTreesToData.R` script).s

## Accompanying paper, outputs and acknowledging this work

The **accompanying paper** (in the `./paper/` directory) describes in detail the data sources and the conversion process.
The paper itself is written in [`R Markdown`](http://rmarkdown.rstudio.com/) and can be compiled to PDF (the primary output in the `family-trees-with-brlength.pdf` file) or HTML (the `family-trees-with-brlength.html` file).

The actual Newick trees with branch lengths are in the `./output/` directory and can be used directly (the file formats are described in the **accompanying paper** but briefly they come as **CSV TAB-separated files** and equivalent **Nexus files** that contain the language family trees in the **Newick format**; the file name gives details about the classification, method and parameters used to compute the topology and branch lengths).

**Note**: when using these trees from `R` the best (and recommended) way to read them is with the function `languageclassification()` (in file `FamilyTrees.R`) which returns an `S3` object of type `languageclassification` containing the list of trees and giving access to various useful things such as pretty printing, collapsing and restoring single nodes, etc. (besides, those trees extend the standard `phylo` class so most usual things should work out-of-the-box). Definitely do *not* use `ape`'s `read.tree()` (as it is known to be pretty fussy especially when it comes to single nodes) and if you must please do use instead `phytools`'s `read.newick()` instead!

If you use (parts of) the `R` scripts and/or the generated Newick trees, please do cite this in your work and provide links to this repository ([https://github.com/ddediu/lgfam-newick](https://github.com/ddediu/lgfam-newick))!


## Releases

"Official" releases can be found in the `./relases` directory.


## Running the `R` code

If you are **trying to run the `R` code yourself**, please note that I have removed some of the large cached intermediary results (in order to save space).
Thus, you must first generate these cached data, as follows.

Run the `./input/distances/WALS/process-wals-distances.R` script to generate the WALS-based distance matrices.

Run the `./input/distances/ASJP/process-asjp16-distances.R` script to generate the ASJP16-based distance matrix.

Run the `./code/StandardizedTrees.R` main `R` script with the following parameters set to `TRUE`: `MATCH_CODES` (compute the equivalences between the ISO, WALS, AUTOTYP and GLOTTOLOG codes and generate the UULIDs), 'PREOPTIMIZE_DISTS' (pre-optimize the distance matrices for fast loading when required), `COMPUTE_GEO_DISTS` (compute the geographic distances between languages).
For later runs (after these data has been generated and cached) these parameters can be safely set to `FALSE` (this pre-processing is computationally very expensive).
The parameters `TRANSFORM_TREES` (transform the trees from their original specific representation to the Newick notation no branch length), `EXPORT_NEXUS` (export the trees to a NEXUS file), `EXPORT_NEXUS_TRANSLATE_BLOCK` (when exporting NEXUS files, generate a TRANSLATE block; useful when using programs such as BayesTraits that have issues parsing complicated taxa names), `EXPORT_CSV` (export the trees to a CSV file) can be left on `TRUE` (except perhaps the first as the tree topologies will probably not change very often in the original databases).
Please note that the first time the Ethnologue tree topologies are transformed to Newick, these will be downloaded from the Ethnologue website and cached locally.
The last two parameters are `COMPUTE_BRLEN` (apply the various branch length methods to the Newick topologies) and `COMPARE_TREES` (compute the distance between equivalent trees).
Finally, `CPU_CORES` controls multi-core processing (using `mclapply` -- might not work on Windows!).
(It is a good idea to leave `quotes="'"`).
Parameters `CLASSIFICATIONS`, `METHODS`, `CONSTANT` and `DISTS.CODES` control which classification, methods and parameters to use for generating the Newick trees.
These are very specific to the current implementation but can be used to extend this work to other classifications of branch length methods.s

## Possible bugs! Please report them!

Please note that even if the `R` code is relatively well-tested there might be bugs or other issues!
So, please use these with caution and any comments, suggestions or bug reports are welcome, either through GitHub's own issue reporting facilities or by e-mail to <Dan.Dediu@mpi.nl>. 


## Thank you

Dan Dediu

The Netherlands

October 2015


