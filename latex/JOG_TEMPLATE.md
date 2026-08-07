# Journal of Glaciology LaTeX files

The manuscript uses the International Glaciological Society review class and
bibliography style:

- `igs.cls` (v4.00, 3 September 2015)
- `igs.bst`
- `igsnatbib.sty`

The files are unchanged copies from commit
`a070bd9bd8b1ced825ff60f62eda05987030fc67` of
[`bueler/perf-model-ism`](https://github.com/bueler/perf-model-ism/tree/a070bd9bd8b1ced825ff60f62eda05987030fc67/paper),
the source archive for a published *Journal of Glaciology* paper. This version
matches the class used by the journal's
[official Overleaf template](https://www.overleaf.com/latex/templates/latex-template-for-journal-of-glaciology-jog/prcsfgbpbckc).
The Cambridge author-instructions link to a standalone class download returned
an error when these files were retrieved on 6 August 2026.

Build the review manuscript from this directory with

```sh
latexmk -pdf manuscript.tex
```
