
For the most recently compiled PDF version of the manuscript, go to
<https://uwmadison.box.com/s/bso8x3z16x2hvybho6296njf2tqxwsgi>



Citations in `__ms.tex`:

1. Go to Google Scholar; click the top left, and select Settings; then in the
   Bibliography manager section, check show links to import citations into BibTeX, 
   then Save it. This step only need to be done once.
2. Back to the Google Scholar homepage and search the paper you want to cite. 
   E.g. "Generalized linear mixed models for phylogenetic analyses of community 
   structure", under the item, click Import into BibTeX; then copy the whole new page 
   (Cmd + A then Cmd + C).

```
@article{ives2011generalized,
  title={Generalized linear mixed models for phylogenetic analyses of community structure},
  author={Ives, Anthony R and Helmus, Matthew R},
  journal={Ecological Monographs},
  volume={81},
  number={3},
  pages={511--525},
  year={2011},
  publisher={Wiley Online Library}
}
```

3. Open the `refs.bib` in this folder, and paste the citation information there.

4. To cite the paper in the manuscript, copy its key (e.g. `ives2011generalized` here),
   and the LaTeX cite command (`\citep{ives2011generalized}` here) to where you
   want to cite it. For multiple papers, use `\citep{key1,key2}`.
   To add prefixes or suffixes, use `\citep[e.g.,][, see?]{ives2011generalized}` to
   generate "(e.g., Ives & Helmus 2011, see?)".

