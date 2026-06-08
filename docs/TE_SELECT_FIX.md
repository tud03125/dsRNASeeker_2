# TE ChIPseeker select() conflict fix

This patch fixes an R namespace conflict where AnnotationDbi/S4 generics mask dplyr::select(), producing:

```
Error: unable to find an inherited method for function ‘select’ for signature ‘x = "data.frame"’
```

Replace `r/te_atena_chipseeker_to_dsRNASeeker_csv.R` with the patched version in this archive.
