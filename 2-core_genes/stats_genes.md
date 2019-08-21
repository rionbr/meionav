# Number of genes

|  Species   | Cell          |   # Genes |
| :----------|:--------------|----------:|
|  HS        | Cyte vs Gonia |     16165 |
|  HS        | Cyte vs Tid   |     13561 |
|  MM        | Cyte vs Gonia |     11806 |
|  MM        | Cyte vs Tid   |     13995 |
|  DM        | Apical Testis |     13700 |
|  DM        | Mid Testis    |     13700 |


# Number of genes differently expressed

|  Species   | Regulation   | Cell          | FDR    |   # Genes |
| :----------|:-------------|:--------------|:-------|----------:|
|  HS        | Up           | Cyte vs Gonia | <=0.01 |       625 |
|  HS        | Up           | Cyte vs Gonia | <=0.05 |      1440 |
|  HS        | Down         | Cyte vs Tid   | <=0.01 |      2338 |
|  HS        | Down         | Cyte vs Tid   | <=0.05 |      3133 |
|  MM        | Up           | Cyte vs Gonia | <=0.01 |      3847 |
|  MM        | Up           | Cyte vs Gonia | <=0.05 |      3848 |
|  MM        | Down         | Cyte vs Tid   | <=0.01 |      2693 |
|  MM        | Down         | Cyte vs Tid   | <=0.05 |      2694 |
|  DM        | Up           | Apical Testis | n/a    |      9683 |
|  DM        | Down         | Mid Testis    | n/a    |      9379 |


# Number of protein-coding genes differently expressed

## HS with FDR=0.01

|                      |   biotype |
|:---------------------|----------:|
| protein_coding       |      1823 |
| lncRNA               |         2 |
| processed_pseudogene |         1 |

## HS with: FDR=0.05

|                                    |   biotype |
|:-----------------------------------|----------:|
| protein_coding                     |      2733 |
| lncRNA                             |         4 |
| transcribed_unprocessed_pseudogene |         1 |
| processed_pseudogene               |         1 |
| transcribed_unitary_pseudogene     |         1 |

## MM with FDR=0.01

|                                    |   biotype |
|:-----------------------------------|----------:|
| protein_coding                     |      5102 |
| lncRNA                             |        20 |
| processed_pseudogene               |         7 |
| unprocessed_pseudogene             |         2 |
| transcribed_processed_pseudogene   |         2 |
| transcribed_unitary_pseudogene     |         2 |
| polymorphic_pseudogene             |         2 |
| TEC                                |         1 |
| transcribed_unprocessed_pseudogene |         1 |

## MM with FDR=0.05

|                                    |   biotype |
|:-----------------------------------|----------:|
| protein_coding                     |      5104 |
| lncRNA                             |        20 |
| processed_pseudogene               |         7 |
| unprocessed_pseudogene             |         2 |
| transcribed_processed_pseudogene   |         2 |
| transcribed_unitary_pseudogene     |         2 |
| polymorphic_pseudogene             |         2 |
| TEC                                |         1 |
| transcribed_unprocessed_pseudogene |         1 |

## DS

|                |   biotype |
|:---------------|----------:|
| protein_coding |     10446 |
| pseudogene     |        14 |
| ncRNA          |         8 |
