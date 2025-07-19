Environmental Heavy Metals Alter Life History
=============================================

Anuran body size impacts the their reproduction, potentially manipulating their sustainability. However, their sustainability may be threatened by the wide occurrence of environmental heavy metals (HMs) in aquatic environment. Here we utilized data mining to analyze the global distribution patterns and the relationship between HMs bioaccumulation and body size. As a result, we find both a hot spot of small body size and a cold spot of high HMs bioaccumulation occur around 20°N, with low correlation to the species difference. According to this difference, the samples were divided into two groups that high bioaccumulation group appears smaller sizes while another one appears larger sizes. This difference may be attributed more to the evolution, not the plasticity. Here we conclude that the high HM loads have altered their sustaining strategies that they distribute more energy to detoxify with the price of shrinking body.<br>

Our repository consists of `Python Codes` and `Fundamental Data`.<br>
------------------------------------------------------------
The `Python Codes` include:<br>
>Preprocess.py*<br>
>Body_Size.py*<br>
>Heavy_Metal.py*<br>
>Quantile.py*<br>

The `Fundamental Data` include:<br>
>raw_data.xlsx*<br>
>Species.nwk*<br>
>Taxonomy.xlsx*<br>
>Continent.xlsx*<br>
>Chordplot.xlsx*<br>
    
## `Preprocess.py`<br>
### Transformation
Before imputation, we have tested that all body size features, including `Snout-vent Length(SVL)`, `Body Mass(BM)` and `Condition Factor(CF)`, and bioaccumulated heavy metals are right-skewed. Therefore we applied logarithm transformation to normalize the data.<br>

### Imputation
Our original data `raw_data.xlsx` with `n` rows were collected from the previous articles. Given that there are many empty cells in this file, we used Phylo-KNN to impute the null cells.<br>
<br>
KNN clusters the current data by calculating the distance $D_{ij}$ between the sample pairs.<br>
<p align="center">$d_{ij} = (\sum |S_i - S_j|^p)^{\frac{1}{p}}$<p><br>
Here $S_i$ and $S_j$ represent two samples. When $i = j$, both samples are the same (that is $D_{ij} = 0$).<br>
After this, we will get a distance matrix $D$<br>
<p align="center">$D = (d_{ij})$<p><br>
However, the way that KNN clusters samples may not capture the phylogenetic difference.<br>
Therefore we develop a Phylo-KNN to capture the phylogenetic features by adding a `n × n` phylogenetic distance matrix $P$.<br>
<p align = "center">$P = (p_{ij})$<p><br>
Here $p_{ij}$ refers to the distance from the phylogenetic root to the common ancestor `MRCA` of species `i` and species `j`.<br>


### Revision
Our research partly focuses on the HMs effects on the Anuran body size. Considering it is well proved that there are some covariance (including latitude, altitude and sex) impacting their body size, <br>

## `Body_Size.py`<br>


## `Heavy_Metal.py`<br>


## `Quantile.py`<br>
