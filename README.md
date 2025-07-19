Environmental Heavy Metals Alter Life History
=============================================

Anuran body size impacts the their reproduction, potentially manipulating their sustainability. However, their sustainability may be threatened by the wide occurrence of environmental heavy metals (HMs) in aquatic environment. Here we utilized data mining to analyze the global distribution patterns and the relationship between HMs bioaccumulation and body size. As a result, we find both a hot spot of small body size and a cold spot of high HMs bioaccumulation occur around 20°N, with low correlation to the species difference. According to this difference, the samples were divided into two groups that high bioaccumulation group appears smaller sizes while another one appears larger sizes. This difference may be attributed more to the evolution, not the plasticity. Here we conclude that the high HM loads have altered their sustaining strategies that they distribute more energy to detoxify with the price of shrinking body.<br>

Our repository consists of `Python Codes` and `Fundamental Data`.<br>
------------------------------------------------------------
The `Python Codes` include:<br>
>Preprocess.py<br>
>Body_Size.py<br>
>Heavy_Metal.py<br>
>Quantile.py<br>

The `Fundamental Data` include:<br>
>raw_data.xlsx<br>
>Species.nwk<br>
>Taxonomy.xlsx<br>
>Continent.xlsx<br>
>Chordplot.xlsx<br>
    
## `Preprocess.py`<br>
### Transformation
Before imputation, we have tested that all body size `BS` features, including `Snout-vent Length(SVL)`, `Body Mass(BM)` and `Condition Factor(CF)`, and bioaccumulated heavy metals are right-skewed. Therefore we applied logarithm transformation to normalize the data.<br>

### Imputation
Our original data `raw_data.xlsx` with `n` rows were collected from the previous articles. Given that there are many empty cells in this file, we used `Phylo-KNN` to impute the null cells.<br>
<br>
* KNN<br>
KNN clusters the current data according to the distance $d_{ij}$ between the sample pairs.<br>
Here $s_i$ and $s_j$ represent two samples. When $i = j$, both samples are the same (that is $d_{ij} = 0$).<br>
After this, we will get a `n×n` distance matrix $D_{ij}$<br>
* Phylo-KNN<br>
However, the way that KNN clusters samples may not capture the phylogenetic difference.<br>
Therefore we develop a `Phylo-KNN` to capture the phylogenetic features by adding a `n×n` phylogenetic distance matrix $P_{ij}$.<br>
Here $p_{ij}$ refers to the distance from the phylogenetic root to the common ancestor `MRCA` of species `i` and species `j`.<br>
After we acquire the `n×n` $D_{ij}$ and the `n×n` $P_{ij}$, a dot product of both matrix will be used as a revised sample distance matrix.<br>
<p align="center">$d_{ij} = (\sum |s_i - s_j|^p)^{\frac{1}{p}}$<p>
<p align="center">$D_{ij} = (S_{ij})$<p>
<p align="center">$P_{ij} = (p_{ij})$<p>

### Revision
Our research partly focuses on the HMs effects on the Anuran body size `BS`. Considering it is well proved that there are some covariance (including `latitude`, `altitude` and `sex`) impacting their body size, we should eliminate the covariance impacts with<br>
<p align="center">$log BS_{Revised} = log BS_{Original} - (β_1 × latitude + β_2 × altitude + β_3 × sex)$<p>

## `Body_Size.py`<br>
In part 1 of Results, we focus on the global latitudinal `BS` distribution patterns.<br>
### Scatter Plot
Firstly we illustrate scatter plots with x and y axes referring to `BS` and `latitude`.<br>
We assume that these distribution patterns should be symmetric along the equator. To capture the symmetric features, we use `Quadratic Regression` to fit the data derived from `Preprocess.py`.<br>

### World Map
Secondly, to visualize to geological distribution of `BS`, we plot the `BS` sactters on a world map.<br>

```
Intraspecific trait variance (ITV) indicates the trait variance (that is, trait plasticity).
We use this parameter to show the plasticity difference among different sample sites, shown in the supplementary materials.
```

## `Heavy_Metal.py`<br>
After reviewing, we only select seven `HMs` that are obserrved they change the body size, including `Cd`, `Cr`, `Cu`, `Fe`, `Hg`, `Mn`, `Pb` and `Zn`.<br>
Parallelly, we apply `Quadratic Rgression` to fit the latitudinal distribution patterns of the selected HMs and plot their correspongding geological distribution map.<br>

## `Quantile.py`<br>

