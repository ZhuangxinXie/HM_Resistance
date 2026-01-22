Anuran life history evolution converges under heavy metal pressure
=============================================

Anuran body size impacts the their reproduction, potentially manipulating their sustainability. However, their sustainability may be threatened by the wide occurrence of environmental heavy metals (HMs) in aquatic environment. Here we utilized data mining to analyze the global distribution patterns and the relationship between HMs bioaccumulation and body size. As a result, we find both a hot spot of small body size and a cold spot of high HMs bioaccumulation occur around 20°N, with low correlation to the species difference. According to this difference, the samples were divided into two groups that high bioaccumulation group appears smaller sizes while another one appears larger sizes. This difference may be attributed more to the evolution, not the plasticity. Here we conclude that the high HM loads have altered their sustaining strategies that they distribute more energy to detoxify with the price of shrinking body.<br>

Our repository consists of `GEE Codes`, `Python Codes` and `Fundamental Data`.<br>
------------------------------------------------------------
The `GEE Codes` include:<br>
>Extract_Altitude<br>
>Extract_Climate<br>

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

## `GEE Codes`<br>
Considering some covariances of the physical environment (such as `Alutitude`, `Latitude`, `Sex`, `Precipitation`, `Temperature` and `Evaporation`) impact the Anuran body size, here we should collect these related values of the sample sites. We need to prepare a table of coordinates by copying from the raw data file with columns including `Latitude`, `Longitude` and `site_name`. Ater this, upload this table to GEE (Google Earth Engine).<br>
### Extract_Altitude
>Load the table<br>
```
var csvFile = ee.FeatureCollection("Your table path on GEE");
```
>Run the code<br>
>Download the csv file from your `Google Drive` manually<br>

### Extract_Climate
`Abatzoglou, J.T., S.Z. Dobrowski, S.A. Parks, K.C. Hegewisch, 2018, Terraclimate, a high-resolution global dataset of monthly climate and climatic water balance from 1958-2015, Scientific Data 5:170191, doi:10.1038/sdata.2017.191`
>Load the table<br>
>Extract all possible covariances statistic values, including `Mean`, `Min`, `Max`, `Median` and `Std` of the variables in the table below.<br>

<table class="MsoNormalTable" border="0" cellspacing="0" style="border-collapse:collapse;width:642.1000pt;margin-left:4.6500pt;
mso-table-layout-alt:fixed;border:none;mso-padding-alt:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;"><tbody><tr style="height:13.5000pt;"><td width="53" valign="center" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Name</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="center" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Units</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="48" valign="center" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Min</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="52" valign="center" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Max</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="43" valign="center" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Scale</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="613" valign="center" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Description</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">aet</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">mm</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">3140</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.1</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Actual evapotranspiration, derived using a one-dimensional soil water balance model</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">def</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">mm</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">4548</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.1</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Climate water deficit, derived using a one-dimensional soil water balance model</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">pdsi</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="text-align:left;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p>&nbsp;</o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">-4317</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">3418</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.01</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Palmer Drought Severity Index</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">pet</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">mm</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">4548</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.1</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Reference evapotranspiration (ASCE Penman-Montieth)</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">pr</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">mm</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">7245</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="text-align:left;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p>&nbsp;</o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Precipitation accumulation</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">ro</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">mm</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">12560</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="text-align:left;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p>&nbsp;</o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Runoff, derived using a one-dimensional soil water balance model</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">soil</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">mm</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">8882</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.1</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Soil moisture, derived using a one-dimensional soil water balance model</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">srad</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">W/m</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;"><font face="Times New Roman">²</font></span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">5477</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.1</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Downward surface shortwave radiation</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">swe</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">mm</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">32767</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="text-align:left;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p>&nbsp;</o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Snow water equivalent, derived using a one-dimensional soil water balance model</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">tmmn</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;"><font face="Times New Roman">°</font></span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">C</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">-770</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">387</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.1</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Minimum temperature</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">tmmx</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;"><font face="Times New Roman">°</font></span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">C</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">-670</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">576</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.1</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Maximum temperature</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">vap</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">kPa</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">14749</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.001</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Vapor pressure</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">vpd</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">kPa</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">1113</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.01</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Vapor pressure deficit</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr><tr style="height:14.2500pt;"><td width="53" valign="top" nowrap="" style="width:39.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">vs</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(55,71,79);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:1.0000pt;"><o:p></o:p></span></b></p></td><td width="46" valign="top" nowrap="" style="width:34.9500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">m/s</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="48" valign="top" nowrap="" style="width:36.1000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="52" valign="top" nowrap="" style="width:39.2500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">2923</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="43" valign="top" nowrap="" style="width:32.3000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">0.01</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td><td width="613" valign="top" nowrap="" style="width:459.7500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:top;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;">Wind-speed at 10m</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-style:normal;font-size:10.5000pt;mso-font-kerning:1.0000pt;"><o:p></o:p></span></p></td></tr></tbody></table>


## `Preprocess.py`<br>
### Log-transformation
Before imputation, we have tested that all body size `BS` features, including `Snout-vent Length(SVL)`, `Body Mass(BM)` and `Condition Factor(CF)`, and bioaccumulated heavy metals are right-skewed. Therefore we applied logarithm transformation to normalize the data.<br>

### Imputation
Our original data `raw_data.xlsx` with `n` rows were collected from the previous articles. Given that there are many empty cells in this file, we used our `Phylo-KNN` to impute the null cells.<br>
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
Our research partially focuses on the HM effects on the Anuran body size `BS`. Considering it is well proved that there are some covariance (including `latitude`, `altitude` and `sex`) impacting their body size, we should eliminate the covariance impacts with<br>
<p align="center">$log_a BS_{Revised} = log_a BS_{Original} - (β_1 × latitude + β_2 × altitude + β_3 × sex)$<p>

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
After reviewing literature, we only select seven `HMs` that are obserrved they change the body size, including `Cd`, `Cr`, `Cu`, `Fe`, `Hg`, `Mn`, `Pb` and `Zn`.<br>
Parallelly, we apply `Quadratic Rgression` to fit the latitudinal distribution patterns of the selected HMs and plot their correspongding geological distribution map.<br>

## `Quantile.py`<br>
### Quantile Regression
To study the body size changing effects of HMs, we apply `quantule regression` to see how the body size in different size ranges are sensitive to the HMs.<br>
Considering the phylogenetic effect may disrupt the regression result, we should use a dot product of the original data and phylogenetic covariance matrix instead of the untreated data.
<p align="center">$Y_{Revised} = Y_{Orginal} · P_{ij}$<p>

### PCA
We may find different body size changing sensitivity to the HMs from the previous subplot of `quantile regression`.<br>
We may get different groups with differentbody size and the corresponding HM bioaccumulative levels.<br>
Consequently, we apply PCA to get four clusters.
> Here we qualitatively classify the body size and HM bioaccumulation as `large size` and `small size`, and `high HMs` and `small HMs`.<br>
> The PCA clusters refer the cross match between those two body size and two HM levels, shown as below.<br>

<table class="MsoTableGrid" border="1" cellspacing="0" style="border-collapse:collapse;border:none;mso-border-left-alt:0.5000pt solid windowtext;
mso-border-top-alt:0.5000pt solid windowtext;mso-border-right-alt:0.5000pt solid windowtext;mso-border-bottom-alt:0.5000pt solid windowtext;
mso-border-insideh:0.5000pt solid windowtext;mso-border-insidev:0.5000pt solid windowtext;mso-padding-alt:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;"><tbody><tr><td width="95" valign="top" style="width:71.6000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></b></p></td><td width="189" valign="top" style="width:142.0500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">High HMs</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></b></p></td><td width="183" valign="top" style="width:137.4500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Low HMs</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></b></p></td></tr><tr style="height:14.6500pt;"><td width="95" valign="top" style="width:71.6000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Small Size</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></b></p></td><td width="189" valign="top" style="width:142.0500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Small Size </span><span style="font-family:宋体;color:rgb(32,33,36);font-weight:normal;
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;"><font face="宋体">× </font></span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">High HMs</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></p></td><td width="183" valign="top" style="width:137.4500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Small Size </span><span style="font-family:宋体;color:rgb(32,33,36);font-weight:normal;
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;"><font face="宋体">× </font></span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Low HMs</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></p></td></tr><tr><td width="95" valign="top" style="width:71.6000pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;background:rgb(241,241,241);"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Large Size</span></b><b><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:bold;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></b></p></td><td width="189" valign="top" style="width:142.0500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Large Size </span><span style="font-family:宋体;color:rgb(32,33,36);font-weight:normal;
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;"><font face="宋体">× </font></span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">High HMs</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></p></td><td width="183" valign="top" style="width:137.4500pt;padding:0.0000pt 5.4000pt 0.0000pt 5.4000pt ;border-left:none;
mso-border-left-alt:none;border-right:none;mso-border-right-alt:none;
border-top:none;mso-border-top-alt:none;border-bottom:none;
mso-border-bottom-alt:none;"><p class="MsoNormal" style="mso-pagination:widow-orphan;text-align:left;vertical-align:middle;"><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Large Size </span><span style="font-family:宋体;color:rgb(32,33,36);font-weight:normal;
font-style:normal;font-size:10.5000pt;mso-font-kerning:0.0000pt;"><font face="宋体">× </font></span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;">Low HMs</span><span style="font-family:'Times New Roman';mso-fareast-font-family:宋体;color:rgb(32,33,36);
font-weight:normal;font-style:normal;font-size:10.5000pt;
mso-font-kerning:0.0000pt;"><o:p></o:p></span></p></td></tr></tbody></table>

### Heatmap
To determine whether the clusters are grouped according to the `trait plasticity` or `evolutionary difference`, we should acquire the correlation matrix between the body size and HM levels.<br>

#### Trait plasticity
The trait plasticity was depicted by ITV (Intraspecific Trait Variance) of each sample site.<br>
#### Evolutionary difference
We use $K_{multi}$ to dipict the evolutionary difference of each group.<br>


_____________________________________________________
## Environment
We executed all python code with Python 3.12.<br>
The related python packages and their versions are listed as follow.<br>

>`numpy` version 1.26.4<br>
>`pandas` version 2.2.2<br>
>`biopython` version 1.85<br>
>`tqdm` version 4.66.5<br>
>`matplotlib` version 3.9.2<br>
>`fancyimpute` version 0.7.0<br>
>`sciki-learn` version 1.5.1<br>
>`geopandas` version 1.0.1<br>
>`seaborn` version 0.11.0<br>
>`pillow` version 10.4.0<br>
>`mpl_chord_diagram` version 0.4.1<br>
>`openpyxl` version 3.1.5<br>
>`scipy` version 1.13.1<br>


For more details, please contact via zxxie@iue.ac.cn



