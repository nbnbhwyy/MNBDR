# MNBDR
README file for Python source code supporting the paper "MNBDR: A Module Network Based Method for Drug Repositioning".
## Requirements
* Python = 3.X
* Sklearn 
* Numpy <br>
Numpy and Sklearn can be installed using pip.
## Detail
### Input files
__PPI_Clusetr__  --dense clusters in the SRTING PPI network, geneId = "Entrez". <br>
Example: 
```
5693 56929 10199 55313	 
522 6390 9551 6391 
51542 9827 25782 636 
29883 485 23019 10605
```
__PPI_Clusetr_network__  --the edge of dense clusters. <br>
Example: 
```
0	1	1	1
1	0	1	1
0	1	0	0
1	0	1	1
```
__Drug_gene__  --Since Disease Data and drug data often come from different databases, gene's order and number <br>
don't usually match. Therefore, we provide a list of drug genes to facilitate the extraction of drug expression information. <br>
Example: 
```
5693 
56929 
10199 
55313	 
29883 
```
__Disease__  --tab delimited, Expression data of disease, the column are Gene_Entrez, Normal group, Disease group,respective (In the test cases, we provide data for breast cancer patients). 
Example: 
```
1	0.10158704	0.161963237
29974	0.01379371	0.013093553
2	242.8467729	115.1052494
144568	0.197790261	2.093597589
53947	6.994426346	6.398611842
```
__AVE_DMSO_MCF7_ctl_vehicle_ALL h__  --Expression data of Cell Lines after drug stimulation (Control group) (In the test cases, we provide data for breast cancer Cell Lines). Example: 
```
0.07991	0.06652	-0.06172	0.04316	
```
__AVE_MCF7_Min__  --Expression data of Cell Lines after drug stimulation (treatment group) (In the test cases, we provide data for breast cancer Cell Lines). Example: 
```
lamotrigine -0.31464 0.2192 -0.0083 -0.05618
glimepiride 0.8858 -0.26607 -0.13603 -0.16573
```
__RUN__ <br> 
In this test case, we used breast cancer as an example and provided a small data set of drug responses. <br> 
If you run "test.py" , you can get a scores for breast cancer and drugs. <br> 
Example: 
```
Disease Drug  Score        
BRCA pentobarbital  46.37012        
BRCA SN-38  45.97034        
BRCA epirubicin  44.385        
BRCA daunorubicin  40.46897         
```
__Custom dataset__
You can also apply our method to your Dataset by replacing the data in DATA. <br> 
For example, you can provide disease data in our format, then you can get scores for new diseases and drugs.<br> 
Similarly, you can also replace drug data or module network data.<br> 
<br> If you have any questions,contact the maintainer of this package: HG_Chen 13247702278@163.com
