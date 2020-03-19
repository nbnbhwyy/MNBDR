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
__Disease_gene and Drug_gene__  --Since Disease Data and drug data often come from different databases, gene's order and number 
don't usually match, so we provide two different gene lists for the code to extract the information. <br>
Example: 
```
5693 
56929 
10199 
55313	 
29883 
```
__BRCA_normal_expression and BRCA_cancer_expression__  --Breast Cancer Disease Samples and normal samples, respectively. 
The preprocessing process is as described in the article. <br>
Example: 
```
0.10158704
0.01379371
242.8467729
0.197790261
6.994426346
0.080435
8.393741522
```
__AVE_DMSO_MCF7_ctl_vehicle_ALL h__  --Expression data of Breast Cancer Cell Lines after drug stimulation (Control group). <br>
Example: 
```
0.07991	0.06652	-0.06172	0.04316	
```
__AVE_MCF7_Min__  --Expression data of Breast Cancer Cell Lines after drug stimulation (treatment group). <br>
Example: 
```
lamotrigine	-0.31464	0.2192	-0.0083	-0.05618
glimepiride	0.8858	-0.26607	-0.13603	-0.16573
```
__RUN__ <br> 
In this test case, we used breast cancer as an example and provided a small data set of drug responses. <br> 
If you run "test.py" , you can get a scores for breast cancer and drugs. <br> 
Example: 
```
BRCA lamotrigine  14.60946508783487        
BRCA glimepiride  7.334471041208792        
BRCA pentobarbital  46.37012152727272        
BRCA darifenacin  4.7773773982351715        
BRCA amiloride  9.592314238066283        
BRCA chlorpromazine  11.531387122847168        
```
You can also apply our method to your Dataset by replacing the data in DATA. <br> 
<br> If you have any questions,contact the maintainer of this package: HG_Chen 13247702278@163.com
