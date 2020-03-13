# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 16:14:37 2019

@author: ppttssy
"""

import numpy as np
from sklearn import preprocessing
class GOcmap(object):  
    Cancer_ce_D=np.zeros(12328)                                 
    Disease_expression_value=np.zeros(116)
    def  __init__(self,Selectd_module_num,Ess_thr,All_module_num):
        self.SELECTED_MODULE_NUM,self.ESS_THR,self.ALL_MODULE_NUM=Selectd_module_num,Ess_thr,All_module_num
        self.Cluster_martix=self.Add_Gonetwork()
        self.Dicts=self.Add_Dict()
        self.Cancer_ce_D=self.Add_Cancell()
        self.Add_Base()
        self.first()
    def  Add_Gonetwork(self):
        Cluster_martix=np.zeros((self.ALL_MODULE_NUM,self.ALL_MODULE_NUM))    # Import module network  
        with open('Data/PPI_Clusetr_network.txt','r') as fs: 
            lines_index=0
            for lines in fs:
                 lines=lines.strip('\n').strip('\t').split('\t')               
                 for index in range(len(lines)):
                     Cluster_martix[lines_index][index]=float(lines[index])
                 lines_index+=1
            matrix_normalized = preprocessing.normalize(Cluster_martix, norm='l1') 
            matrix_normalized=np.transpose(matrix_normalized)                     
            return  matrix_normalized 
    def  Add_Dict(self):  # Create a dict ,gene_id as key ,module_id which gene belongs to as Disease_expression_value
      Dicts={}          
      with open('Data/PPI_Clusetr.txt','r') as f:       
        lines_index=0                                    
        for line in f:
            line=line.strip('\n').strip('\t').split('\t') 
            for index in range(0,len(line)):
              if not Dicts.has_key(line[index]):          
                  Dicts[line[index]]=str(lines_index)
              else:
                  Dicts[line[index]]+=','+str(lines_index)     
            lines_index+=1
        return Dicts
    def Add_Cancell(self):          #Import the cancer cell_line data
        Cancer_ce_D=np.random.rand(12328)
        return Cancer_ce_D
    '''
    Add_Base:
        Averaged gene expression data were imported (disease and control groups)
        Calculate the (Imp1,Imp2,...,Impn) on this disease for each denser module      
    '''
    def Add_Base(self):
            dict1mid={}        
            Normal=np.random.rand(12328)
            Disease=np.random.rand(12328)
            for nn in range(self.ALL_MODULE_NUM):
              dict1mid[nn]=[]
            with open('Data/Genenanme12328.txt', 'r') as gene:  #基因名
                lines_index=0
                for g in gene:             
                  g=g.strip('\n')
                  if self.Dicts.has_key(g):
                     for s in range(len(self.Dicts[g].split(','))):                                        
                         dict1mid[int(self.Dicts[g].split(',')[s])].append(np.log2(Disease[lines_index]+1)-np.log2(Normal[lines_index]+1))
                  lines_index+=1  
                for index in range(self.ALL_MODULE_NUM):
                    if dict1mid[index]!=[] :
                        if np.max(np.array(dict1mid[index]))>0 and np.min(np.array(dict1mid[index]))>0:
                            self.Disease_expression_value[index]=np.max(np.array(dict1mid[index]))
                        elif np.max(np.array(dict1mid[index]))<0 and np.min(np.array(dict1mid[index]))<0:
                            self.Disease_expression_value[index]=abs(np.min(np.array(dict1mid[index])))
                        else:                                                                  
                                 self.Disease_expression_value[index]=(abs(np.max(np.array(dict1mid[index])))  +   abs(np.min(np.array(dict1mid[index]))))                                                                                                           
                    else:
                         self.Disease_expression_value[index]=0
    def addgenename(self,Drug_in_D): 
      Drug_expression_value=np.zeros(self.ALL_MODULE_NUM)   
      with open('Data/Genenanme12328.txt', 'r') as gene:  #基因名
         line=0
         dict2mid={}
         for nn in range(self.ALL_MODULE_NUM):
             dict2mid[nn]=[]
         for j in range(self.ALL_MODULE_NUM):
             Drug_expression_value[j]=0       
         for g in gene:          
            g=g.strip('\n')
            if self.Dicts.has_key(g):
               for s in range(len(self.Dicts[g].split(','))):
                        dict2mid[int(self.Dicts[g].split(',')[s])].append((Drug_in_D[line]-self.Cancer_ce_D[line]))
            line+=1         
         for qq in range(self.ALL_MODULE_NUM):
             if dict2mid[qq]!=[] :
                            if np.max(np.array(dict2mid[qq]))>0 and np.min(np.array(dict2mid[qq]))>0:
                                Drug_expression_value[qq]=np.max(np.array(dict2mid[qq]))
                            elif np.max(np.array(dict2mid[qq]))<0 and np.min(np.array(dict2mid[qq]))<0:
                                Drug_expression_value[qq]=abs(np.min(np.array(dict2mid[qq])))
                            else:                            
                                Drug_expression_value[qq]=(abs(np.max(np.array(dict2mid[qq])))  +   abs(np.min(np.array(dict2mid[qq]))))   
             else:
                                  Drug_expression_value[qq]=0   
         return Drug_expression_value
    def Import_count(self,Dis_drugdict,Dis_drugdict2,Uplist): 
        score=0                                          
        for j in range(len(Uplist)):
          if Dis_drugdict.has_key(Uplist[j]):
                    score+=Dis_drugdict[Uplist[j]]/(abs(Dis_drugdict2[Uplist[j]]-(j))+1)  
        return score
    def first(self):  
        temp=np.zeros(self.ALL_MODULE_NUM)
        back=np.zeros(self.ALL_MODULE_NUM)
        for i in range(self.ALL_MODULE_NUM):
          temp[i]=self.Disease_expression_value[i]   
          back[i]=i
        dist=1
        while dist>=0.0000000001:               #Combination disease Imp,modular network with Pagerank algorithm
            ago=temp                              #Then,Sorting the new Imp_value to pick the important 
            temp=0.85*np.dot(self.Cluster_martix,temp)+0.15*self.Disease_expression_value
            dist=numpy.linalg.norm(ago-temp)
        back=back[np.argsort(-temp)]               
        temp=temp[np.argsort(-temp)]
        Import_list=[]       
        with open('Result.txt','w') as result:
            back_dc=np.zeros(self.ALL_MODULE_NUM)
            Dis_drugdict={}    
            Dis_drugdict2={} 
            Drug_in_D=np.zeros(12328)       
            ac=Drug_in_D-self.Cancer_ce_D
            ac=ac[np.argsort(ac)]  
            sums=0
            for index in range(50):                       #Drug expression with Ess filtration
                   sums+=(abs(ac[index])+abs(ac[12327-index]))
            if sums<=self.ESS_THR:
                Drug_expression_value=self.addgenename(Drug_in_D) #Calculate the Imp_value of the drug
                temp_dc=np.zeros(self.ALL_MODULE_NUM)     
                for index in range(self.ALL_MODULE_NUM):
                       temp_dc[index]=1.0/self.ALL_MODULE_NUM
                       back_dc[index]=index
                temp_dc=Drug_expression_value
                back_dc=back_dc[np.argsort(-temp_dc)]             
                temp_dc=temp_dc[np.argsort(-temp_dc)]       
                for s in range(self.ALL_MODULE_NUM):                
                   Dis_drugdict2[back_dc[s]]=s
                   Dis_drugdict[back_dc[s]]=temp_dc[s]         
                del Import_list[:]
                for count in range(self.SELECTED_MODULE_NUM):                        
                  Import_list.append(float(back[count]))
                S=self.Import_count(Dis_drugdict,Dis_drugdict2,Import_list)  #Calculate the S between the drug-disease       
                result.write('disease_name  '+'drug_name  '+str(S)+"        \n")



          
 


