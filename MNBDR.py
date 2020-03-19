# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 16:14:37 2019

@author: ppttssy
"""

import numpy as np
from sklearn import preprocessing
class MNBDR_Map(object):                               
    def  __init__(self,Selectd_module_num,Ess_thr=148):
        self.SELECTED_MODULE_NUM,self.ESS_THR=Selectd_module_num,Ess_thr
        self.ALL_MODULE_NUM=self.Add_Module_num()
        self.Cluster_martix=self.Add_Gonetwork()
        self.Disease_expression_value=np.zeros(self.ALL_MODULE_NUM)
        self.Dicts=self.Add_Dict()
        self.Cancer_ce_D,self.DRUG_EXP_LENGTH=self.Add_Cancell()
        self.Add_Base()
        self.first()
    def  Add_Module_num(self):
        with open('Data/PPI_Cluster_network.txt','r') as fs: 
            for lines in fs:
                 lines=lines.strip('\n').strip('\t').split('\t')   
                 return len(lines)
    def  Add_Gonetwork(self):
        Cluster_martix=np.zeros((self.ALL_MODULE_NUM,self.ALL_MODULE_NUM))    # Import module network  
        with open('Data/PPI_Cluster_network.txt','r') as fs: 
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
              if line[index] not  in Dicts:          
                  Dicts[line[index]]=str(lines_index)
              else:
                  Dicts[line[index]]+=','+str(lines_index)     
            lines_index+=1
        return Dicts
    def Add_Cancell(self):          #Import the cancer cell_line data
        with open('Data/AVE_DMSO_MCF7_ctl_vehicle_ALL h.txt','r') as r:        
            for lines in r:
                lines=lines.strip('\n').split()
                Drug_exp_length=len(lines)
                break
        Cancer_ce_D=np.zeros(Drug_exp_length)
        with open('Data/AVE_DMSO_MCF7_ctl_vehicle_ALL h.txt','r') as r:        
            for lines in r:
                lines=lines.strip('\n').split()
                for index in range(len(lines)):
                    Cancer_ce_D[index]=float(lines[index])
        return Cancer_ce_D,Drug_exp_length
                
    '''
    Add_Base:
        Averaged gene expression data were imported (disease and control groups)
        Calculate the (Imp1,Imp2,...,Impn) on this disease for each denser module      
    '''
    def Add_Base(self):
            Normal,Disease,dict1mid=[],[],{}        
#            Normal=np.zeros(14477)
#            Disease=np.zeros(14477)
#            with open ('Data/BRCA_cancer_expression.txt', 'r') as cancer:
#                 with open ('Data/BRCA_normal_expression.txt', 'r') as normal:
#                     for index,(lines_cancer,line_normal) in enumerate(zip(cancer,normal)):
#                         Disease[index],Normal[index]=float(lines_cancer.strip('\n')),float(line_normal.strip('\n'))
            with open ('Data/Disease.txt','r') as disease:
                     for index,lines in enumerate(disease):
                         Disease.append(float(lines.strip('\n').split('\t')[2])),Normal.append(float(lines.strip('\n').split('\t')[1]))
            Normal,Disease=np.array(Normal),np.array(Disease)
            for nn in range(self.ALL_MODULE_NUM):
              dict1mid[nn]=[]
            with open('Data/Disease.txt', 'r') as gene:  #基因名
                lines_index=0
                for g in gene:             
                  g=g.strip('\n').split('\t')
                  if g[0] in self.Dicts:
                     for s in range(len(self.Dicts[g[0]].split(','))):                                        
                         dict1mid[int(self.Dicts[g[0]].split(',')[s])].append(np.log2(Disease[lines_index]+1)-np.log2(Normal[lines_index]+1))
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
      with open('Data/Drug_gene.txt', 'r') as gene:  #基因名
         line=0
         dict2mid={}
         for nn in range(self.ALL_MODULE_NUM):
             dict2mid[nn]=[]
         for j in range(self.ALL_MODULE_NUM):
             Drug_expression_value[j]=0       
         for g in gene:          
            g=g.strip('\n')
            if g in self.Dicts:
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
          if Uplist[j] in Dis_drugdict:
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
            dist=np.linalg.norm(ago-temp)
        back=back[np.argsort(-temp)]               
        temp=temp[np.argsort(-temp)]
        Import_list=[]       
        with open('Result.txt','w') as result:
              result.write('Disease '+'Drug  '+"Score        \n")
        with open('Result.txt','a') as result:
          All_score=[]
          with open('Data/AVE_MCF7_Min.txt','r') as disease: 
                Drug_in_D=np.zeros(self.DRUG_EXP_LENGTH)    
                for lines in disease:   
                      lines=lines.strip('\n').split()
                      for index in range(1,len(lines)):
                          Drug_in_D[index-1]=float(lines[index])
                      back_dc=np.zeros(self.ALL_MODULE_NUM)
                      Dis_drugdict={}    
                      Dis_drugdict2={} 
                      ac=Drug_in_D-self.Cancer_ce_D
                      ac=ac[np.argsort(ac)]  
                      sums=0
                      for index in range(50):                       #Drug expression with Ess filtration
                           sums+=(abs(ac[index])+abs(ac[self.DRUG_EXP_LENGTH-1-index]))
                      if sums>=self.ESS_THR:
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
                            All_score.append(('BRCA',lines[0],round(S,5)))
#                            result.write('BRCA '+lines[0]+'  '+str(S)+"        \n")
          All_score.sort(key=lambda x:x[2],reverse=True)
          for index in range(len(All_score)):
              result.write(All_score[index][0]+' '+All_score[index][1]+'  '+str(All_score[index][2])+"        \n")

          
 


