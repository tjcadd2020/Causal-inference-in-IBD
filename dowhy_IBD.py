#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 09:50:07 2019

@author: gaosheng
"""
import numpy as np
import pandas as pd

import dowhy
from dowhy.do_why import CausalModel
import dowhy.datasets

#data = dowhy.datasets.linear_dataset(beta=10,
#        num_common_causes=5,
#        num_instruments = 2,
#        num_samples=10000,
#        treatment_is_binary=True)
#df = data["df"]
#print(df.head())
#print(data["dot_graph"])
#print("\n")
#print(data["gml_graph"])




# =============================================================================
# Causal Inference
# =============================================================================
def write_dot(graph_edge): ### convert edge csv into .dot 
    dot_str = 'digraph {'
    for i in range(graph_edge.shape[0]):
        edge = ' ' + graph_edge.iloc[i][0] + '->' + graph_edge.iloc[i][1] + ';'
        dot_str = dot_str + edge
    dot_str = dot_str + '}'
    return dot_str



CD_graph_edge = pd.read_csv('CD_graph_edges.csv')
CD_causal_inference_df = pd.read_csv('CD_causal_inference_df.csv')# row was sample, col was gene, the last col was group/condition info
CD_dot = write_dot(CD_graph_edge)  ## convert graghfile to dot
CD_DEG = list(CD_graph_edge.loc[CD_graph_edge['target']=='IBD']['source'])




# =============================================================================
# Significance Evaluation of Causal Inference
# =============================================================================

CD_ci_df_temp = CD_causal_inference_df.copy()


import random

class permutation_test:
    def __init__(self,causal_inference_df,dot,DEG):
        self.causal_inference_df = causal_inference_df             
        self.dot = dot
        self.index = DEG
        self.outcome = self.causal_inference_df['IBD'].copy()



    
    def Causal_Inference_class(self):   ######## causal inference process
        Gene = []
        cev_list = []
        i = 0
        for item in self.index:
            print(i)
            i = i + 1
            focus_gene = item
            # With graph
            model = CausalModel(
              data = self.causal_inference_df,
              treatment = focus_gene,
              outcome = 'IBD',
              graph = self.dot
            )
            
            identified_estimand = model.identify_effect()
            
            causal_estimate = model.estimate_effect(identified_estimand,
            method_name="backdoor.linear_regression") ##  LogisitcRegression
            Gene.append(focus_gene)
            
            value = causal_estimate.value[0]  ## 因果推理参数值
            print(value)
            
            expr = causal_estimate.realized_estimand_expr   ##### extract 构造表达式，用于置换检验
            
            
            
            ######## Permutation 
            random_cev = self.permutation_func(expr)
            random_cev.append(value) 
            random_cev = pd.Series(random_cev)
            
            cev_list.append(random_cev)

            
        permutation_df = pd.DataFrame(cev_list,index = self.index).T  #### 所有随机结果CEV 构成df,默认按行叠加
            

        return permutation_df
    
    def permutation_func(self,expr):  ### 蒙特卡洛方法随机999次
        
        random_cev = []

        for i in range(999):
            outcome = np.random.permutation(self.outcome) #### 随机最后一列disease状态值
            features = self.causal_inference_df[expr.split('~')[1].split('+')]
            from sklearn import linear_model
            model = linear_model.LogisticRegression()
            model.fit(features, outcome)
            coefficients = model.coef_
            value = coefficients[0,0]
            random_cev.append(value)
        return random_cev
    
    def compute_significance(self,CEV_series): #### 根据每个series计算pvalue
        CEV_list = list(CEV_series)
        bias = [item for item in CEV_list if item >= CEV_list[-1]]
        sig = len(bias)/1000
        return sig
    
    def applyfunc(self,CEV_df):
        sig_df = CEV_df.apply(self.compute_significance,axis = 'rows') ### 根据每列,apply函数计算所有gene的pvalue
        Causal_Inference = pd.DataFrame({'Gene':self.index, 'Causal_Estimated_Value':list(CEV_df.iloc[-1,:]), 'Pvalue':list(sig_df)} )
        
        return Causal_Inference
    
import time   
time1 = time.time()


######### permutation in CD
import time   
time1 = time.time()   
PT_CD = permutation_test(CD_ci_df_temp,CD_dot,CD_DEG)
CD_CEV_df = PT_CD.Causal_Inference_class()
CD_CI_result = PT_CD.applyfunc(CD_CEV_df)
time2 = time.time()
time = (time2 - time1)/3600
print('\n\n\n it has cost %s Hrs \n\n\n' %time)
CD_CI_result.to_csv('CD_CI_result.csv',header = True, index = True)




