# -*- coding: utf-8 -*-
 

import ABCO as fyp 
import numpy as np


import matplotlib.pyplot as plt

"""This function runs the algorithm 50 times and stores the time taken, best solutions and error rate for each run. It stores this information in a text file """
def graph(pop_s, step_s, step_exp, step_explt, step_tum, dim, lb, ub, threshold, s, k, weight_mode, variable,val1,val2,unchangedthreshold):
    arr=[]
    arrsol=[]
    arrerr=[]
    for i in range(50):
        fyp.BCO(pop_s, step_s, step_exp, step_explt, step_tum, dim, lb, ub, threshold, s, k, weight_mode, variable,val1,val2,unchangedthreshold) #step_s, step_exp, steps_explt, step_tum, dim, lb, ub, threshold, s, k, 
        arr.insert(i,fyp.BCO.final_time)
        arrsol.insert(i, fyp.BCO.bestsolution)
        arrerr.insert(i, fyp.BCO.error)
    #print("Average time for this algorithm on function",sum(arr)/10)
    #print("Average solution for this algorithm:",sum(arrsol)/10)
    #print(arrsol)
    with open(f'testing.txt' , 'a') as f:  
        #print(f'#threshold={threshold}, popsplit={s}, kneighbours={k}',file=f)                                     
        print(f'{variable}times.append({arr})' ,file=f)
        #print(f'{variable}error= np.round({arrerr},decimals=4)',file=f)
        print(f'{variable}.append(np.round({arrerr},decimals=4))',file=f)

functionVariables=['Holders','Holderstimes','Easom','Easomtimes','Goldstein','Goldsteintimes','Rosenbrock','Rosenbrocktimes'
                  ,'Ackley','Ackleytimes','Schaffer','Schaffertimes','Rastrigin','Rastrigintimes','Sphere','Spheretimes','Booth','Boothtimes'
                  ,'Himmelblau','Himmelblautimes']

#for i in functionVariables:
    
#    with open(f'testingB.txt' , 'a') as f: 
#        print(f'{i}=[]',file=f)
#global param
#param={}

testDict= {'popsize':[25,50,75], 'stepsize':[1,2,3], 'stepexplore':[4], 'stepexploit':[1,2,3],
           'tumblesteps':[2,3,4], 'threshold':[0.05,0.1,0.3,0.4], 'popsplit':[0.8],'kN': [2]}

"""This function runs the experiments using the dictionary value above """
def testing():
        #index=0
        #for t in thrs:
        #for exp in testDict['stepexplore']: 
         
        #for t in testDict['threshold']:

        for s in testDict['popsplit']:
            for exp in testDict['stepexplore']: 
                    for kn in testDict['kN']:
                            #print(f'testing when threshold={t}, popsplit={s},kn={k}')
                            
                            
                            #print(f'testing when threshold={t}, popsplit={s},kn={k}')
                            #param.update({index: f'when,stepexplore={exp},kn={kn}'}) #popsplit={s},kn={kn}'})
                            #print(param[index])
                            with open(f'testing.txt' , 'a') as f:  
                                    print(f'',file=f)
                                    print(f'#When stepexplore={exp}and popsize=25  knn  results ',file=f)
                                   
                            graph(25, 1, exp, 1, 1, 2, -10, 10, 0.4, s, kn,'min','Holders',8.05502,9.66459, 80)#
     
                            graph(25, 1, exp, 1, 1, 2, -100, 100, 0.05, s, kn,'min','Easom',fyp.pi,fyp.pi, 80)#
     
                            graph(25, 1, exp, 1, 1, 2, -2, 2, 0.05, s, kn,'min','Goldstein',0,-1, 80)#
     
                            graph(25, 1, exp, 1, 1, 2, -5, 10, 0.4, s, kn,'min','Rosenbrock',1,1, 80)#
     
                            graph(25,1, exp, 1, 1, 2, -5, 5, 0.05, s, kn,'min','Ackley',0,0, 80)#
     
                            graph(25, 1, exp, 1, 1, 2, -100, 100, 0.3, s, kn,'min','Schaffer',0,0, 80)#
     
                            graph(25, 1, exp, 1, 1, 1, -5.12, 5.12, 0.05,s, kn,'min','Rastrigin',0,fyp.pi, 80)#
     
                            graph(25, 1, 1, 1, 1, 1, -100, 100, 0.4, 0.5, 15,'min','Sphere',0,0, 80)# def knn=15 and popsize=25 or knn=2 when pop=15
     
                            graph(25, 1, exp, 1, 1, 2, -10, 10, 0.3, s, kn,'min','Booth',1,3, 80)#
     
                            graph(25, 1, exp, 1, 1, 2, -5, 5, 0.05, 0.3, 5,'min','Himmelblau',3.0,2.0, 80)# popsplit=0.3
                            print("#done")
                            index+=1
                                
        
testing() 


 
