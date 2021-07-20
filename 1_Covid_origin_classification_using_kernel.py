#!/usr/bin/env python
# coding: utf-8

# In[54]:


import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.decomposition import TruncatedSVD
import random
# import seaborn as sns
import os.path as path
import os
import matplotlib
import matplotlib.font_manager
import matplotlib.pyplot as plt # graphs plotting
from Bio import SeqIO # some BioPython that will come in handy
#matplotlib inline
import numpy
import csv 

from matplotlib import rc

from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso, LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.preprocessing import StandardScaler

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics
from sklearn import svm

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics

import pandas as pd
import sklearn
from sklearn import preprocessing
from sklearn.model_selection import train_test_split 
from sklearn.preprocessing import StandardScaler  
from sklearn.neural_network import MLPClassifier 
from sklearn.metrics import classification_report, confusion_matrix 

from sklearn.neighbors import KNeighborsClassifier

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix

from pandas import DataFrame

from sklearn.model_selection import KFold 
from sklearn.model_selection import RepeatedKFold

from sklearn.metrics import confusion_matrix

from numpy import mean


from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import roc_auc_score

import statistics

from wordcloud import WordCloud, STOPWORDS
from itertools import chain
from random import randint
import random

import statistics

from sklearn.cluster import KMeans

from sklearn.datasets import load_digits
from sklearn.decomposition import KernelPCA

import math
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import confusion_matrix

from sklearn.metrics import f1_score


# for Arial typefont
matplotlib.rcParams['font.family'] = 'Arial'


## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
# matplotlib.rcParams['mathtext.fontset'] = 'cm'

## for LaTeX typefont
# matplotlib.rcParams['mathtext.fontset'] = 'stix'
# matplotlib.rcParams['font.family'] = 'STIXGeneral'

## for another LaTeX typefont
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

# rc('text', usetex = True)

print("Packages imported")


# In[55]:


# # path = "D:/University/RA/covid_origin/Dataset/variant/"
# path = "D:/University/RA/JCB_Viral_Host_Classification/Dataset/"

# # file_name = "spikeprot0702"
# # input_file_name = path + file_name + ".fasta"

# file_name = "protein_sequences_alligned"
# input_file_name = path + file_name + ".txt"

# aa = []
# with open(input_file_name) as infile:
#     for line in infile:
#         aa.append(line)
        
# allinged_prot_sequences = aa


# In[56]:


# import mne
path = "D:/University/RA/covid_origin/Dataset/variant/"
# path = "D:/University/RA/JCB_Viral_Host_Classification/Dataset/"

file_name = "spikeprot0702"
input_file_name = path + file_name + ".fasta"

# file_name = "my_output"
# input_file_name = path + file_name + ".meg"


sequences_dictionary = {sequence.id : sequence.seq for sequence in SeqIO.parse(input_file_name,'fasta')}

# sequences_dictionary = mne.io.read_raw_fif(input_file_name)

print("Data Loaded")


# In[57]:


(sequences_dictionary)
# key value pairs
for i in sequences_dictionary:
    aa = sequences_dictionary[i]
    aaa = aa.split('|')
    print(i,aaa[0])


# In[58]:



path = "D:/University/RA/covid_origin/Dataset/variant/"

file_name = "variant_surveillance"
input_file_name = path + file_name + ".tsv"

data_variant = pd.read_csv (input_file_name, sep = '\t')


# In[59]:


deflines = [entry for entry in sequences_dictionary.keys()]             # create a list of deflines
protein_sequences = [entry for entry in sequences_dictionary.values()]  # create a list of protein sequences 


# In[60]:


data_variant_1 = (data_variant.loc[:,["Accession ID","Pango lineage","Type","Host"]])
data_variant_1


# # Check for Missing values and remove the

# In[61]:


# missing_val_check = []
# for i in range(0, len(deflines)):
#     length_check = len(deflines[i])
#     if length_check<50:
#         print(i,', Length = ', len(deflines[i]))
#         missing_val_check.append(i)


# In[62]:


# #removing missing values rows
# for i in range(0, len(missing_val_check)):
#     deflines.remove(deflines[missing_val_check[i]-i])
#     protein_sequences.remove(protein_sequences[missing_val_check[i]-i])
# #     targets.remove(targets[missing_val_check[i]-i])


# In[63]:


variant_id = []
variant_name = []
variant_host_name = []
# variant_name_full = []
for tt in range(len(data_variant_1)):
    variant_tmp = data_variant_1.iloc[tt]
    variant_id.append(variant_tmp[0])
    variant_name.append(variant_tmp[1])
    variant_host_name.append(variant_tmp[3])


# In[70]:


idx = pd.Index(variant_name) # creates an index which allows counting the entries easily
print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
aq = idx.value_counts()
print(aq[0:60])


# In[65]:


tmp_host_name = []
for i in range(0, len(deflines)):
#     print("i = ",i,"/",len(deflines))
    aqw = deflines[i].split('|')
    if(len(deflines[i].split("|"))>=7):
        tmp_host_name.append(aqw[6])
#     else:
#         print(deflines[i])
        if(aqw[6]=='Environment'):
            print(deflines[i])


# In[66]:


idx = pd.Index(tmp_host_name) # creates an index which allows counting the entries easily
print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
aq = idx.value_counts()
print(aq)


# In[68]:


len(tmp_host_name)


# In[ ]:


# deflines_red = []
# for i in range(0, len(deflines)):
# #     print("i = ",i,"/",len(deflines))
#     aqw = deflines[i].split('|')
#     if(len(deflines[i].split("|"))>=7):
# #         if(len(protein_sequences[i])==1274):
#         if(aqw[6]=='Environment'):
#             deflines_red.append(deflines[i])


# In[ ]:


# deflines_red2 = []
# for i in range(0, len(deflines)):
# #     print("i = ",i,"/",len(deflines))
#     aqw = deflines[i].split('|')
#     if(len(deflines[i].split("|"))>=7):
# #         if(len(protein_sequences[i])==1274):
#         if(aqw[6]=='Panthera'):
#             deflines_red2.append(deflines[i])
            
# deflines_red2


# In[ ]:


# total_rand_nums = 10
# random_int_vals = random.sample(range(0, 30000), total_rand_nums)
# random_int_vals.sort()
# random_int_vals

# B.1.1.7       967672
# B.1.617.2     113835
# B.1.2          95191
# B.1            83781
# B.1.177        74273
# B.1.1          47835
# B.1.526        46784
# P.1            46761
# B.1.429        34088
# B.1.351        26831


# B.1.160        25476
# B.1.1.519      20914
# B.1.1.214      17109
# B.1.427        14700
# B.1.258        13914
# B.1.221        13351
# B.1.177.21     13048
# D.2            12745
# B.1.243        11597
# B.1.596        10170


# In[160]:


cities = []
protein_sequences_old = []
deflines_old = []
variant_lst = []

city_names_only = []
# city_name_orig_save_2 = []

total_rand_nums = 20000
random_int_vals = random.sample(range(0, 700000), total_rand_nums)
random_int_vals.sort()

# first variant
counter = 1
counter_2 = 0

uk_counter = 0
sec_counter = 0

#       967672
#      113835
#           95191
#             83781
#         74273
#           47835
#         46784
#             46761
#         34088
#         26831

for i in range(0, len(deflines)):
    aqw = deflines[i].split('|')
    if(len(deflines[i].split("|"))>=7):
        if(len(protein_sequences[i])==1274):
            if(counter<total_rand_nums and random_int_vals[counter]==counter_2): #get random human sequences
                tmp_cariant = (deflines[i].split("|")[3])
                find_variant_index = [ir for ir, nb in enumerate(variant_id) if nb == tmp_cariant] # List comprehension
                if(len(find_variant_index)>0): 
#                     if(variant_name[find_variant_index[0]]=='B.1.160'):
                    if(variant_name[find_variant_index[0]]=='B.1.1.7'):
                        if(uk_counter<1500):
                            print("counter = ",counter, ", First Variant :",variant_name[find_variant_index[0]])
                            variant_lst.append(variant_name[find_variant_index[0]])
                            aaa = (deflines[i].split("|")[1])
                            aaaa = aaa.split("/")[1]
                            cities.append(aaaa)
                            protein_s1 = protein_sequences[i]
                            protein_sequences_old.append(protein_s1)  # whole spike protein
                            deflines_old.append(deflines[i])
                            temp = aqw[5]
                            temp2 = temp.split('^')
                            temp3 = temp2[2]
                            city_names_only.append(temp3)
                            counter = counter+1

                            uk_counter = uk_counter+1
#                     elif(variant_name[find_variant_index[0]]=='B.1.1.519'):
                    elif(variant_name[find_variant_index[0]]=='B.1.617.2'):
                        if(sec_counter<1500):
                            print("counter = ",counter, ", second Variant :",variant_name[find_variant_index[0]])
                            variant_lst.append(variant_name[find_variant_index[0]])
                            aaa = (deflines[i].split("|")[1])
                            aaaa = aaa.split("/")[1]
                            cities.append(aaaa)
                            protein_s1 = protein_sequences[i]
                            protein_sequences_old.append(protein_s1)  # whole spike protein
                            deflines_old.append(deflines[i])
                            temp = aqw[5]
                            temp2 = temp.split('^')
                            temp3 = temp2[2]
                            city_names_only.append(temp3)
                            counter = counter+1

                            sec_counter = sec_counter+1
#                     elif(variant_name[find_variant_index[0]]=='B.1.1.214'):
                    elif(variant_name[find_variant_index[0]]=='B.1.2'):
                        print("counter = ",counter, ", third Variant :",variant_name[find_variant_index[0]])
                        variant_lst.append(variant_name[find_variant_index[0]])
                        aaa = (deflines[i].split("|")[1])
                        aaaa = aaa.split("/")[1]
                        cities.append(aaaa)
                        protein_s1 = protein_sequences[i]
                        protein_sequences_old.append(protein_s1)  # whole spike protein
                        deflines_old.append(deflines[i])
                        temp = aqw[5]
                        temp2 = temp.split('^')
                        temp3 = temp2[2]
                        city_names_only.append(temp3)
                        counter = counter+1
#                     elif(variant_name[find_variant_index[0]]=='B.1.427'):
                    elif(variant_name[find_variant_index[0]]=='B.1'):
                        print("counter = ",counter, ", fourth Variant :",variant_name[find_variant_index[0]])
                        variant_lst.append(variant_name[find_variant_index[0]])
                        aaa = (deflines[i].split("|")[1])
                        aaaa = aaa.split("/")[1]
                        cities.append(aaaa)
                        protein_s1 = protein_sequences[i]
                        protein_sequences_old.append(protein_s1)  # whole spike protein
                        deflines_old.append(deflines[i])
                        temp = aqw[5]
                        temp2 = temp.split('^')
                        temp3 = temp2[2]
                        city_names_only.append(temp3)
                        counter = counter+1
#                     elif(variant_name[find_variant_index[0]]=='B.1.258'):
                    elif(variant_name[find_variant_index[0]]=='B.1.177'):
                        print("counter = ",counter, ", fifth Variant :",variant_name[find_variant_index[0]])
                        variant_lst.append(variant_name[find_variant_index[0]])
                        aaa = (deflines[i].split("|")[1])
                        aaaa = aaa.split("/")[1]
                        cities.append(aaaa)
                        protein_s1 = protein_sequences[i]
                        protein_sequences_old.append(protein_s1)  # whole spike protein
                        deflines_old.append(deflines[i])
                        temp = aqw[5]
                        temp2 = temp.split('^')
                        temp3 = temp2[2]
                        city_names_only.append(temp3)
                        counter = counter+1
#                     elif(variant_name[find_variant_index[0]]=='B.1.221'):
                    elif(variant_name[find_variant_index[0]]=='B.1.1'):
                        print("counter = ",counter, ", sixth Variant :",variant_name[find_variant_index[0]])
                        variant_lst.append(variant_name[find_variant_index[0]])
                        aaa = (deflines[i].split("|")[1])
                        aaaa = aaa.split("/")[1]
                        cities.append(aaaa)
                        protein_s1 = protein_sequences[i]
                        protein_sequences_old.append(protein_s1)  # whole spike protein
                        deflines_old.append(deflines[i])
                        temp = aqw[5]
                        temp2 = temp.split('^')
                        temp3 = temp2[2]
                        city_names_only.append(temp3)
                        counter = counter+1
#                     elif(variant_name[find_variant_index[0]]=='B.1.177.21'):
                    elif(variant_name[find_variant_index[0]]=='B.1.526'):
                        print("counter = ",counter, ", seventh Variant :",variant_name[find_variant_index[0]])
                        variant_lst.append(variant_name[find_variant_index[0]])
                        aaa = (deflines[i].split("|")[1])
                        aaaa = aaa.split("/")[1]
                        cities.append(aaaa)
                        protein_s1 = protein_sequences[i]
                        protein_sequences_old.append(protein_s1)  # whole spike protein
                        deflines_old.append(deflines[i])
                        temp = aqw[5]
                        temp2 = temp.split('^')
                        temp3 = temp2[2]
                        city_names_only.append(temp3)
                        counter = counter+1
#                     elif(variant_name[find_variant_index[0]]=='D.2'):
                    elif(variant_name[find_variant_index[0]]=='P.1'):
                        print("counter = ",counter, ", eight Variant :",variant_name[find_variant_index[0]])
                        variant_lst.append(variant_name[find_variant_index[0]])
                        aaa = (deflines[i].split("|")[1])
                        aaaa = aaa.split("/")[1]
                        cities.append(aaaa)
                        protein_s1 = protein_sequences[i]
                        protein_sequences_old.append(protein_s1)  # whole spike protein
                        deflines_old.append(deflines[i])
                        temp = aqw[5]
                        temp2 = temp.split('^')
                        temp3 = temp2[2]
                        city_names_only.append(temp3)
                        counter = counter+1
#                     elif(variant_name[find_variant_index[0]]=='B.1.243'):
                    elif(variant_name[find_variant_index[0]]=='B.1.429'):
                        print("counter = ",counter, ", ninth Variant :",variant_name[find_variant_index[0]])
                        variant_lst.append(variant_name[find_variant_index[0]])
                        aaa = (deflines[i].split("|")[1])
                        aaaa = aaa.split("/")[1]
                        cities.append(aaaa)
                        protein_s1 = protein_sequences[i]
                        protein_sequences_old.append(protein_s1)  # whole spike protein
                        deflines_old.append(deflines[i])
                        temp = aqw[5]
                        temp2 = temp.split('^')
                        temp3 = temp2[2]
                        city_names_only.append(temp3)
                        counter = counter+1
#                     elif(variant_name[find_variant_index[0]]=='B.1.596'):
                    elif(variant_name[find_variant_index[0]]=='B.1.351'):
                        print("counter = ",counter, ", tenth Variant :",variant_name[find_variant_index[0]])
                        variant_lst.append(variant_name[find_variant_index[0]])
                        aaa = (deflines[i].split("|")[1])
                        aaaa = aaa.split("/")[1]
                        cities.append(aaaa)
                        protein_s1 = protein_sequences[i]
                        protein_sequences_old.append(protein_s1)  # whole spike protein
                        deflines_old.append(deflines[i])
                        temp = aqw[5]
                        temp2 = temp.split('^')
                        temp3 = temp2[2]
                        city_names_only.append(temp3)
                        counter = counter+1
            else:
                counter_2 = counter_2+1
                        
print(cities)


# In[162]:


# np.unique(variant_lst)

idx = pd.Index(variant_lst) # creates an index which allows counting the entries easily
print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
aq = idx.value_counts()
print(aq)


# In[ ]:


# cities = []
# protein_sequences_old = []
# deflines_old = []
# variant_lst = []

# city_names_only = []
# # city_name_orig_save_2 = []

# total_rand_nums = 1500
# random_int_vals = random.sample(range(0, 30000), total_rand_nums)
# random_int_vals.sort()

# counter = 1
# counter_2 = 0
# for i in range(0, len(deflines)):
# #     print("i = ",i,"/",len(deflines))
#     aqw = deflines[i].split('|')
#     if(len(deflines[i].split("|"))>=7):
# #         if(len(protein_sequences[i])==1274):
#         if(aqw[6]=='Human'):
#             if(counter<total_rand_nums and random_int_vals[counter]==counter_2): #get random human sequences
#                 tmp_cariant = (deflines[i].split("|")[3])
#                 find_variant_index = [ir for ir, nb in enumerate(variant_id) if nb == tmp_cariant] # List comprehension
#                 if(len(find_variant_index)>0): 
#                     if(variant_name[find_variant_index[0]]!='nan'):
#                         if(variant_name[find_variant_index[0]]!='None'):
#                             print("counter = ",counter)
#                             counter = counter+1

#                             variant_lst.append(variant_name[find_variant_index[0]])


#                             aaa = (deflines[i].split("|")[1])
#                             aaaa = aaa.split("/")[1]

#                             cities.append(aaaa)
#                             protein_s1 = protein_sequences[i]
#                             protein_sequences_old.append(protein_s1[14:686])  # S1 spike protein
# #                             protein_sequences_old.append(protein_s1)  # whole spike protein
#                             deflines_old.append(deflines[i])

#                             temp = aqw[5]
#                             temp2 = temp.split('^')
#                             temp3 = temp2[2]

#                             city_names_only.append(temp3)
#             else:
#                 counter_2 = counter_2+1
                        
#         else:
# #             if(len(protein_sequences[i])==1274):
#             tmp_cariant = (deflines[i].split("|")[3])
#             find_variant_index = [ir for ir, nb in enumerate(variant_id) if nb == tmp_cariant] # List comprehension
#             if(len(find_variant_index)>0): 
#                 if(variant_name[find_variant_index[0]]!='nan'):
#                     if(variant_name[find_variant_index[0]]!='None'):
# #                         print("counter = ",counter)
# #                         counter = counter+1

#                         variant_lst.append(variant_name[find_variant_index[0]])


#                         aaa = (deflines[i].split("|")[1])
#                         aaaa = aaa.split("/")[1]

#                         cities.append(aaaa)
#                         protein_s1 = protein_sequences[i]
#                         protein_sequences_old.append(protein_s1)  # S1 spike protein
#                         deflines_old.append(deflines[i])

#                         temp = aqw[5]
#                         temp2 = temp.split('^')
#                         temp3 = temp2[2]

#                         city_names_only.append(temp3)
            
            


# # protein_sequences = protein_sequences_old
# print(cities)


# In[79]:


(deflines_old)

tmp_host_name_2 = []
for i in range(0, len(deflines_old)):
    aqw = deflines_old[i].split('|')
    if(len(deflines_old[i].split("|"))>=7):
        tmp_host_name_2.append(aqw[6])
#         if(aqw[6]=='Human'):


# In[99]:


idx = pd.Index(tmp_host_name_2[0:10000]) # creates an index which allows counting the entries easily
print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
aq = idx.value_counts()
print(aq)


# In[ ]:


final_host_names = []

final_variant_names = []
final_cities = []
final_deflines = []
final_proteins = []

for w in range(len(protein_sequences_old)):
    if(len(protein_sequences_old[w])>0):
        if(tmp_host_name_2[w]=="Human" or tmp_host_name_2[w]=="Neovison" or tmp_host_name_2[w]=="Environment" or tmp_host_name_2[w]=="Felis" or tmp_host_name_2[w]=="Panthera" or tmp_host_name_2[w]=="Canis"):
            final_cities.append(cities[w])
            final_deflines.append(deflines_old[w])
            final_proteins.append(protein_sequences_old[w])
            final_host_names.append(tmp_host_name_2[w])


# In[ ]:


# ofile = open("D:/University/RA/JCB_Viral_Host_Classification/Dataset/my_fasta.fasta", "w")

# from Bio import SeqIO

# az = (zip(final_proteins & final_deflines))
# SeqIO.write(az), "D:/University/RA/JCB_Viral_Host_Classification/Dataset/example.fasta", "fasta")


# with open("D:/University/RA/JCB_Viral_Host_Classification/Dataset/example.fasta", "w") as output_handle:
#     SeqIO.write(az, output_handle, "fasta")


# In[81]:


idx = pd.Index(final_host_names) # creates an index which allows counting the entries easily
print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
aq = (idx.value_counts())
print(aq)


# In[163]:


protein_sequences_final = []
for r in range(len(protein_sequences_old)):
    protein_sequences_final.append(list(protein_sequences_old[r]))


# In[164]:


len(protein_sequences_final)


# In[87]:


len(protein_sequences_final[0])

new_protein_seq = []

for i in range(len(protein_sequences_final)):
    tmp_lst = []
    for j in range(1277):
        tmp_lst.append("X")
    new_protein_seq.append(tmp_lst)


# In[88]:


new_protein_seq_2 = []

for i in range(len(protein_sequences_final)):
    tt = protein_sequences_final[i]
    tt_2 = new_protein_seq[i]
    tt_2[0:len(tt)] = tt
    ls = [x if (x!="*") else "X" for x in tt_2]
    new_protein_seq_2.append(ls)
#     new_protein_seq_2[0:len(tt)] = tt


# In[137]:


len(protein_sequences_final)


# In[96]:


# len(new_protein_seq_2)

flatten_list = list(chain.from_iterable(new_protein_seq_2))

len(np.unique(flatten_list))


# In[165]:


# dataset_save_path = "D:/University/RA/JCB_Viral_Host_Classification/Dataset/"
dataset_save_path = "D:/University/RA/covid_origin/Dataset/"
temp_name_tmp = "First_dataset_Alligned_Host_protein_data_with_variants"
Trial = 1

with open(dataset_save_path + temp_name_tmp + ".csv", 'w', newline='') as file:
    writer = csv.writer(file)
    for i in range(0,len(protein_sequences_old)):
        ccv = protein_sequences_old[i]
#         ccv.replace('\"','')
        writer.writerow(ccv)


# In[145]:


(deflines_old[0]).split("|")[2]


# In[168]:


variant_true_names = []
countries_true_names = []
city_true_name = []
date_true_val = []

for i in range(0,len(new_protein_seq_2)):
    variant_true_names.append(variant_lst[i])
    countries_true_names.append(cities[i])
    city_true_name.append(city_names_only[i])
    date_true_val.append((deflines_old[i]).split("|")[2])
    
    


# In[169]:


# Convert the dictionary into DataFrame
df = pd.DataFrame(variant_true_names)
  
# Using 'Address' as the column name and equating it to the list
df2 = df.assign(country = countries_true_names)

df3 = df2.assign(city = city_true_name)
df4 = df3.assign(date = date_true_val)

df4.columns = ['Variant', 'Country', 'City', 'Date']


# In[170]:


np.unique(variant_true_names)


# In[171]:


idx = pd.Index(city_true_name) # creates an index which allows counting the entries easily
print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
aq = idx.value_counts()
print(aq)


# In[172]:


# dataset_save_path = "D:/University/RA/JCB_Viral_Host_Classification/Dataset/"
dataset_save_path = "D:/University/RA/covid_origin/Dataset/"

temp_name_tmp = "First_dataset_dataframe_Host_names_data"
Trial = 1

with open(dataset_save_path + temp_name_tmp + ".csv", 'w', newline='') as file:
    writer = csv.writer(file)
    for i in range(0,len(df4)):
        writer.writerow(list(df4.iloc[i]))


# # Writing Complete Reduced Data

# In[74]:


final_variant_names = []
final_cities = []
final_deflines = []
final_proteins = []

for w in range(len(protein_sequences_old)):
    if(len(protein_sequences_final[w])>0):
        if(variant_lst[w]=="B.1.1.7" or variant_lst[w]=="B.1.2" or variant_lst[w]=="P.1" or 
           variant_lst[w]=="B.1.526" or variant_lst[w]=="B.1.617.2" or variant_lst[w]=="B.1.429" or 
           variant_lst[w]=="B.1.351" or variant_lst[w]=="B.1" or variant_lst[w]=="B.1.1.519" or variant_lst[w]=="B.1.427"):
            final_variant_names.append(variant_lst[w])
            final_cities.append(cities[w])
            final_deflines.append(deflines_old[w])
            final_proteins.append(protein_sequences_old[w])

# Here are all of the viral species in the dataset: 
#  8536 entries in total
# B.1.1.7      5626
# B.1.2         261
# P.1           259
# B.1.526       259
# B.1.617.2     186
# B.1.429       170
# B.1.351       133
# B.1           128
# B.1.1.519     112
# B.1.427        76


# In[78]:


protein_sequences_final = []
for r in range(len(final_proteins)):
    protein_sequences_final.append(list(final_proteins[r]))


# In[79]:


dataset_save_path = "D:/University/RA/covid_origin/Dataset/variant/preprocessed/"
temp_name_tmp = "New_reduced_data"
Trial = 1

with open(dataset_save_path + temp_name_tmp + "_covid_protein_sequences_trial_" + str(Trial) + ".csv", 'w', newline='') as file:
    writer = csv.writer(file)
    for i in range(0,len(protein_sequences_final)):
        if(len(protein_sequences_final[i])>0):
            writer.writerow(protein_sequences_final[i])


# In[233]:


len(protein_sequences_final[0])


# In[84]:


idx = pd.Index(final_variant_names) # creates an index which allows counting the entries easily
print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
aq = (idx.value_counts())
print(aq)


# # Now run the Java code for Kernel Generation, then run the below code!!!!!!

# # Reading kernel matrix

# In[226]:


# read_path = "D:/University/RA/covid_origin/Code/Kernel_approximation/new_Kernel_k4_m0_Alphabet_Size21_trial_1.txt"
read_path = "D:/University/RA/covid_origin/Code/Kernel_approximation/new_Kernel_k6_m3_Alphabet_Size21_trial_1.txt"

seq_file_small = []
with open(read_path) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        seq_file_small.append(', '.join(row))


# In[227]:


final_jernel_mat = []
mat_2 = []
for y in range(0,len(seq_file_small)):
    aa = seq_file_small[y].split(" ")
    tmp_list = []
    for yy in range(0,len(aa)-1):
        tmp_list.append(int(aa[yy]))
    
    final_jernel_mat.append(tmp_list)


# In[244]:


len(final_jernel_mat[0])


# # Agglomerative Clustering

# In[229]:


################ Kernel Trick ######################
mat_2 = np.zeros((len(final_jernel_mat),len(final_jernel_mat)))

for i in range(0,len(final_jernel_mat)):
    tmp_2 = final_jernel_mat[i][i]
    for j in range(0,len(final_jernel_mat)):
        denom = math.sqrt(tmp_2 * final_jernel_mat[j][j])
        temp = final_jernel_mat[i][j] / (denom)
#         mat_2[i][j] = 1 - temp
        mat_2[i][j] = final_jernel_mat[i][i] + final_jernel_mat[j][j] - (2*(final_jernel_mat[i][j])) # distance kernel trick
################ Kernel Trick ######################

################ Agglomerative Clustering ######################
model = AgglomerativeClustering(affinity='precomputed', n_clusters=6, linkage='average').fit(mat_2)
# print(model.labels_)
clusters_ids = model.labels_
################ Agglomerative Clustering ######################


# # Kernel PCA

# In[225]:


#############  Kernel PCA ################
transformer = KernelPCA(n_components=20, kernel='precomputed')
X_transformed = transformer.fit_transform(final_jernel_mat)
X_transformed.shape
#############  Kernel PCA ################

#############  K-Means ################
kmeans = KMeans(n_clusters=5, random_state=0).fit(X_transformed)
clusters_ids = kmeans.labels_
#############  K-Means ################


# In[230]:


#############  Confusion Matrix ################

unique_labels = np.unique(final_variant_names)
labels_integers = []
for i in range(len(final_variant_names)):
    aa = (np.where(unique_labels==final_variant_names[i]))
    labels_integers.append(aa[0][0])

# confusion_matrix(clusters_ids,labels_integers)

confuse = confusion_matrix(labels_integers,clusters_ids)

confuse_matrix = []
for i in range(0,len(confuse)):
    confuse_matrix.append(list(confuse[i]))
#############  Confusion Matrix ################
pd.DataFrame(confuse_matrix)


# In[224]:


from sklearn import metrics

def purity_score(y_true, y_pred):
    # compute contingency matrix (also called confusion matrix)
    contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
    # return purity
    return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix) 


# purity_score(labels_integers,clusters_ids)


# In[231]:


for qq in range(len(np.unique(clusters_ids))):
    names_holder = []
    for ee in range(len(clusters_ids)):
        if(clusters_ids[ee]==qq):
            names_holder.append(final_variant_names[ee])

    mode_val = statistics.mode(names_holder)
    print("Cluster Name = ",mode_val, ", Cluster ID = ",qq)

#     print("Cluster = ",clusters_ids)
    idx = pd.Index(names_holder) # creates an index which allows counting the entries easily
    print(len(idx),"entries in total")
    aq = (idx.value_counts())
    print(aq,"\n")


# In[166]:


idx = pd.Index(final_variant_names) # creates an index which allows counting the entries easily
print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
aq = (idx.value_counts())
print(aq)


# In[219]:


cluster_number_tmp = 2
cluster_name_tmp = "B.1.1.7"

names_holder = []
for ee in range(len(clusters_ids)):
    if(clusters_ids[ee]==cluster_number_tmp):
        names_holder.append(final_variant_names[ee])

mode_val = statistics.mode(names_holder)
# print("Cluster Name = ",mode_val, ", Cluster ID = ",qq)

#     print("Cluster = ",clusters_ids)
idx = pd.Index(names_holder) # creates an index which allows counting the entries easily
print(len(idx),"entries in total")
aq = (idx.value_counts())
print(aq,"\n")

#############  Confusion Matrix ################

unique_labels = np.unique(names_holder)
labels_integers = []
for i in range(len(names_holder)):
    aa = (np.where(unique_labels==names_holder[i]))
    labels_integers.append(aa[0][0])

pred_labels = []
for ch in range(len(names_holder)):
    pred_labels.append(cluster_name_tmp)

confuse = confusion_matrix(names_holder,pred_labels)

confuse_matrix = []
for i in range(0,len(confuse)):
    confuse_matrix.append(list(confuse[i]))
#############  Confusion Matrix ################

print("F1 macro",f1_score(names_holder, pred_labels, average='macro'))
print("F1 micro",f1_score(names_holder, pred_labels, average='micro'))
print("F1 weighted",f1_score(names_holder, pred_labels, average='weighted'))


# In[213]:


cluster_number_tmp = 2
cluster_name_tmp = "P.1"

names_holder = []
for ee in range(len(clusters_ids)):
    if(clusters_ids[ee]==cluster_number_tmp):
        names_holder.append(final_variant_names[ee])

mode_val = statistics.mode(names_holder)
# print("Cluster Name = ",mode_val, ", Cluster ID = ",qq)

#     print("Cluster = ",clusters_ids)
idx = pd.Index(names_holder) # creates an index which allows counting the entries easily
print(len(idx),"entries in total")
aq = (idx.value_counts())
print(aq,"\n")

#############  Confusion Matrix ################

unique_labels = np.unique(names_holder)
labels_integers = []
for i in range(len(names_holder)):
    aa = (np.where(unique_labels==names_holder[i]))
    labels_integers.append(aa[0][0])

pred_labels = []
for ch in range(len(names_holder)):
    pred_labels.append(cluster_name_tmp)

confuse = confusion_matrix(names_holder,pred_labels)

confuse_matrix = []
for i in range(0,len(confuse)):
    confuse_matrix.append(list(confuse[i]))
#############  Confusion Matrix ################

print("F1 macro",f1_score(names_holder, pred_labels, average='macro'))
print("F1 micro",f1_score(names_holder, pred_labels, average='micro'))
print("F1 weighted",f1_score(names_holder, pred_labels, average='weighted'))


# # Classification Approach

# In[245]:


final_freq_vec = []
flat_list = [item for sublist in protein_sequences_final for item in sublist]

uniqueWords = np.unique(flat_list)
uniqueWords_2 = list(uniqueWords)


# In[249]:


uniqueWords


# In[250]:


for ind_freq in range(len(protein_sequences_final)):
    print(ind_freq,"/",len(protein_sequences_final))
    temp_freq_vec = [0]* len(uniqueWords)
    aaa = protein_sequences_final[ind_freq]
    for i in range(len(aaa)):
        ind = uniqueWords_2.index(aaa[i])
        temp_freq_vec[ind] = temp_freq_vec[ind] + 1
    final_freq_vec.append(temp_freq_vec)


# In[265]:


X_train, X_test, y_train, y_test = train_test_split(
    np.array(final_jernel_mat),final_variant_names,
    test_size=0.3,
    random_state=1,
shuffle=True)
# final_freq_vec


# In[ ]:


svm_return = svm_fun(X_train,y_train,X_test,y_test)


# In[ ]:


svm_table_final_no_DL = DataFrame(svm_table_no_DL, columns=["Accuracy","Precision","Recall",
                                                "F1 (weighted)","F1 (Macro)","F1 (Micro)","ROC AUC"])


# In[266]:


##########################  SVM Classifier  ################################
def svm_fun(X_train,y_train,X_test,y_test):
    #Create a svm Classifier
    clf = svm.SVC(kernel='linear') # Linear Kernel

    #Train the model using the training sets
    clf.fit(X_train, y_train)

    #Predict the response for test dataset
    y_pred = clf.predict(X_test)
    
    svm_acc = metrics.accuracy_score(y_test, y_pred)
#     print("SVM Accuracy:",svm_acc)
    
    svm_prec = metrics.precision_score(y_test, y_pred,average='weighted')
#     print("SVM Precision:",svm_prec)
    
    svm_recall = metrics.recall_score(y_test, y_pred,average='weighted')
#     print("SVM Recall:",svm_recall)

    svm_f1_weighted = metrics.f1_score(y_test, y_pred,average='weighted')
#     print("SVM F1 Weighted:",svm_f1_weighted)
    
    svm_f1_macro = metrics.f1_score(y_test, y_pred,average='macro')
#     print("SVM F1 macro:",svm_f1_macro)
    
    svm_f1_micro = metrics.f1_score(y_test, y_pred,average='micro')
#     print("SVM F1 micro:",svm_f1_micro)
    
    confuse = confusion_matrix(y_test, y_pred)
    print("Confusion Matrix SVM : \n", confuse)
    print("SVM Class Wise Accuracy : ",confuse.diagonal()/confuse.sum(axis=1))
    ######################## Compute ROC curve and ROC area for each class ################
    y_prob = y_pred
    macro_roc_auc_ovo = roc_auc_score_multiclass(y_test, y_prob, average='macro')
#    print(macro_roc_auc_ovo[1])
    check = [svm_acc,svm_prec,svm_recall,svm_f1_weighted,svm_f1_macro,svm_f1_micro,macro_roc_auc_ovo[1]]
    return(check)


# In[261]:


top_k_val = 3
predicted_labels = []

for test_ind in range(len(X_test)):
    temp_dot_prod = []
    for i in range(len(X_train)):
        temp_dot_prod.append(np.dot(X_test[test_ind],X_train[i]))

    top_2_idx = np.argsort(temp_dot_prod)[-top_k_val:] #index of top k values
#     top_2_idx = np.argsort(temp_dot_prod)[:top_k_val] #index of smallest k values
    top_2_values = [temp_dot_prod[i] for i in top_2_idx] #top k values

    train_val_check = [y_train[i] for i in top_2_idx] #label index of top k 
    final_label_check = statistics.mode(train_val_check) #mode label for top k
    predicted_labels.append(final_label_check)


# In[262]:


# In[4]
def roc_auc_score_multiclass(actual_class, pred_class, average = "macro"):
    #creating a set of all the unique classes using the actual class list
    unique_class = set(actual_class)
    roc_auc_dict = {}
    for per_class in unique_class:
        #creating a list of all the classes except the current class 
        other_class = [x for x in unique_class if x != per_class]

        #marking the current class as 1 and all other classes as 0
        new_actual_class = [0 if x in other_class else 1 for x in actual_class]
        new_pred_class = [0 if x in other_class else 1 for x in pred_class]

        #using the sklearn metrics method to calculate the roc_auc_score
        roc_auc = roc_auc_score(new_actual_class, new_pred_class, average = average)
        roc_auc_dict[per_class] = roc_auc


    check = pd.DataFrame(roc_auc_dict.items())
    return mean(check)


# In[263]:


LR_acc = metrics.accuracy_score(y_test, predicted_labels)
LR_prec = metrics.precision_score(y_test, predicted_labels,average='weighted')
LR_recall = metrics.recall_score(y_test, predicted_labels,average='weighted')
LR_f1_weighted = metrics.f1_score(y_test, predicted_labels,average='weighted')
LR_f1_macro = metrics.f1_score(y_test, predicted_labels,average='macro')
LR_f1_micro = metrics.f1_score(y_test, predicted_labels,average='micro')
macro_roc_auc_ovo = roc_auc_score_multiclass(y_test, predicted_labels, average='macro')

confuse = confusion_matrix(y_test, predicted_labels)
print("Confusion Matrix LR : \n", confuse)
print("LR Class Wise Accuracy : ",confuse.diagonal()/confuse.sum(axis=1))


# In[264]:


check = [LR_acc,LR_prec,LR_recall,LR_f1_weighted,LR_f1_macro,LR_f1_micro,macro_roc_auc_ovo[1]]

aa = DataFrame(check)
aaa = np.transpose(aa)

final_mat = DataFrame(np.transpose(aa))
final_mat.columns =["Accuracy","Precision","Recall",
                                                "F1 (weighted)","F1 (Macro)","F1 (Micro)","ROC AUC"]
final_mat


# In[ ]:





# In[ ]:




