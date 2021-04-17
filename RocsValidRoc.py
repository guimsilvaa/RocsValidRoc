import subprocess
import shutil
import os
import sys
import smtplib
import pandas as pd
import matplotlib.pyplot as plt
from numpy import argmax
from time import sleep
from fnmatch import fnmatch
from datetime import datetime
from sklearn.metrics import *

def leiaInt(msg):
    while True:
        try:
            n = int(input(msg))
        except (ValueError, TypeError):
            print('  o.O  this option does not exists. Try again!')
            continue
        except (KeyboardInterrupt):
            print('User prefered not to type a thing.. what?!')
            return 0
        else:
            return n

def linha(tam = 80):
    return '~' * tam

def cabeçalho(txt):
    print(linha())
    print(txt.center(80))
    print(linha())

def menu(lista):
    cabeçalho('R o c s V a l i d R o c')
    c = 1
    for item in lista:
        print('{} - {}'.format(c, item))
        c += 1
    print(linha())
    opc = leiaInt('\nYour option: ')
    return opc
    print('\n')

def findfiles_ByExtension(path, extension):
    filenames = os.listdir(path)
    return [filename for filename in filenames if filename.endswith(extension)]

def findfiles_ByPrefix(path, prefix):
    filenames = os.listdir(path)
    return [filename for filename in filenames if filename.startswith(prefix)]

def WalkDirs_ListPattern(pattern1, pattern2, pattern3, path):
    c = 1
    return [name for root, dirs, files in os.walk(os.path.join(path))
            for name in files if
            fnmatch(name, pattern1) or
            fnmatch(name, pattern2) or
            fnmatch(name, pattern3)]

def CreateDir(name_dir, parent_dir):
    path_out = os.path.join(parent_dir, name_dir)
    return [os.makedirs(i, exist_ok=True) for i in [path_out]]

############################ MENU ##############################
# Menu options
while True:
    resposta = menu(['Help', 'Check Databases', 'Check Queries', 'Run RocsValidRoc!', 'Run only ROCS', 'Build single OR multi query ROC curve from ROCS results', 'Exit'])

############################ Answer for option 1 ##############################
# option 1 opens HELP
    if resposta == 1:
        print('\n~~ Help\n'
              '\n~~ What does this script do?\n'
              '\nRocsValidRoc fully automatize in an interactive way the validation procedure of using one (or various) molecule(s) as query(ies) in 3D shape overlaping (similarity search) using ROCS, within databases of known actives and decoys. For this, it automatically builds ROC curve(s) allowing you to check its metrics and to assess your query(ies) performance in ROCS.' 
              '\nIt does the same as vROCS, however in a more efficient manner, since it runs in command line (not in GUI). Moreover, it allows one to use a stronger computer processing power when setting the desired number of threads in MPI (mpi_np). In addition, vROCS only allows you to validate 1 query at a time, here we can appreciate RocsValidRoc function of automatically generating merged ROC curves to simultaneously validate multi-queries (with much more beautiful ROC graphics, let´s face it!).' 
              '\nReminder: validation (performance assessment) is a must when it comes to check the efficiency of a methodology such as ROCS in a VS campaign.'
              '\nRocsValidRoc is strongly recommended to be used in conjunction with, for instance, DUD-E databases http://dude.docking.org/ You might want to check their ready to run active/decoys databases for several targets.'
              '\nFurthermore, it is very useful if you only wish to use several databases for a given query(ies) in a ROCS VS campaign. You will not have to type every single code line to run ROCS, for each database, as it is usually required in traditional procedure (for this, try ´Run only ROCS´).\n' 
              '\n~~ What do you need?\n'
              '\n> OpenEye applications installed and regularly working in your machine (including, obviously, ROCS), with appropriate OpenEye valid license (see https://docs.eyesopen.com/ for more info)'
              '\n> Python modules: subprocess, shutil, os, sys, smtplib, pandas,  matplotlib, numpy, time, fnmatch, datetime, scikit-learn. In general, they all can be easily installed by typing pip install.\n'
              '\n~~ Place in folder:\n'
              '\n>> RocsValidRoc.py script'
              '\n>> 1 ´actives´ and 1 ´decoys´ databases files previously processed by OMEGA conformer generation that you wish to run ROCS validation. We highly recommend using default names for databases: [omega]actives.oeb.gz and [omega]decoys.oeb.gz'
              '\n>> query file (in sdf or mol2 or sq). If you wish to use more than 1 query, you should merge them all into one single molecule file (in sdf or mol2 or sq) and inform it as the chosen query. Please note that maximum of 10 multi-queries (within 1 query file) is allowed in this current version.'
              '\n* if you´re running only ´Build single OR multi query ROC curve from ROCS results´ you should place this script in folder containing ROCS results, specially .rpt files, and then run it.\n'
              '\n~~ Some keynotes:\n'
              '\n>>> Avoid running it in folders (path or directories) named with SPACES on it!'
              '\n>>> RocsValidRoc will only list DATABASEs files containing the terms *actives* or *decoys* in corresponding filenames, within the current working directory (and sub-directories)! Moreover, it should contain term in brackets, e.g. [omega]'
              '\n>>> If you´re using multi-QUERY (into one single molecule file, in sdf or mol2), results for each query will be in same order as molecules in such file, i.e. first molecule in query-file will be hits_1 and first curve, second hits_2 and second curve, so on... ROC curve output txt file will also help you to track your queries.'
              '\n>>> RocsValidRoc will create a sub-directory that you set as your project name at your current working directory, and it´ll send all results there. In addition, it will generate time.txt file with start/end of ROCS run, plus output_rocsvalidroc.txt with all that you wanna know about your built ROC curves.'
              '\n>>> Besides AUC values, we provide the optimal threshold in a given ROC curve, which is calculated following https://en.wikipedia.org/wiki/Youden%27s_J_statistic'
              '\n>>> RocsValidRoc must be used with databases previously prepared by OMEGA conformer generation.'
              '\n>>> Also, it is compatible to be used after a regular ROCS run (for this, try ´Build single OR multi query ROC curve from ROCS results´). However, RocsValidRoc already comes with ROCS function implemented. ** Just remember, if you use ROCS separetely/previous to RocsValidRoc, name your OMEGA processed databases with terms ´active´ and ´decoys´ and then you may proceed to run RocsValidRoc **\n'
              '\n~~ Author: Guilherme M. Silva (silvagm@usp.br)'
              '\n   *adapted from https://github.com/lbfederico/OpenEye/tree/main/RocsEon_Run by Leonardo B. Federico'
              '\n   Computational Laboratory of Pharmaceutical Chemistry (LCQF)' 
              '\n   FCFRP, University of São Paulo, Brazil'
              '\n   https://sites.usp.br/lcqf/ \n'
              '\n Returning to menu... \n')
                

############################ Answer for option 2 ##############################
# answer 2 list databases
# search for keywords active or inactive or decoys matching within filenames present in directory and list´em
    elif resposta == 2:
        parent_dir = os.getcwd()
        name1 = WalkDirs_ListPattern('*active*', '*inactive*', '*decoy*', os.path.join(parent_dir) )
        c = 1
        for lst in name1:
            print(f'{c} - {lst}')
            c += 1
            sleep(0.01)
        print('\n Returning to menu...\n')
                

############################ Answer for option 3 ##############################
# answer 3 list queries
# search files mol2 sdf present in directory and list´em
    elif resposta == 3: 
        parent_dir = os.getcwd()
        name2 = WalkDirs_ListPattern('*.mol2*', '*.sdf*', '*.sq*', os.path.join(parent_dir) )
        c = 1
        for lst in name2:
            print(f'{c} - {lst}')
            c += 1
            sleep(0.01)
        print('\n Returning to menu...\n')
        
############################ Answer for option 4 ##############################
# let´s AUTO RUN ROCSVALIDROC
    elif resposta == 4:    
        cabeçalho('Yes! Let´s run RocsValidRoc all at once!')        
        parent_dir = os.getcwd()
        print(str('We are working at'), parent_dir, str('\n'))
        nome_projeto = str(input('\nSet the name of your RocsValidRoc project, that u want to be printed in output folder\n>>'))
                
# DATABASES PART
        print(linha()) 
        print('Now check the available Actives,Decoys DATABASES: \n')

# list databases containing mol2 or sdf or oeb in filenames and present the list to user
# opens the list with databases names and address
        pattern1 = '*active*'
        pattern2 = '*inactive*'
        pattern3 = '*decoy*'
        parent_dir = os.getcwd()
        c = 1
        name1_base = []
        end1_base = []
        for root, dirs, files in os.walk(os.path.join(parent_dir)):
            for name1 in files:
                if fnmatch(name1, pattern1) or fnmatch(name1, pattern2) or fnmatch(name1, pattern3):
                    print(f'{c} - {name1}')
                    c += 1
                    sleep(0.01)

# gather results in name lists and addresses
                    name1_base.append(name1)
                    end1_base.append(os.path.join(root,name1))
                    #print(name1_base) # optional/if necessary

# select database within list
        print('\nSelect which Actives,Decoys DATABASES u wanna use (separate them by comma, eg, 1,2)')
        #base_proj = [] # optional/if necessary
        end_proj1 = []
        w = input()
        #print(w.split(',')) # optional/if necessary
        print(linha()) 
        print('You have selected the following DATABASES: \n')
        list = w.split(',')
        for i in list:
            x = int(i) - 1
            print(name1_base[x])
            #print(end1_base[x]) # optional/if necessary

# gather selected databases (name_base) in one list (base_proj)
            end_proj1.append(end1_base[x])
            #print(end_proj1) # optional/if necessary
                                  
# ask user to check databases
        print(linha())
        k = str(input('Is that correct? Shall we proceed? [y/n] '))
        print()
        if k == 'n':
            cabeçalho('Damn! Let´s try again! ')
            sleep(3)
            continue

        elif k == 'y':
                           
# QUERIES PART
            print(linha()) 
            print('Now let´s go to the QUERIES part. See what´s available here: \n')
            parent_dir = os.getcwd()
            name2 = WalkDirs_ListPattern('*.mol2*', '*.sdf*', '*.sq*', os.path.join(parent_dir) )
            for lst in name2:
                print(f'{lst}')
                sleep(0.01)
    
            y = os.getcwd()
            z = str(input('\nType the exact full filename of QUERY (or QUERIES if in one single molecule file) with extension (Eg: malato.mol2)\n** if your query file is not here, type findquery **\n>>'))
            print()
            if z == 'findquery':
                y = str(input('Please confirm the path of directory containing QUERY (Eg: C:\Gui\RocsValidRoc\queries)\n>>'))
                zz = str(input('\nInsert the exact full filename of QUERY with extension (Eg: malato.mol2)\n>>'))
                z = zz
            else:
                print()
                    
# stablishing other ROCS flags and running it        
            print(linha()) 
            print('Some other ROCS settings:')
            
            h = str(input('\nChoose extension/format for output results ´oformat´ (Eg: sdf, sdf.gz, oeb.gz, oeb)\n** do not use . dot before/after the extension name **\n>>'))
            mpinp = str(input('\nLast, but not least, specify the number of processors to run ROCS in MPI mode\n** if this doesn´t matter for you, just type 0 **\n>>'))
           
            print(linha())
            print('Ok, ROCS is ready to go.. now let´s just define a few more things for ROC plotting')
            sleep(2)
            print('...')
            sleep(1)
            print('......')
            print('Alrite!\nIf OpenEye logo comes up it is everything ok...'
                      '\nLet´s do this!')
            sleep(2)

# Criar diretório para os resultados do ROCS na pasta do usuário
            directory = '[OutputRocs]' + str(nome_projeto)  # Pode colocar a data
            parent_dir = os.getcwd()
            path_out_rocs = os.path.join(parent_dir, directory)
            [os.makedirs(i, exist_ok=True) for i in [path_out_rocs]]

# starting at date and time
            now = datetime.now()
            day_inicio_proces = now.strftime('%d.%m.%Y')
            hour_inicio_proces = now.strftime('%H:%M:%S')
            
# databases/query loop and further settings
            for base in end_proj1:
                
                file_name1 = os.path.basename(os.path.join(base))
                index_of_dot1 = file_name1.index('.')
                i = file_name1[:index_of_dot1]
                dbname = i.split(']')[1]
                
                queryname = z.split('.')[0]
                
# show date and time in terminal for each ROCS run
                now = datetime.now()
                day_inicio = now.strftime('%d.%m.%Y')
                hour_inicio = now.strftime('%H:%M:%S')
                print(linha())
                print()
                print('Running ROCS for database', i)
                print('on {} at {}'.format(day_inicio, hour_inicio))
                print()
                print(linha())

# create txt file to write date and time of each ROCS run
                arquivo = open('time_rocs.txt', 'a')
                arquivo.writelines(['\nRunning ROCS for database: ', i,
                                    '\n on {} at {}'.format(day_inicio, hour_inicio),
                                    '\n',
                                    ])
                arquivo.close()
# ROCS run
                subprocess.run(['rocs.bat',
                                '-dbase', os.path.join(base),
                                '-query', os.path.join(y, z),
                                '-outputdir', os.path.join(path_out_rocs),
                                '-prefix', 'OutRocs_' + dbname + '_' + queryname,
                                '-oformat', h,
                                '-besthits', '0',
                                '-mpi_np', mpinp
                                ])
# end of loop Rocs
# Final message
            now = datetime.now()
            day_final = now.strftime('%d.%m.%Y')
            hour_final = now.strftime('%H:%M:%S')
            print(linha())
            print('ROCS has finished running!')
            print('Started {} at {}'.format(day_inicio_proces, hour_inicio_proces))
            print('Ended {} at {}'.format(day_final, hour_final))
            print('Thank you for using AutomaROCS :-)')

# writting final message in txt file
            arquivo = open('time_rocs.txt', 'a')
            arquivo.writelines(['\nROCS has finished running!',
                                '\nStarted {} at {}'.format(day_inicio_proces, hour_inicio_proces),
                                '\nEnded {} at {}'.format(day_final, hour_final),
                                '\nThank you for using AutomaROCS :-)',])
            arquivo.close()

                  
                  
#
# STARTING AUTOMATIC ROC CURVE GENERATION FROM ROCS RESULTS
            print(linha())
            print('Starting ROC curve(s) generation')
            sleep(1)
            print('..')
            sleep(1)
            print('...')
            sleep(1)
            raU = os.path.join(path_out_rocs)
            os.chdir(os.path.join(raU))
            print("Current working directory: {0}".format(os.getcwd()))

# RETRIEVING RESULTS OF QUERY (1) FOR ACTIVE_1 (A1) and DECOYS_1 (D1)
            parentdirectory = os.getcwd()
            for fileA1 in os.listdir(parentdirectory):
                if fnmatch(fileA1, '*active*') and fileA1.endswith("_1.rpt"):
                        with open(fileA1) as fin, open('fileA1.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA1 = pd.read_csv("fileA1.csv") 
                        A1data = rawA1.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A1data.head()
                        A1data.insert(2, "Outcome", "1")
                        A1data.head()
                        A1data.to_csv("A1file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD1 in os.listdir(parentdirectory):
                if fnmatch(fileD1, '*decoy*') and fileD1.endswith("_1.rpt"):
                        with open(fileD1) as fin, open('fileD1.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD1 = pd.read_csv("fileD1.csv") 
                        D1data = rawD1.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D1data.head()
                        D1data.insert(2, "Outcome", "0")
                        D1data.head()
                        D1data.to_csv("D1file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA1 = str('fileA1.csv')
            existencialA1 = os.path.exists(os.path.join(parentdirectory, existA1))
            if existencialA1 is True:
                file1A1 = open("A1file.csv", "a")
                file2D1 = open("D1file.csv", "r")
                for line in file2D1:
                    file1A1.write(line)
                file1A1.close()
                file2D1.close()
                sortdataA1 = pd.read_csv("A1file.csv") 
                sortdataaA1 = sortdataA1.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA1 = sortdataaA1.drop(sortdataaA1.index[0])
                sortdataaaaA1 = sortdataaaA1.rename(columns={"TanimotoCombo": "TanimotoCombo1", "Outcome": "Outcome1"}) #attention
                sortdataaaaaA1 = sortdataaaaA1.reset_index()
                sortdataaaaaA1.insert(0, 'New_ID', 0)
                sortdataaaaaA1['New_ID'] = sortdataaaaaA1.index + 0
                sortdataaaaaA1.to_csv("A1D1MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A1file.csv")
                os.remove("fileA1.csv")
                os.remove("fileD1.csv")
                os.remove("D1file.csv")
            else:
                print(str('The existence of a query A1/D1 is'), existencialA1)

                
# RETRIEVING RESULTS OF QUERY (2) FOR ACTIVE_2 (A2) and DECOYS_2 (D2)
            parentdirectory = os.getcwd()
            for fileA2 in os.listdir(parentdirectory):
                if fnmatch(fileA2, '*active*') and fileA2.endswith("_2.rpt"):
                        with open(fileA2) as fin, open('fileA2.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA2 = pd.read_csv("fileA2.csv") 
                        A2data = rawA2.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A2data.head()
                        A2data.insert(2, "Outcome", "1")
                        A2data.head()
                        A2data.to_csv("A2file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD2 in os.listdir(parentdirectory):
                if fnmatch(fileD2, '*decoy*') and fileD2.endswith("_2.rpt"):
                        with open(fileD2) as fin, open('fileD2.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD2 = pd.read_csv("fileD2.csv") 
                        D2data = rawD2.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D2data.head()
                        D2data.insert(2, "Outcome", "0")
                        D2data.head()
                        D2data.to_csv("D2file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA2 = str('fileA2.csv')
            existencialA2 = os.path.exists(os.path.join(parentdirectory, existA2))
            if existencialA2 is True:
                file1A2 = open("A2file.csv", "a")
                file2D2 = open("D2file.csv", "r")
                for line in file2D2:
                    file1A2.write(line)
                file1A2.close()
                file2D2.close()
                sortdataA2 = pd.read_csv("A2file.csv") 
                sortdataaA2 = sortdataA2.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA2 = sortdataaA2.drop(sortdataaA2.index[0])
                sortdataaaaA2 = sortdataaaA2.rename(columns={"TanimotoCombo": "TanimotoCombo2", "Outcome": "Outcome2"}) #attention
                sortdataaaaaA2 = sortdataaaaA2.reset_index()
                sortdataaaaaA2.insert(0, 'New_ID', 0)
                sortdataaaaaA2['New_ID'] = sortdataaaaaA2.index + 0
                sortdataaaaaA2.to_csv("A2D2MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A2file.csv")
                os.remove("fileA2.csv")
                os.remove("fileD2.csv")
                os.remove("D2file.csv")
            else:
                print(str('The existence of a second query A2/D2 is'), existencialA2)

# RETRIEVING RESULTS OF QUERY (3) FOR ACTIVE_3 (A3) and DECOYS_3 (D3)
            parentdirectory = os.getcwd()
            for fileA3 in os.listdir(parentdirectory):
                if fnmatch(fileA3, '*active*') and fileA3.endswith("_3.rpt"):
                        with open(fileA3) as fin, open('fileA3.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA3 = pd.read_csv("fileA3.csv") 
                        A3data = rawA3.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A3data.head()
                        A3data.insert(2, "Outcome", "1")
                        A3data.head()
                        A3data.to_csv("A3file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD3 in os.listdir(parentdirectory):
                if fnmatch(fileD3, '*decoy*') and fileD3.endswith("_3.rpt"):
                        with open(fileD3) as fin, open('fileD3.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD3 = pd.read_csv("fileD3.csv") 
                        D3data = rawD3.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D3data.head()
                        D3data.insert(2, "Outcome", "0")
                        D3data.head()
                        D3data.to_csv("D3file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA3 = str('fileA3.csv')
            existencialA3 = os.path.exists(os.path.join(parentdirectory, existA3))
            if existencialA3 is True:
                file1A3 = open("A3file.csv", "a")
                file2D3 = open("D3file.csv", "r")
                for line in file2D3:
                    file1A3.write(line)
                file1A3.close()
                file2D3.close()
                sortdataA3 = pd.read_csv("A3file.csv") 
                sortdataaA3 = sortdataA3.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA3 = sortdataaA3.drop(sortdataaA3.index[0])
                sortdataaaaA3 = sortdataaaA3.rename(columns={"TanimotoCombo": "TanimotoCombo3", "Outcome": "Outcome3"}) #attention
                sortdataaaaaA3 = sortdataaaaA3.reset_index()
                sortdataaaaaA3.insert(0, 'New_ID', 0)
                sortdataaaaaA3['New_ID'] = sortdataaaaaA3.index + 0
                sortdataaaaaA3.to_csv("A3D3MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A3file.csv")
                os.remove("fileA3.csv")
                os.remove("fileD3.csv")
                os.remove("D3file.csv")
            else:
                print(str('The existence of a third query A3/D3 is'), existencialA3)

# RETRIEVING RESULTS OF QUERY (4) FOR ACTIVE_4 (A4) and DECOYS_4 (D4)
            parentdirectory = os.getcwd()
            for fileA4 in os.listdir(parentdirectory):
                if fnmatch(fileA4, '*active*') and fileA4.endswith("_4.rpt"):
                        with open(fileA4) as fin, open('fileA4.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA4 = pd.read_csv("fileA4.csv") 
                        A4data = rawA4.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A4data.head()
                        A4data.insert(2, "Outcome", "1")
                        A4data.head()
                        A4data.to_csv("A4file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD4 in os.listdir(parentdirectory):
                if fnmatch(fileD4, '*decoy*') and fileD4.endswith("_4.rpt"):
                        with open(fileD4) as fin, open('fileD4.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD4 = pd.read_csv("fileD4.csv") 
                        D4data = rawD4.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D4data.head()
                        D4data.insert(2, "Outcome", "0")
                        D4data.head()
                        D4data.to_csv("D4file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA4 = str('fileA4.csv')
            existencialA4 = os.path.exists(os.path.join(parentdirectory, existA4))
            if existencialA4 is True:
                file1A4 = open("A4file.csv", "a")
                file2D4 = open("D4file.csv", "r")
                for line in file2D4:
                    file1A4.write(line)
                file1A4.close()
                file2D4.close()
                sortdataA4 = pd.read_csv("A4file.csv") 
                sortdataaA4 = sortdataA4.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA4 = sortdataaA4.drop(sortdataaA4.index[0])
                sortdataaaaA4 = sortdataaaA4.rename(columns={"TanimotoCombo": "TanimotoCombo4", "Outcome": "Outcome4"}) #attention
                sortdataaaaaA4 = sortdataaaaA4.reset_index()
                sortdataaaaaA4.insert(0, 'New_ID', 0)
                sortdataaaaaA4['New_ID'] = sortdataaaaaA4.index + 0
                sortdataaaaaA4.to_csv("A4D4MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A4file.csv")
                os.remove("fileA4.csv")
                os.remove("fileD4.csv")
                os.remove("D4file.csv")
            else:
                print(str('The existence of a fourth query A4/D4 is'), existencialA4)

# RETRIEVING RESULTS OF QUERY (5) FOR ACTIVE_5 (A5) and DECOYS_5 (D5)
            parentdirectory = os.getcwd()
            for fileA5 in os.listdir(parentdirectory):
                if fnmatch(fileA5, '*active*') and fileA5.endswith("_5.rpt"):
                        with open(fileA5) as fin, open('fileA5.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA5 = pd.read_csv("fileA5.csv") 
                        A5data = rawA5.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A5data.head()
                        A5data.insert(2, "Outcome", "1")
                        A5data.head()
                        A5data.to_csv("A5file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD5 in os.listdir(parentdirectory):
                if fnmatch(fileD5, '*decoy*') and fileD5.endswith("_5.rpt"):
                        with open(fileD5) as fin, open('fileD5.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD5 = pd.read_csv("fileD5.csv") 
                        D5data = rawD5.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D5data.head()
                        D5data.insert(2, "Outcome", "0")
                        D5data.head()
                        D5data.to_csv("D5file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA5 = str('fileA5.csv')
            existencialA5 = os.path.exists(os.path.join(parentdirectory, existA5))
            if existencialA5 is True:
                file1A5 = open("A5file.csv", "a")
                file2D5 = open("D5file.csv", "r")
                for line in file2D5:
                    file1A5.write(line)
                file1A5.close()
                file2D5.close()
                sortdataA5 = pd.read_csv("A5file.csv") 
                sortdataaA5 = sortdataA5.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA5 = sortdataaA5.drop(sortdataaA5.index[0])
                sortdataaaaA5 = sortdataaaA5.rename(columns={"TanimotoCombo": "TanimotoCombo5", "Outcome": "Outcome5"}) #attention
                sortdataaaaaA5 = sortdataaaaA5.reset_index()
                sortdataaaaaA5.insert(0, 'New_ID', 0)
                sortdataaaaaA5['New_ID'] = sortdataaaaaA5.index + 0
                sortdataaaaaA5.to_csv("A5D5MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A5file.csv")
                os.remove("fileA5.csv")
                os.remove("fileD5.csv")
                os.remove("D5file.csv")
            else:
                print(str('The existence of a fifth query A5/D5 is'), existencialA5)

# RETRIEVING RESULTS OF QUERY (6) FOR ACTIVE_6 (A6) and DECOYS_6 (D6)
            parentdirectory = os.getcwd()
            for fileA6 in os.listdir(parentdirectory):
                if fnmatch(fileA6, '*active*') and fileA6.endswith("_6.rpt"):
                        with open(fileA6) as fin, open('fileA6.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA6 = pd.read_csv("fileA6.csv") 
                        A6data = rawA6.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A6data.head()
                        A6data.insert(2, "Outcome", "1")
                        A6data.head()
                        A6data.to_csv("A6file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD6 in os.listdir(parentdirectory):
                if fnmatch(fileD6, '*decoy*') and fileD6.endswith("_6.rpt"):
                        with open(fileD6) as fin, open('fileD6.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD6 = pd.read_csv("fileD6.csv") 
                        D6data = rawD6.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D6data.head()
                        D6data.insert(2, "Outcome", "0")
                        D6data.head()
                        D6data.to_csv("D6file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA6 = str('fileA6.csv')
            existencialA6 = os.path.exists(os.path.join(parentdirectory, existA6))
            if existencialA6 is True:
                file1A6 = open("A6file.csv", "a")
                file2D6 = open("D6file.csv", "r")
                for line in file2D6:
                    file1A6.write(line)
                file1A6.close()
                file2D6.close()
                sortdataA6 = pd.read_csv("A6file.csv") 
                sortdataaA6 = sortdataA6.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA6 = sortdataaA6.drop(sortdataaA6.index[0])
                sortdataaaaA6 = sortdataaaA6.rename(columns={"TanimotoCombo": "TanimotoCombo6", "Outcome": "Outcome6"}) #attention
                sortdataaaaaA6 = sortdataaaaA6.reset_index()
                sortdataaaaaA6.insert(0, 'New_ID', 0)
                sortdataaaaaA6['New_ID'] = sortdataaaaaA6.index + 0
                sortdataaaaaA6.to_csv("A6D6MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A6file.csv")
                os.remove("fileA6.csv")
                os.remove("fileD6.csv")
                os.remove("D6file.csv")
            else:
                print(str('The existence of a sixth query A6/D6 is'), existencialA6)

# RETRIEVING RESULTS OF QUERY (7) FOR ACTIVE_7 (A7) and DECOYS_7 (D7)
            parentdirectory = os.getcwd()
            for fileA7 in os.listdir(parentdirectory):
                if fnmatch(fileA7, '*active*') and fileA7.endswith("_7.rpt"):
                        with open(fileA7) as fin, open('fileA7.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA7 = pd.read_csv("fileA7.csv") 
                        A7data = rawA7.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A7data.head()
                        A7data.insert(2, "Outcome", "1")
                        A7data.head()
                        A7data.to_csv("A7file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD7 in os.listdir(parentdirectory):
                if fnmatch(fileD7, '*decoy*') and fileD7.endswith("_7.rpt"):
                        with open(fileD7) as fin, open('fileD7.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD7 = pd.read_csv("fileD7.csv") 
                        D7data = rawD7.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D7data.head()
                        D7data.insert(2, "Outcome", "0")
                        D7data.head()
                        D7data.to_csv("D7file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA7 = str('fileA7.csv')
            existencialA7 = os.path.exists(os.path.join(parentdirectory, existA7))
            if existencialA7 is True:
                file1A7 = open("A7file.csv", "a")
                file2D7 = open("D7file.csv", "r")
                for line in file2D7:
                    file1A7.write(line)
                file1A7.close()
                file2D7.close()
                sortdataA7 = pd.read_csv("A7file.csv") 
                sortdataaA7 = sortdataA7.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA7 = sortdataaA7.drop(sortdataaA7.index[0])
                sortdataaaaA7 = sortdataaaA7.rename(columns={"TanimotoCombo": "TanimotoCombo7", "Outcome": "Outcome7"}) #attention
                sortdataaaaaA7 = sortdataaaaA7.reset_index()
                sortdataaaaaA7.insert(0, 'New_ID', 0)
                sortdataaaaaA7['New_ID'] = sortdataaaaaA7.index + 0
                sortdataaaaaA7.to_csv("A7D7MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A7file.csv")
                os.remove("fileA7.csv")
                os.remove("fileD7.csv")
                os.remove("D7file.csv")
            else:
                print(str('The existence of a seventh query A7/D7 is'), existencialA7)

# RETRIEVING RESULTS OF QUERY (8) FOR ACTIVE_8 (A8) and DECOYS_8 (D8)
            parentdirectory = os.getcwd()
            for fileA8 in os.listdir(parentdirectory):
                if fnmatch(fileA8, '*active*') and fileA8.endswith("_8.rpt"):
                        with open(fileA8) as fin, open('fileA8.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA8 = pd.read_csv("fileA8.csv") 
                        A8data = rawA8.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A8data.head()
                        A8data.insert(2, "Outcome", "1")
                        A8data.head()
                        A8data.to_csv("A8file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD8 in os.listdir(parentdirectory):
                if fnmatch(fileD8, '*decoy*') and fileD8.endswith("_8.rpt"):
                        with open(fileD8) as fin, open('fileD8.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD8 = pd.read_csv("fileD8.csv") 
                        D8data = rawD8.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D8data.head()
                        D8data.insert(2, "Outcome", "0")
                        D8data.head()
                        D8data.to_csv("D8file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA8 = str('fileA8.csv')
            existencialA8 = os.path.exists(os.path.join(parentdirectory, existA8))
            if existencialA8 is True:
                file1A8 = open("A8file.csv", "a")
                file2D8 = open("D8file.csv", "r")
                for line in file2D8:
                    file1A8.write(line)
                file1A8.close()
                file2D8.close()
                sortdataA8 = pd.read_csv("A8file.csv") 
                sortdataaA8 = sortdataA8.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA8 = sortdataaA8.drop(sortdataaA8.index[0])
                sortdataaaaA8 = sortdataaaA8.rename(columns={"TanimotoCombo": "TanimotoCombo8", "Outcome": "Outcome8"}) #attention
                sortdataaaaaA8 = sortdataaaaA8.reset_index()
                sortdataaaaaA8.insert(0, 'New_ID', 0)
                sortdataaaaaA8['New_ID'] = sortdataaaaaA8.index + 0
                sortdataaaaaA8.to_csv("A8D8MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A8file.csv")
                os.remove("fileA8.csv")
                os.remove("fileD8.csv")
                os.remove("D8file.csv")
            else:
                print(str('The existence of a eigth query A8/D8 is'), existencialA8)

# RETRIEVING RESULTS OF QUERY (9) FOR ACTIVE_9 (A9) and DECOYS_9 (D9)
            parentdirectory = os.getcwd()
            for fileA9 in os.listdir(parentdirectory):
                if fnmatch(fileA9, '*active*') and fileA9.endswith("_9.rpt"):
                        with open(fileA9) as fin, open('fileA9.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA9 = pd.read_csv("fileA9.csv") 
                        A9data = rawA9.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A9data.head()
                        A9data.insert(2, "Outcome", "1")
                        A9data.head()
                        A9data.to_csv("A9file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD9 in os.listdir(parentdirectory):
                if fnmatch(fileD9, '*decoy*') and fileD9.endswith("_9.rpt"):
                        with open(fileD9) as fin, open('fileD9.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD9 = pd.read_csv("fileD9.csv") 
                        D9data = rawD9.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D9data.head()
                        D9data.insert(2, "Outcome", "0")
                        D9data.head()
                        D9data.to_csv("D9file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA9 = str('fileA9.csv')
            existencialA9 = os.path.exists(os.path.join(parentdirectory, existA9))
            if existencialA9 is True:
                file1A9 = open("A9file.csv", "a")
                file2D9 = open("D9file.csv", "r")
                for line in file2D9:
                    file1A9.write(line)
                file1A9.close()
                file2D9.close()
                sortdataA9 = pd.read_csv("A9file.csv") 
                sortdataaA9 = sortdataA9.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA9 = sortdataaA9.drop(sortdataaA9.index[0])
                sortdataaaaA9 = sortdataaaA9.rename(columns={"TanimotoCombo": "TanimotoCombo9", "Outcome": "Outcome9"}) #attention
                sortdataaaaaA9 = sortdataaaaA9.reset_index()
                sortdataaaaaA9.insert(0, 'New_ID', 0)
                sortdataaaaaA9['New_ID'] = sortdataaaaaA9.index + 0
                sortdataaaaaA9.to_csv("A9D9MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A9file.csv")
                os.remove("fileA9.csv")
                os.remove("fileD9.csv")
                os.remove("D9file.csv")
            else:
                print(str('The existence of a ninth query A9/D9 is'), existencialA9)

# RETRIEVING RESULTS OF QUERY (10) FOR ACTIVE_10 (A10) and DECOYS_10 (D10)
            parentdirectory = os.getcwd()
            for fileA10 in os.listdir(parentdirectory):
                if fnmatch(fileA10, '*active*') and fileA10.endswith("_10.rpt"):
                        with open(fileA10) as fin, open('fileA10.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA10 = pd.read_csv("fileA10.csv") 
                        A10data = rawA10.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A10data.head()
                        A10data.insert(2, "Outcome", "1")
                        A10data.head()
                        A10data.to_csv("A10file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD10 in os.listdir(parentdirectory):
                if fnmatch(fileD10, '*decoy*') and fileD10.endswith("_10.rpt"):
                        with open(fileD10) as fin, open('fileD10.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD10 = pd.read_csv("fileD10.csv") 
                        D10data = rawD10.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D10data.head()
                        D10data.insert(2, "Outcome", "0")
                        D10data.head()
                        D10data.to_csv("D10file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA10 = str('fileA10.csv')
            existencialA10 = os.path.exists(os.path.join(parentdirectory, existA10))
            if existencialA10 is True:
                file1A10 = open("A10file.csv", "a")
                file2D10 = open("D10file.csv", "r")
                for line in file2D10:
                    file1A10.write(line)
                file1A10.close()
                file2D10.close()
                sortdataA10 = pd.read_csv("A10file.csv") 
                sortdataaA10 = sortdataA10.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA10 = sortdataaA10.drop(sortdataaA10.index[0])
                sortdataaaaA10 = sortdataaaA10.rename(columns={"TanimotoCombo": "TanimotoCombo10", "Outcome": "Outcome10"}) #attention
                sortdataaaaaA10 = sortdataaaaA10.reset_index()
                sortdataaaaaA10.insert(0, 'New_ID', 0)
                sortdataaaaaA10['New_ID'] = sortdataaaaaA10.index + 0
                sortdataaaaaA10.to_csv("A10D10MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A10file.csv")
                os.remove("fileA10.csv")
                os.remove("fileD10.csv")
                os.remove("D10file.csv")
            else:
                print(str('The existence of a tenth query A10/D10 is'), existencialA10)



# merging tables of different queries
            parentdirectory = os.getcwd()
            gerouA1 = str('A1D1MergedActInactTCOutcome.csv')
            gerouA2 = str('A2D2MergedActInactTCOutcome.csv')
            gerouA3 = str('A3D3MergedActInactTCOutcome.csv')
            gerouA4 = str('A4D4MergedActInactTCOutcome.csv')
            gerouA5 = str('A5D5MergedActInactTCOutcome.csv')
            gerouA6 = str('A6D6MergedActInactTCOutcome.csv')
            gerouA7 = str('A7D7MergedActInactTCOutcome.csv')
            gerouA8 = str('A8D8MergedActInactTCOutcome.csv')
            gerouA9 = str('A9D9MergedActInactTCOutcome.csv')
            gerouA10 = str('A10D10MergedActInactTCOutcome.csv')
            existegerA1 = os.path.exists(os.path.join(parentdirectory, gerouA1))
            existegerA2 = os.path.exists(os.path.join(parentdirectory, gerouA2))
            existegerA3 = os.path.exists(os.path.join(parentdirectory, gerouA3))
            existegerA4 = os.path.exists(os.path.join(parentdirectory, gerouA4))
            existegerA5 = os.path.exists(os.path.join(parentdirectory, gerouA5))
            existegerA6 = os.path.exists(os.path.join(parentdirectory, gerouA6))
            existegerA7 = os.path.exists(os.path.join(parentdirectory, gerouA7))
            existegerA8 = os.path.exists(os.path.join(parentdirectory, gerouA8))
            existegerA9 = os.path.exists(os.path.join(parentdirectory, gerouA9))
            existegerA10 = os.path.exists(os.path.join(parentdirectory, gerouA10))
            if existegerA1 and existegerA2 is True:
                dfmA1 = pd.read_csv("A1D1MergedActInactTCOutcome.csv")
                dfmA2 = pd.read_csv("A2D2MergedActInactTCOutcome.csv")
                dfmA1A2 = (dfmA2.merge(dfmA1, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                #print(dfmA1A2)
                dfmA1A2data = dfmA1A2.drop(['index', 'index_'], axis=1)
                dfmA1A2data.to_csv("dfmA1A2.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('...'))

            gerouA1A2 = str('dfmA1A2.csv')
            existA1A2 = os.path.exists(os.path.join(parentdirectory, gerouA1A2))
            if existA1A2 and existegerA3 is True:
                df2mA1A2 = pd.read_csv("dfmA1A2.csv")
                dfmA3 = pd.read_csv("A3D3MergedActInactTCOutcome.csv")
                dfmA1A2A3 = (dfmA3.merge(df2mA1A2, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3.to_csv("dfmA1A2A3.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('....'))
                    
            gerouA1A2A3 = str('dfmA1A2A3.csv')
            existA1A2A3 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3))
            if existA1A2A3 and existegerA4 is True:
                df2mA1A2A3 = pd.read_csv("dfmA1A2A3.csv")
                dfmA4 = pd.read_csv("A4D4MergedActInactTCOutcome.csv")
                dfmA1A2A3A4 = (dfmA4.merge(df2mA1A2A3, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4.to_csv("dfmA1A2A3A4.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('.....'))

            gerouA1A2A3A4 = str('dfmA1A2A3A4.csv')
            existA1A2A3A4 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4))
            if existA1A2A3A4 and existegerA5 is True:
                df2mA1A2A3A4 = pd.read_csv("dfmA1A2A3A4.csv")
                dfmA5 = pd.read_csv("A5D5MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5 = (dfmA5.merge(df2mA1A2A3A4, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5.to_csv("dfmA1A2A3A4A5.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('......'))   
                
            gerouA1A2A3A4A5 = str('dfmA1A2A3A4A5.csv')
            existA1A2A3A4A5 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5))
            if existA1A2A3A4A5 and existegerA6 is True:
                df2mA1A2A3A4A5 = pd.read_csv("dfmA1A2A3A4A5.csv")
                dfmA6 = pd.read_csv("A6D6MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6 = (dfmA6.merge(df2mA1A2A3A4A5, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6.to_csv("dfmA1A2A3A4A5A6.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('.......'))
                
            gerouA1A2A3A4A5A6 = str('dfmA1A2A3A4A5A6.csv')
            existA1A2A3A4A5A6 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6))
            if existA1A2A3A4A5A6 and existegerA7 is True:
                df2mA1A2A3A4A5A6 = pd.read_csv("dfmA1A2A3A4A5A6.csv")
                dfmA7 = pd.read_csv("A7D7MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6A7 = (dfmA7.merge(df2mA1A2A3A4A5A6, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6A7.to_csv("dfmA1A2A3A4A5A6A7.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('........'))
                
            gerouA1A2A3A4A5A6A7 = str('dfmA1A2A3A4A5A6A7.csv')
            existA1A2A3A4A5A6A7 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6A7))
            if existA1A2A3A4A5A6A7 and existegerA8 is True:
                df2mA1A2A3A4A5A6A7 = pd.read_csv("dfmA1A2A3A4A5A6A7.csv")
                dfmA8 = pd.read_csv("A8D8MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6A7A8 = (dfmA8.merge(df2mA1A2A3A4A5A6A7, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6A7A8.to_csv("dfmA1A2A3A4A5A6A7A8.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('.........'))
                
            gerouA1A2A3A4A5A6A7A8 = str('dfmA1A2A3A4A5A6A7A8.csv')
            existA1A2A3A4A5A6A7A8 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6A7A8))
            if existA1A2A3A4A5A6A7A8 and existegerA9 is True:
                df2mA1A2A3A4A5A6A7A8 = pd.read_csv("dfmA1A2A3A4A5A6A7A8.csv")
                dfmA9 = pd.read_csv("A9D9MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6A7A8A9 = (dfmA9.merge(df2mA1A2A3A4A5A6A7A8, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6A7A8A9.to_csv("dfmA1A2A3A4A5A6A7A8A9.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('..........'))
                
            gerouA1A2A3A4A5A6A7A8A9 = str('dfmA1A2A3A4A5A6A7A8A9.csv')
            existA1A2A3A4A5A6A7A8A9 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6A7A8A9))
            if existA1A2A3A4A5A6A7A8A9 and existegerA10 is True:
                df2mA1A2A3A4A5A6A7A8A9 = pd.read_csv("dfmA1A2A3A4A5A6A7A8A9.csv")
                dfmA10 = pd.read_csv("A10D10MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6A7A8A9A10 = (dfmA10.merge(df2mA1A2A3A4A5A6A7A8A9, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6A7A8A9A10.to_csv("dfmA1A2A3A4A5A6A7A8A9A10.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('...........'))
    
# final adjustments and ploting ROC curve(S)
            gerouA1A2A3A4A5A6A7A8A9A10 = str('dfmA1A2A3A4A5A6A7A8A9A10.csv')
            existA1A2A3A4A5A6A7A8A9A10 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6A7A8A9A10))
            if existA1A2A3A4A5A6A7A8A9A10 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6A7A8A9A10.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                fpr7, tpr7, thresholds7 = roc_curve(rocdata['Outcome7'], rocdata['TanimotoCombo7'])
                fpr8, tpr8, thresholds8 = roc_curve(rocdata['Outcome8'], rocdata['TanimotoCombo8'])
                fpr9, tpr9, thresholds9 = roc_curve(rocdata['Outcome9'], rocdata['TanimotoCombo9'])
                fpr10, tpr10, thresholds10 = roc_curve(rocdata['Outcome10'], rocdata['TanimotoCombo10'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                auc7 = auc(fpr7,tpr7)
                auc8 = auc(fpr8,tpr8)
                auc9 = auc(fpr9,tpr9)
                auc10 = auc(fpr10,tpr10)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                J7 = tpr7 - fpr7
                ix7 = argmax(J7)
                best_thresh7 = thresholds7[ix7]
                J8 = tpr8 - fpr8
                ix8 = argmax(J8)
                best_thresh8 = thresholds8[ix8]
                J9 = tpr9 - fpr9
                ix9 = argmax(J9)
                best_thresh9 = thresholds9[ix9]
                J10 = tpr10 - fpr10
                ix10 = argmax(J10)
                best_thresh10 = thresholds10[ix10]                
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                plt.plot(fpr7, tpr7, color='green',
             lw=lw, label='AUC = %0.3f' % auc7)
                plt.plot(fpr8, tpr8, color='saddlebrown',
             lw=lw, label='AUC = %0.3f' % auc8)
                plt.plot(fpr9, tpr9, color='darkblue',
             lw=lw, label='AUC = %0.3f' % auc9)
                plt.plot(fpr10, tpr10, color='darkviolet',
             lw=lw, label='AUC = %0.3f' % auc10)
                # retrieving query mol2 info to report results
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 10 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n7 >> query/curve in green, with AUC = {} and best threshold = {}'.format(auc7, best_thresh7),
                                          '\n8 >> query/curve in brown, with AUC = {} and best threshold = {}'.format(auc8, best_thresh8),
                                          '\n9 >> query/curve in darkblue, with AUC = {} and best threshold = {}'.format(auc9, best_thresh9),
                                          '\n10 >> query/curve in purple, with AUC = {} and best threshold = {}'.format(auc10, best_thresh10),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                #cleaning
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("A7D7MergedActInactTCOutcome.csv")
                os.remove("A8D8MergedActInactTCOutcome.csv")
                os.remove("A9D9MergedActInactTCOutcome.csv")
                os.remove("A10D10MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") 
                os.remove("dfmA1A2A3A4A5A6.csv") 
                os.remove("dfmA1A2A3A4A5A6A7.csv") 
                os.remove("dfmA1A2A3A4A5A6A7A8.csv")
                os.remove("dfmA1A2A3A4A5A6A7A8A9.csv")
                os.remove("dfmA1A2A3A4A5A6A7A8A9A10.csv") # add a # at the beginning if u wanna keep these files
                
            elif existA1A2A3A4A5A6A7A8A9 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6A7A8A9.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                fpr7, tpr7, thresholds7 = roc_curve(rocdata['Outcome7'], rocdata['TanimotoCombo7'])
                fpr8, tpr8, thresholds8 = roc_curve(rocdata['Outcome8'], rocdata['TanimotoCombo8'])
                fpr9, tpr9, thresholds9 = roc_curve(rocdata['Outcome9'], rocdata['TanimotoCombo9'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                auc7 = auc(fpr7,tpr7)
                auc8 = auc(fpr8,tpr8)
                auc9 = auc(fpr9,tpr9)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                J7 = tpr7 - fpr7
                ix7 = argmax(J7)
                best_thresh7 = thresholds7[ix7]
                J8 = tpr8 - fpr8
                ix8 = argmax(J8)
                best_thresh8 = thresholds8[ix8]
                J9 = tpr9 - fpr9
                ix9 = argmax(J9)
                best_thresh9 = thresholds9[ix9]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                plt.plot(fpr7, tpr7, color='green',
             lw=lw, label='AUC = %0.3f' % auc7)
                plt.plot(fpr8, tpr8, color='saddlebrown',
             lw=lw, label='AUC = %0.3f' % auc8)
                plt.plot(fpr9, tpr9, color='darkblue',
             lw=lw, label='AUC = %0.3f' % auc9)
                
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 9 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n7 >> query/curve in green, with AUC = {} and best threshold = {}'.format(auc7, best_thresh7),
                                          '\n8 >> query/curve in brown, with AUC = {} and best threshold = {}'.format(auc8, best_thresh8),
                                          '\n9 >> query/curve in darkblue, with AUC = {} and best threshold = {}'.format(auc9, best_thresh9),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("A7D7MergedActInactTCOutcome.csv")
                os.remove("A8D8MergedActInactTCOutcome.csv")
                os.remove("A9D9MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") 
                os.remove("dfmA1A2A3A4A5A6.csv") 
                os.remove("dfmA1A2A3A4A5A6A7.csv") 
                os.remove("dfmA1A2A3A4A5A6A7A8.csv")
                os.remove("dfmA1A2A3A4A5A6A7A8A9.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4A5A6A7A8 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6A7A8.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                fpr7, tpr7, thresholds7 = roc_curve(rocdata['Outcome7'], rocdata['TanimotoCombo7'])
                fpr8, tpr8, thresholds8 = roc_curve(rocdata['Outcome8'], rocdata['TanimotoCombo8'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                auc7 = auc(fpr7,tpr7)
                auc8 = auc(fpr8,tpr8)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                J7 = tpr7 - fpr7
                ix7 = argmax(J7)
                best_thresh7 = thresholds7[ix7]
                J8 = tpr8 - fpr8
                ix8 = argmax(J8)
                best_thresh8 = thresholds8[ix8]
                
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                plt.plot(fpr7, tpr7, color='green',
             lw=lw, label='AUC = %0.3f' % auc7)
                plt.plot(fpr8, tpr8, color='saddlebrown',
             lw=lw, label='AUC = %0.3f' % auc8)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 8 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n7 >> query/curve in green, with AUC = {} and best threshold = {}'.format(auc7, best_thresh7),
                                          '\n8 >> query/curve in brown, with AUC = {} and best threshold = {}'.format(auc8, best_thresh8),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("A7D7MergedActInactTCOutcome.csv")
                os.remove("A8D8MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") 
                os.remove("dfmA1A2A3A4A5A6.csv") 
                os.remove("dfmA1A2A3A4A5A6A7.csv") 
                os.remove("dfmA1A2A3A4A5A6A7A8.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4A5A6A7 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6A7.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                fpr7, tpr7, thresholds7 = roc_curve(rocdata['Outcome7'], rocdata['TanimotoCombo7'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                auc7 = auc(fpr7,tpr7)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                J7 = tpr7 - fpr7
                ix7 = argmax(J7)
                best_thresh7 = thresholds7[ix7]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                plt.plot(fpr7, tpr7, color='green',
             lw=lw, label='AUC = %0.3f' % auc7)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 7 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n7 >> query/curve in green, with AUC = {} and best threshold = {}'.format(auc7, best_thresh7),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("A7D7MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") 
                os.remove("dfmA1A2A3A4A5A6.csv") 
                os.remove("dfmA1A2A3A4A5A6A7.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4A5A6 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 6 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv")
                os.remove("dfmA1A2A3A4A5A6.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4A5 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 5 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 4 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3 is True:
                rocdata = pd.read_csv("dfmA1A2A3.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running! :-) ',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 3 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2 is True:
                rocdata = pd.read_csv("dfmA1A2.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 2 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv") # add a # at the beginning if u wanna keep these files
            else:
                if existA1A2 is False:
                    rocdata = pd.read_csv("A1D1MergedActInactTCOutcome.csv")
                    plt.figure()
                    lw = 2
                    fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                    auc1 = auc(fpr1,tpr1)
                    J1 = tpr1 - fpr1
                    ix1 = argmax(J1)
                    best_thresh1 = thresholds1[ix1]
                    plt.plot(fpr1, tpr1, color='darkorange',
                 lw=lw, label='AUC = %0.3f' % auc1)
                    parentdirectory = os.getcwd()
                    for fileQ1 in os.listdir(parentdirectory):
                        if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                            nmfileQ1 = fileQ1.split('_')[2]
                            #print(nmfileQ1)
                            # saving results and writting it in txt file
                            archv = open('out_rocsvalidroc.txt', 'w+')
                            archv.writelines(['\nRocsValidRoc has finished running!',
                                              '\n',
                                              '\nHere are your obtained results:',
                                              '\n',
                                              '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 1 query for which corresponding results are summarized next:',
                                              '\n',
                                              '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                              '\n',
                                              '\nThank you for using RocsValidRoc :-)',
                                              '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                            archv.close()
                    plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                    plt.xlim([-0.02, 1.02])
                    plt.ylim([-0.02, 1.02])
                    plt.xlabel('False Positive Rate')
                    plt.ylabel('True Positive Rate')
                    plt.title('ROC Curve')
                    plt.legend(loc="lower right")
                    plt.show()
                    os.remove("A1D1MergedActInactTCOutcome.csv")
                    #print(str('Oops.. something went wrong! Try again!'))

            print(str('Done! Thank you for using RocsValidRoc :-) '))
            print('Author: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil\n')
            # uuuufffffffff
            # THE END of option4
                    
############################ Answer for option 5 ##############################
# Only ROCS mode
    elif resposta == 5:    
        cabeçalho('Awesome! Let´s overlap some shapes...')        
        parent_dir = os.getcwd()
        print(str('We are working at'), parent_dir, str('\n'))
        nome_projeto = str(input('\nSet the name of your ROCS project, that u want to be printed in output folder\n>>'))
                
# DATABASES PART
        print(linha()) 
        print('Now check the available DATABASES for running ROCS: \n')

# list databases containing mol2 or sdf or oeb in filenames and present the list to user
# opens the list with databases names and address
        pattern1 = '*active*'
        pattern2 = '*inactive*'
        pattern3 = '*decoy*'
        parent_dir = os.getcwd()
        c = 1
        name1_base = []
        end1_base = []
        for root, dirs, files in os.walk(os.path.join(parent_dir)):
            for name1 in files:
                if fnmatch(name1, pattern1) or fnmatch(name1, pattern2) or fnmatch(name1, pattern3):
                    print(f'{c} - {name1}')
                    c += 1
                    sleep(0.01)

# gather results in name lists and addresses
                    name1_base.append(name1)
                    end1_base.append(os.path.join(root,name1))
                    #print(name1_base) # optional/if necessary

# select database within list
        print('\nSelect which DATABASES u wanna use (separate them by comma, eg, 1,2)')
        #base_proj = [] # optional/if necessary
        end_proj1 = []
        w = input()
        #print(w.split(',')) # optional/if necessary
        print(linha()) 
        print('You have selected the following DATABASES: \n')
        list = w.split(',')
        for i in list:
            x = int(i) - 1
            print(name1_base[x])
            #print(end1_base[x]) # optional/if necessary

# gather selected databases (name_base) in one list (base_proj)
            end_proj1.append(end1_base[x])
            #print(end_proj1) # optional/if necessary
                                  
# ask user to check databases
        print(linha())
        k = str(input('Is that correct? Shall we proceed? [y/n] '))
        print()
        if k == 'n':
            cabeçalho('Damn! Let´s try again! ')
            sleep(3)
            continue

        elif k == 'y':
                           
# QUERIES PART
            print(linha()) 
            print('Now let´s go to the QUERIES part. See what´s available here: \n')
            parent_dir = os.getcwd()
            name2 = WalkDirs_ListPattern('*.mol2*', '*.sdf*', '*.sq*', os.path.join(parent_dir) )
            for lst in name2:
                print(f'{lst}')
                sleep(0.01)
    
            y = os.getcwd()
            z = str(input('\nType the exact full filename of QUERY (or QUERIES if in one single molecule file) with extension (Eg: malato.mol2)\n** if your query file is not here, type findquery **\n>>'))
            print()
            if z == 'findquery':
                y = str(input('Please confirm the path of directory containing QUERY (Eg: C:\Gui\RocsValidRoc\queries)\n>>'))
                zz = str(input('\nInsert the exact full filename of QUERY with extension (Eg: malato.mol2)\n>>'))
                z = zz
            else:
                print()
                    
# stablishing other ROCS flags and running it        
            print(linha()) 
            print('Some other ROCS settings:')
            
            h = str(input('\nChoose extension/format for output results ´oformat´ (Eg: sdf, sdf.gz, oeb.gz, oeb)\n** do not use . dot before/after the extension name **\n>>'))
            hitsbest = str(input('\nNow choose how many hits (besthits) do you want to obtain in each ROCS run\n** if you are just running ROCS in advance to ROC curve validation use 0 **\n>>'))
            mpinp = str(input('\nLast, but not least, specify the number of processors to run ROCS in MPI mode\n** if this doesn´t matter for you, just type 0 **\n>>'))
           
            print(linha())
            print('Alrite!\nIf OpenEye logo comes up it is everything ok...'
                      '\nLet´s do this!')
            sleep(5)

# Criar diretório para os resultados do ROCS na pasta do usuário
            directory = '[OutputRocs]' + str(nome_projeto)  # Pode colocar a data
            parent_dir = os.getcwd()
            path_out_rocs = os.path.join(parent_dir, directory)
            [os.makedirs(i, exist_ok=True) for i in [path_out_rocs]]

# starting at date and time
            now = datetime.now()
            day_inicio_proces = now.strftime('%d.%m.%Y')
            hour_inicio_proces = now.strftime('%H:%M:%S')
            
# databases/query loop and further settings
            for base in end_proj1:
                
                file_name1 = os.path.basename(os.path.join(base))
                index_of_dot1 = file_name1.index('.')
                i = file_name1[:index_of_dot1]
                dbname = i.split(']')[1]
                
                queryname = z.split('.')[0]
                
# show date and time in terminal for each ROCS run
                now = datetime.now()
                day_inicio = now.strftime('%d.%m.%Y')
                hour_inicio = now.strftime('%H:%M:%S')
                print(linha())
                print()
                print('Running ROCS for database', i)
                print('on {} at {}'.format(day_inicio, hour_inicio))
                print()
                print(linha())

# create txt file to write date and time of each ROCS run
                arquivo = open('time_rocs.txt', 'a')
                arquivo.writelines(['\nRunning ROCS for database: ', i,
                                    '\n on {} at {}'.format(day_inicio, hour_inicio),
                                    '\n',
                                    ])
                arquivo.close()
                                   
# ROCS run
                subprocess.run(['rocs.bat',
                                '-dbase', os.path.join(base),
                                '-query', os.path.join(y, z),
                                '-outputdir', os.path.join(path_out_rocs),
                                '-prefix', 'OutRocs_' + dbname + '_' + queryname,
                                '-oformat', h,
                                '-besthits', hitsbest,
                                '-mpi_np', mpinp
                                ])

# end of loop Rocs
# Final message
            now = datetime.now()
            day_final = now.strftime('%d.%m.%Y')
            hour_final = now.strftime('%H:%M:%S')
            print(linha())
            print('ROCS has finished running!')
            print('Started {} at {}'.format(day_inicio_proces, hour_inicio_proces))
            print('Ended {} at {}'.format(day_final, hour_final))
            print('Thank you for using AutomaROCS :-)')
            print('Author: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil\n')

# writting final message in txt file
            arquivo = open('time_rocs.txt', 'a')
            arquivo.writelines(['\nROCS has finished running!',
                                '\nStarted {} at {}'.format(day_inicio_proces, hour_inicio_proces),
                                '\nEnded {} at {}'.format(day_final, hour_final),
                                '\nThank you for using AutomaROCS :-)',
                                '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
            arquivo.close()

            
############################ Answer for option 6 ##############################
    elif resposta == 6:
            cabeçalho('Yey let´s do some ROC curve!')
            
            ra = str(input('Please inform us the path of directory containing rpt files, ie, folder with output results of ROCS (Eg: C:\Gui\RocsValidRoc\[OutputRocs]test)\n>>'))  
            #print(ra)
            os.chdir(os.path.join(ra))
            print("Current working directory: {0}".format(os.getcwd()))
            
# RETRIEVING RESULTS OF QUERY (1) FOR ACTIVE_1 (A1) and DECOYS_1 (D1)
            parentdirectory = os.getcwd()
            for fileA1 in os.listdir(parentdirectory):
                if fnmatch(fileA1, '*active*') and fileA1.endswith("_1.rpt"):
                        with open(fileA1) as fin, open('fileA1.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA1 = pd.read_csv("fileA1.csv") 
                        A1data = rawA1.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A1data.head()
                        A1data.insert(2, "Outcome", "1")
                        A1data.head()
                        A1data.to_csv("A1file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD1 in os.listdir(parentdirectory):
                if fnmatch(fileD1, '*decoy*') and fileD1.endswith("_1.rpt"):
                        with open(fileD1) as fin, open('fileD1.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD1 = pd.read_csv("fileD1.csv") 
                        D1data = rawD1.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D1data.head()
                        D1data.insert(2, "Outcome", "0")
                        D1data.head()
                        D1data.to_csv("D1file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA1 = str('fileA1.csv')
            existencialA1 = os.path.exists(os.path.join(parentdirectory, existA1))
            if existencialA1 is True:
                file1A1 = open("A1file.csv", "a")
                file2D1 = open("D1file.csv", "r")
                for line in file2D1:
                    file1A1.write(line)
                file1A1.close()
                file2D1.close()
                sortdataA1 = pd.read_csv("A1file.csv") 
                sortdataaA1 = sortdataA1.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA1 = sortdataaA1.drop(sortdataaA1.index[0])
                sortdataaaaA1 = sortdataaaA1.rename(columns={"TanimotoCombo": "TanimotoCombo1", "Outcome": "Outcome1"}) #attention
                sortdataaaaaA1 = sortdataaaaA1.reset_index()
                sortdataaaaaA1.insert(0, 'New_ID', 0)
                sortdataaaaaA1['New_ID'] = sortdataaaaaA1.index + 0
                sortdataaaaaA1.to_csv("A1D1MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A1file.csv")
                os.remove("fileA1.csv")
                os.remove("fileD1.csv")
                os.remove("D1file.csv")
            else:
                print(str('The existence of a query A1/D1 is'), existencialA1)

                
# RETRIEVING RESULTS OF QUERY (2) FOR ACTIVE_2 (A2) and DECOYS_2 (D2)
            parentdirectory = os.getcwd()
            for fileA2 in os.listdir(parentdirectory):
                if fnmatch(fileA2, '*active*') and fileA2.endswith("_2.rpt"):
                        with open(fileA2) as fin, open('fileA2.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA2 = pd.read_csv("fileA2.csv") 
                        A2data = rawA2.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A2data.head()
                        A2data.insert(2, "Outcome", "1")
                        A2data.head()
                        A2data.to_csv("A2file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD2 in os.listdir(parentdirectory):
                if fnmatch(fileD2, '*decoy*') and fileD2.endswith("_2.rpt"):
                        with open(fileD2) as fin, open('fileD2.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD2 = pd.read_csv("fileD2.csv") 
                        D2data = rawD2.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D2data.head()
                        D2data.insert(2, "Outcome", "0")
                        D2data.head()
                        D2data.to_csv("D2file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA2 = str('fileA2.csv')
            existencialA2 = os.path.exists(os.path.join(parentdirectory, existA2))
            if existencialA2 is True:
                file1A2 = open("A2file.csv", "a")
                file2D2 = open("D2file.csv", "r")
                for line in file2D2:
                    file1A2.write(line)
                file1A2.close()
                file2D2.close()
                sortdataA2 = pd.read_csv("A2file.csv") 
                sortdataaA2 = sortdataA2.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA2 = sortdataaA2.drop(sortdataaA2.index[0])
                sortdataaaaA2 = sortdataaaA2.rename(columns={"TanimotoCombo": "TanimotoCombo2", "Outcome": "Outcome2"}) #attention
                sortdataaaaaA2 = sortdataaaaA2.reset_index()
                sortdataaaaaA2.insert(0, 'New_ID', 0)
                sortdataaaaaA2['New_ID'] = sortdataaaaaA2.index + 0
                sortdataaaaaA2.to_csv("A2D2MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A2file.csv")
                os.remove("fileA2.csv")
                os.remove("fileD2.csv")
                os.remove("D2file.csv")
            else:
                print(str('The existence of a second query A2/D2 is'), existencialA2)

# RETRIEVING RESULTS OF QUERY (3) FOR ACTIVE_3 (A3) and DECOYS_3 (D3)
            parentdirectory = os.getcwd()
            for fileA3 in os.listdir(parentdirectory):
                if fnmatch(fileA3, '*active*') and fileA3.endswith("_3.rpt"):
                        with open(fileA3) as fin, open('fileA3.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA3 = pd.read_csv("fileA3.csv") 
                        A3data = rawA3.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A3data.head()
                        A3data.insert(2, "Outcome", "1")
                        A3data.head()
                        A3data.to_csv("A3file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD3 in os.listdir(parentdirectory):
                if fnmatch(fileD3, '*decoy*') and fileD3.endswith("_3.rpt"):
                        with open(fileD3) as fin, open('fileD3.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD3 = pd.read_csv("fileD3.csv") 
                        D3data = rawD3.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D3data.head()
                        D3data.insert(2, "Outcome", "0")
                        D3data.head()
                        D3data.to_csv("D3file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA3 = str('fileA3.csv')
            existencialA3 = os.path.exists(os.path.join(parentdirectory, existA3))
            if existencialA3 is True:
                file1A3 = open("A3file.csv", "a")
                file2D3 = open("D3file.csv", "r")
                for line in file2D3:
                    file1A3.write(line)
                file1A3.close()
                file2D3.close()
                sortdataA3 = pd.read_csv("A3file.csv") 
                sortdataaA3 = sortdataA3.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA3 = sortdataaA3.drop(sortdataaA3.index[0])
                sortdataaaaA3 = sortdataaaA3.rename(columns={"TanimotoCombo": "TanimotoCombo3", "Outcome": "Outcome3"}) #attention
                sortdataaaaaA3 = sortdataaaaA3.reset_index()
                sortdataaaaaA3.insert(0, 'New_ID', 0)
                sortdataaaaaA3['New_ID'] = sortdataaaaaA3.index + 0
                sortdataaaaaA3.to_csv("A3D3MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A3file.csv")
                os.remove("fileA3.csv")
                os.remove("fileD3.csv")
                os.remove("D3file.csv")
            else:
                print(str('The existence of a third query A3/D3 is'), existencialA3)

# RETRIEVING RESULTS OF QUERY (4) FOR ACTIVE_4 (A4) and DECOYS_4 (D4)
            parentdirectory = os.getcwd()
            for fileA4 in os.listdir(parentdirectory):
                if fnmatch(fileA4, '*active*') and fileA4.endswith("_4.rpt"):
                        with open(fileA4) as fin, open('fileA4.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA4 = pd.read_csv("fileA4.csv") 
                        A4data = rawA4.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A4data.head()
                        A4data.insert(2, "Outcome", "1")
                        A4data.head()
                        A4data.to_csv("A4file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD4 in os.listdir(parentdirectory):
                if fnmatch(fileD4, '*decoy*') and fileD4.endswith("_4.rpt"):
                        with open(fileD4) as fin, open('fileD4.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD4 = pd.read_csv("fileD4.csv") 
                        D4data = rawD4.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D4data.head()
                        D4data.insert(2, "Outcome", "0")
                        D4data.head()
                        D4data.to_csv("D4file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA4 = str('fileA4.csv')
            existencialA4 = os.path.exists(os.path.join(parentdirectory, existA4))
            if existencialA4 is True:
                file1A4 = open("A4file.csv", "a")
                file2D4 = open("D4file.csv", "r")
                for line in file2D4:
                    file1A4.write(line)
                file1A4.close()
                file2D4.close()
                sortdataA4 = pd.read_csv("A4file.csv") 
                sortdataaA4 = sortdataA4.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA4 = sortdataaA4.drop(sortdataaA4.index[0])
                sortdataaaaA4 = sortdataaaA4.rename(columns={"TanimotoCombo": "TanimotoCombo4", "Outcome": "Outcome4"}) #attention
                sortdataaaaaA4 = sortdataaaaA4.reset_index()
                sortdataaaaaA4.insert(0, 'New_ID', 0)
                sortdataaaaaA4['New_ID'] = sortdataaaaaA4.index + 0
                sortdataaaaaA4.to_csv("A4D4MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A4file.csv")
                os.remove("fileA4.csv")
                os.remove("fileD4.csv")
                os.remove("D4file.csv")
            else:
                print(str('The existence of a fourth query A4/D4 is'), existencialA4)

# RETRIEVING RESULTS OF QUERY (5) FOR ACTIVE_5 (A5) and DECOYS_5 (D5)
            parentdirectory = os.getcwd()
            for fileA5 in os.listdir(parentdirectory):
                if fnmatch(fileA5, '*active*') and fileA5.endswith("_5.rpt"):
                        with open(fileA5) as fin, open('fileA5.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA5 = pd.read_csv("fileA5.csv") 
                        A5data = rawA5.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A5data.head()
                        A5data.insert(2, "Outcome", "1")
                        A5data.head()
                        A5data.to_csv("A5file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD5 in os.listdir(parentdirectory):
                if fnmatch(fileD5, '*decoy*') and fileD5.endswith("_5.rpt"):
                        with open(fileD5) as fin, open('fileD5.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD5 = pd.read_csv("fileD5.csv") 
                        D5data = rawD5.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D5data.head()
                        D5data.insert(2, "Outcome", "0")
                        D5data.head()
                        D5data.to_csv("D5file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA5 = str('fileA5.csv')
            existencialA5 = os.path.exists(os.path.join(parentdirectory, existA5))
            if existencialA5 is True:
                file1A5 = open("A5file.csv", "a")
                file2D5 = open("D5file.csv", "r")
                for line in file2D5:
                    file1A5.write(line)
                file1A5.close()
                file2D5.close()
                sortdataA5 = pd.read_csv("A5file.csv") 
                sortdataaA5 = sortdataA5.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA5 = sortdataaA5.drop(sortdataaA5.index[0])
                sortdataaaaA5 = sortdataaaA5.rename(columns={"TanimotoCombo": "TanimotoCombo5", "Outcome": "Outcome5"}) #attention
                sortdataaaaaA5 = sortdataaaaA5.reset_index()
                sortdataaaaaA5.insert(0, 'New_ID', 0)
                sortdataaaaaA5['New_ID'] = sortdataaaaaA5.index + 0
                sortdataaaaaA5.to_csv("A5D5MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A5file.csv")
                os.remove("fileA5.csv")
                os.remove("fileD5.csv")
                os.remove("D5file.csv")
            else:
                print(str('The existence of a fifth query A5/D5 is'), existencialA5)

# RETRIEVING RESULTS OF QUERY (6) FOR ACTIVE_6 (A6) and DECOYS_6 (D6)
            parentdirectory = os.getcwd()
            for fileA6 in os.listdir(parentdirectory):
                if fnmatch(fileA6, '*active*') and fileA6.endswith("_6.rpt"):
                        with open(fileA6) as fin, open('fileA6.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA6 = pd.read_csv("fileA6.csv") 
                        A6data = rawA6.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A6data.head()
                        A6data.insert(2, "Outcome", "1")
                        A6data.head()
                        A6data.to_csv("A6file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD6 in os.listdir(parentdirectory):
                if fnmatch(fileD6, '*decoy*') and fileD6.endswith("_6.rpt"):
                        with open(fileD6) as fin, open('fileD6.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD6 = pd.read_csv("fileD6.csv") 
                        D6data = rawD6.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D6data.head()
                        D6data.insert(2, "Outcome", "0")
                        D6data.head()
                        D6data.to_csv("D6file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA6 = str('fileA6.csv')
            existencialA6 = os.path.exists(os.path.join(parentdirectory, existA6))
            if existencialA6 is True:
                file1A6 = open("A6file.csv", "a")
                file2D6 = open("D6file.csv", "r")
                for line in file2D6:
                    file1A6.write(line)
                file1A6.close()
                file2D6.close()
                sortdataA6 = pd.read_csv("A6file.csv") 
                sortdataaA6 = sortdataA6.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA6 = sortdataaA6.drop(sortdataaA6.index[0])
                sortdataaaaA6 = sortdataaaA6.rename(columns={"TanimotoCombo": "TanimotoCombo6", "Outcome": "Outcome6"}) #attention
                sortdataaaaaA6 = sortdataaaaA6.reset_index()
                sortdataaaaaA6.insert(0, 'New_ID', 0)
                sortdataaaaaA6['New_ID'] = sortdataaaaaA6.index + 0
                sortdataaaaaA6.to_csv("A6D6MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A6file.csv")
                os.remove("fileA6.csv")
                os.remove("fileD6.csv")
                os.remove("D6file.csv")
            else:
                print(str('The existence of a sixth query A6/D6 is'), existencialA6)

# RETRIEVING RESULTS OF QUERY (7) FOR ACTIVE_7 (A7) and DECOYS_7 (D7)
            parentdirectory = os.getcwd()
            for fileA7 in os.listdir(parentdirectory):
                if fnmatch(fileA7, '*active*') and fileA7.endswith("_7.rpt"):
                        with open(fileA7) as fin, open('fileA7.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA7 = pd.read_csv("fileA7.csv") 
                        A7data = rawA7.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A7data.head()
                        A7data.insert(2, "Outcome", "1")
                        A7data.head()
                        A7data.to_csv("A7file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD7 in os.listdir(parentdirectory):
                if fnmatch(fileD7, '*decoy*') and fileD7.endswith("_7.rpt"):
                        with open(fileD7) as fin, open('fileD7.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD7 = pd.read_csv("fileD7.csv") 
                        D7data = rawD7.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D7data.head()
                        D7data.insert(2, "Outcome", "0")
                        D7data.head()
                        D7data.to_csv("D7file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA7 = str('fileA7.csv')
            existencialA7 = os.path.exists(os.path.join(parentdirectory, existA7))
            if existencialA7 is True:
                file1A7 = open("A7file.csv", "a")
                file2D7 = open("D7file.csv", "r")
                for line in file2D7:
                    file1A7.write(line)
                file1A7.close()
                file2D7.close()
                sortdataA7 = pd.read_csv("A7file.csv") 
                sortdataaA7 = sortdataA7.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA7 = sortdataaA7.drop(sortdataaA7.index[0])
                sortdataaaaA7 = sortdataaaA7.rename(columns={"TanimotoCombo": "TanimotoCombo7", "Outcome": "Outcome7"}) #attention
                sortdataaaaaA7 = sortdataaaaA7.reset_index()
                sortdataaaaaA7.insert(0, 'New_ID', 0)
                sortdataaaaaA7['New_ID'] = sortdataaaaaA7.index + 0
                sortdataaaaaA7.to_csv("A7D7MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A7file.csv")
                os.remove("fileA7.csv")
                os.remove("fileD7.csv")
                os.remove("D7file.csv")
            else:
                print(str('The existence of a seventh query A7/D7 is'), existencialA7)

# RETRIEVING RESULTS OF QUERY (8) FOR ACTIVE_8 (A8) and DECOYS_8 (D8)
            parentdirectory = os.getcwd()
            for fileA8 in os.listdir(parentdirectory):
                if fnmatch(fileA8, '*active*') and fileA8.endswith("_8.rpt"):
                        with open(fileA8) as fin, open('fileA8.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA8 = pd.read_csv("fileA8.csv") 
                        A8data = rawA8.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A8data.head()
                        A8data.insert(2, "Outcome", "1")
                        A8data.head()
                        A8data.to_csv("A8file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD8 in os.listdir(parentdirectory):
                if fnmatch(fileD8, '*decoy*') and fileD8.endswith("_8.rpt"):
                        with open(fileD8) as fin, open('fileD8.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD8 = pd.read_csv("fileD8.csv") 
                        D8data = rawD8.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D8data.head()
                        D8data.insert(2, "Outcome", "0")
                        D8data.head()
                        D8data.to_csv("D8file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA8 = str('fileA8.csv')
            existencialA8 = os.path.exists(os.path.join(parentdirectory, existA8))
            if existencialA8 is True:
                file1A8 = open("A8file.csv", "a")
                file2D8 = open("D8file.csv", "r")
                for line in file2D8:
                    file1A8.write(line)
                file1A8.close()
                file2D8.close()
                sortdataA8 = pd.read_csv("A8file.csv") 
                sortdataaA8 = sortdataA8.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA8 = sortdataaA8.drop(sortdataaA8.index[0])
                sortdataaaaA8 = sortdataaaA8.rename(columns={"TanimotoCombo": "TanimotoCombo8", "Outcome": "Outcome8"}) #attention
                sortdataaaaaA8 = sortdataaaaA8.reset_index()
                sortdataaaaaA8.insert(0, 'New_ID', 0)
                sortdataaaaaA8['New_ID'] = sortdataaaaaA8.index + 0
                sortdataaaaaA8.to_csv("A8D8MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A8file.csv")
                os.remove("fileA8.csv")
                os.remove("fileD8.csv")
                os.remove("D8file.csv")
            else:
                print(str('The existence of a eigth query A8/D8 is'), existencialA8)

# RETRIEVING RESULTS OF QUERY (9) FOR ACTIVE_9 (A9) and DECOYS_9 (D9)
            parentdirectory = os.getcwd()
            for fileA9 in os.listdir(parentdirectory):
                if fnmatch(fileA9, '*active*') and fileA9.endswith("_9.rpt"):
                        with open(fileA9) as fin, open('fileA9.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA9 = pd.read_csv("fileA9.csv") 
                        A9data = rawA9.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A9data.head()
                        A9data.insert(2, "Outcome", "1")
                        A9data.head()
                        A9data.to_csv("A9file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD9 in os.listdir(parentdirectory):
                if fnmatch(fileD9, '*decoy*') and fileD9.endswith("_9.rpt"):
                        with open(fileD9) as fin, open('fileD9.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD9 = pd.read_csv("fileD9.csv") 
                        D9data = rawD9.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D9data.head()
                        D9data.insert(2, "Outcome", "0")
                        D9data.head()
                        D9data.to_csv("D9file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA9 = str('fileA9.csv')
            existencialA9 = os.path.exists(os.path.join(parentdirectory, existA9))
            if existencialA9 is True:
                file1A9 = open("A9file.csv", "a")
                file2D9 = open("D9file.csv", "r")
                for line in file2D9:
                    file1A9.write(line)
                file1A9.close()
                file2D9.close()
                sortdataA9 = pd.read_csv("A9file.csv") 
                sortdataaA9 = sortdataA9.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA9 = sortdataaA9.drop(sortdataaA9.index[0])
                sortdataaaaA9 = sortdataaaA9.rename(columns={"TanimotoCombo": "TanimotoCombo9", "Outcome": "Outcome9"}) #attention
                sortdataaaaaA9 = sortdataaaaA9.reset_index()
                sortdataaaaaA9.insert(0, 'New_ID', 0)
                sortdataaaaaA9['New_ID'] = sortdataaaaaA9.index + 0
                sortdataaaaaA9.to_csv("A9D9MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A9file.csv")
                os.remove("fileA9.csv")
                os.remove("fileD9.csv")
                os.remove("D9file.csv")
            else:
                print(str('The existence of a ninth query A9/D9 is'), existencialA9)

# RETRIEVING RESULTS OF QUERY (10) FOR ACTIVE_10 (A10) and DECOYS_10 (D10)
            parentdirectory = os.getcwd()
            for fileA10 in os.listdir(parentdirectory):
                if fnmatch(fileA10, '*active*') and fileA10.endswith("_10.rpt"):
                        with open(fileA10) as fin, open('fileA10.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawA10 = pd.read_csv("fileA10.csv") 
                        A10data = rawA10.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        A10data.head()
                        A10data.insert(2, "Outcome", "1")
                        A10data.head()
                        A10data.to_csv("A10file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue
        
            parentdirectory = os.getcwd()
            for fileD10 in os.listdir(parentdirectory):
                if fnmatch(fileD10, '*decoy*') and fileD10.endswith("_10.rpt"):
                        with open(fileD10) as fin, open('fileD10.csv', 'w+') as fout:
                            for line in fin:
                                fout.write(line.replace('\t', ','))
                        rawD10 = pd.read_csv("fileD10.csv") 
                        D10data = rawD10.drop(['ShapeQuery', 'Rank', 'ShapeTanimoto', 'ColorTanimoto', 'FitTverskyCombo', 'FitTversky', 'FitColorTversky', 'RefTverskyCombo', 'RefTversky', 'RefColorTversky', 'ColorScore', 'Overlap'], axis=1)
                        D10data.head()
                        D10data.insert(2, "Outcome", "0")
                        D10data.head()
                        D10data.to_csv("D10file.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                else:
        
                        continue

            parentdirectory = os.getcwd()
            existA10 = str('fileA10.csv')
            existencialA10 = os.path.exists(os.path.join(parentdirectory, existA10))
            if existencialA10 is True:
                file1A10 = open("A10file.csv", "a")
                file2D10 = open("D10file.csv", "r")
                for line in file2D10:
                    file1A10.write(line)
                file1A10.close()
                file2D10.close()
                sortdataA10 = pd.read_csv("A10file.csv") 
                sortdataaA10 = sortdataA10.sort_values(by=['TanimotoCombo'], ascending=False)
                sortdataaaA10 = sortdataaA10.drop(sortdataaA10.index[0])
                sortdataaaaA10 = sortdataaaA10.rename(columns={"TanimotoCombo": "TanimotoCombo10", "Outcome": "Outcome10"}) #attention
                sortdataaaaaA10 = sortdataaaaA10.reset_index()
                sortdataaaaaA10.insert(0, 'New_ID', 0)
                sortdataaaaaA10['New_ID'] = sortdataaaaaA10.index + 0
                sortdataaaaaA10.to_csv("A10D10MergedActInactTCOutcome.csv", sep=',', encoding='utf-8', index=False, mode='w+')
                os.remove("A10file.csv")
                os.remove("fileA10.csv")
                os.remove("fileD10.csv")
                os.remove("D10file.csv")
            else:
                print(str('The existence of a tenth query A10/D10 is'), existencialA10)

# merging tables of different queries
            parentdirectory = os.getcwd()
            gerouA1 = str('A1D1MergedActInactTCOutcome.csv')
            gerouA2 = str('A2D2MergedActInactTCOutcome.csv')
            gerouA3 = str('A3D3MergedActInactTCOutcome.csv')
            gerouA4 = str('A4D4MergedActInactTCOutcome.csv')
            gerouA5 = str('A5D5MergedActInactTCOutcome.csv')
            gerouA6 = str('A6D6MergedActInactTCOutcome.csv')
            gerouA7 = str('A7D7MergedActInactTCOutcome.csv')
            gerouA8 = str('A8D8MergedActInactTCOutcome.csv')
            gerouA9 = str('A9D9MergedActInactTCOutcome.csv')
            gerouA10 = str('A10D10MergedActInactTCOutcome.csv')
            existegerA1 = os.path.exists(os.path.join(parentdirectory, gerouA1))
            existegerA2 = os.path.exists(os.path.join(parentdirectory, gerouA2))
            existegerA3 = os.path.exists(os.path.join(parentdirectory, gerouA3))
            existegerA4 = os.path.exists(os.path.join(parentdirectory, gerouA4))
            existegerA5 = os.path.exists(os.path.join(parentdirectory, gerouA5))
            existegerA6 = os.path.exists(os.path.join(parentdirectory, gerouA6))
            existegerA7 = os.path.exists(os.path.join(parentdirectory, gerouA7))
            existegerA8 = os.path.exists(os.path.join(parentdirectory, gerouA8))
            existegerA9 = os.path.exists(os.path.join(parentdirectory, gerouA9))
            existegerA10 = os.path.exists(os.path.join(parentdirectory, gerouA10))
            if existegerA1 and existegerA2 is True:
                dfmA1 = pd.read_csv("A1D1MergedActInactTCOutcome.csv")
                dfmA2 = pd.read_csv("A2D2MergedActInactTCOutcome.csv")
                dfmA1A2 = (dfmA2.merge(dfmA1, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                #print(dfmA1A2)
                dfmA1A2data = dfmA1A2.drop(['index', 'index_'], axis=1)
                dfmA1A2data.to_csv("dfmA1A2.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('...'))

            gerouA1A2 = str('dfmA1A2.csv')
            existA1A2 = os.path.exists(os.path.join(parentdirectory, gerouA1A2))
            if existA1A2 and existegerA3 is True:
                df2mA1A2 = pd.read_csv("dfmA1A2.csv")
                dfmA3 = pd.read_csv("A3D3MergedActInactTCOutcome.csv")
                dfmA1A2A3 = (dfmA3.merge(df2mA1A2, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3.to_csv("dfmA1A2A3.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('....'))
                    
            gerouA1A2A3 = str('dfmA1A2A3.csv')
            existA1A2A3 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3))
            if existA1A2A3 and existegerA4 is True:
                df2mA1A2A3 = pd.read_csv("dfmA1A2A3.csv")
                dfmA4 = pd.read_csv("A4D4MergedActInactTCOutcome.csv")
                dfmA1A2A3A4 = (dfmA4.merge(df2mA1A2A3, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4.to_csv("dfmA1A2A3A4.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('.....'))

            gerouA1A2A3A4 = str('dfmA1A2A3A4.csv')
            existA1A2A3A4 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4))
            if existA1A2A3A4 and existegerA5 is True:
                df2mA1A2A3A4 = pd.read_csv("dfmA1A2A3A4.csv")
                dfmA5 = pd.read_csv("A5D5MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5 = (dfmA5.merge(df2mA1A2A3A4, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5.to_csv("dfmA1A2A3A4A5.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('......'))   
                
            gerouA1A2A3A4A5 = str('dfmA1A2A3A4A5.csv')
            existA1A2A3A4A5 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5))
            if existA1A2A3A4A5 and existegerA6 is True:
                df2mA1A2A3A4A5 = pd.read_csv("dfmA1A2A3A4A5.csv")
                dfmA6 = pd.read_csv("A6D6MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6 = (dfmA6.merge(df2mA1A2A3A4A5, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6.to_csv("dfmA1A2A3A4A5A6.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('.......'))
                
            gerouA1A2A3A4A5A6 = str('dfmA1A2A3A4A5A6.csv')
            existA1A2A3A4A5A6 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6))
            if existA1A2A3A4A5A6 and existegerA7 is True:
                df2mA1A2A3A4A5A6 = pd.read_csv("dfmA1A2A3A4A5A6.csv")
                dfmA7 = pd.read_csv("A7D7MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6A7 = (dfmA7.merge(df2mA1A2A3A4A5A6, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6A7.to_csv("dfmA1A2A3A4A5A6A7.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('........'))
                
            gerouA1A2A3A4A5A6A7 = str('dfmA1A2A3A4A5A6A7.csv')
            existA1A2A3A4A5A6A7 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6A7))
            if existA1A2A3A4A5A6A7 and existegerA8 is True:
                df2mA1A2A3A4A5A6A7 = pd.read_csv("dfmA1A2A3A4A5A6A7.csv")
                dfmA8 = pd.read_csv("A8D8MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6A7A8 = (dfmA8.merge(df2mA1A2A3A4A5A6A7, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6A7A8.to_csv("dfmA1A2A3A4A5A6A7A8.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('.........'))
                
            gerouA1A2A3A4A5A6A7A8 = str('dfmA1A2A3A4A5A6A7A8.csv')
            existA1A2A3A4A5A6A7A8 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6A7A8))
            if existA1A2A3A4A5A6A7A8 and existegerA9 is True:
                df2mA1A2A3A4A5A6A7A8 = pd.read_csv("dfmA1A2A3A4A5A6A7A8.csv")
                dfmA9 = pd.read_csv("A9D9MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6A7A8A9 = (dfmA9.merge(df2mA1A2A3A4A5A6A7A8, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6A7A8A9.to_csv("dfmA1A2A3A4A5A6A7A8A9.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('..........'))
                
            gerouA1A2A3A4A5A6A7A8A9 = str('dfmA1A2A3A4A5A6A7A8A9.csv')
            existA1A2A3A4A5A6A7A8A9 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6A7A8A9))
            if existA1A2A3A4A5A6A7A8A9 and existegerA10 is True:
                df2mA1A2A3A4A5A6A7A8A9 = pd.read_csv("dfmA1A2A3A4A5A6A7A8A9.csv")
                dfmA10 = pd.read_csv("A10D10MergedActInactTCOutcome.csv")
                dfmA1A2A3A4A5A6A7A8A9A10 = (dfmA10.merge(df2mA1A2A3A4A5A6A7A8A9, left_on='New_ID', right_on='New_ID', suffixes=('','_')))
                dfmA1A2A3A4A5A6A7A8A9A10.to_csv("dfmA1A2A3A4A5A6A7A8A9A10.csv", sep=',', encoding='utf-8', index=False, mode='w+')
            else:
                print(str('...........'))

# final adjustments and ploting ROC curve(S)
            gerouA1A2A3A4A5A6A7A8A9A10 = str('dfmA1A2A3A4A5A6A7A8A9A10.csv')
            existA1A2A3A4A5A6A7A8A9A10 = os.path.exists(os.path.join(parentdirectory, gerouA1A2A3A4A5A6A7A8A9A10))
            if existA1A2A3A4A5A6A7A8A9A10 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6A7A8A9A10.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                fpr7, tpr7, thresholds7 = roc_curve(rocdata['Outcome7'], rocdata['TanimotoCombo7'])
                fpr8, tpr8, thresholds8 = roc_curve(rocdata['Outcome8'], rocdata['TanimotoCombo8'])
                fpr9, tpr9, thresholds9 = roc_curve(rocdata['Outcome9'], rocdata['TanimotoCombo9'])
                fpr10, tpr10, thresholds10 = roc_curve(rocdata['Outcome10'], rocdata['TanimotoCombo10'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                auc7 = auc(fpr7,tpr7)
                auc8 = auc(fpr8,tpr8)
                auc9 = auc(fpr9,tpr9)
                auc10 = auc(fpr10,tpr10)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                J7 = tpr7 - fpr7
                ix7 = argmax(J7)
                best_thresh7 = thresholds7[ix7]
                J8 = tpr8 - fpr8
                ix8 = argmax(J8)
                best_thresh8 = thresholds8[ix8]
                J9 = tpr9 - fpr9
                ix9 = argmax(J9)
                best_thresh9 = thresholds9[ix9]
                J10 = tpr10 - fpr10
                ix10 = argmax(J10)
                best_thresh10 = thresholds10[ix10]                
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                plt.plot(fpr7, tpr7, color='green',
             lw=lw, label='AUC = %0.3f' % auc7)
                plt.plot(fpr8, tpr8, color='saddlebrown',
             lw=lw, label='AUC = %0.3f' % auc8)
                plt.plot(fpr9, tpr9, color='darkblue',
             lw=lw, label='AUC = %0.3f' % auc9)
                plt.plot(fpr10, tpr10, color='darkviolet',
             lw=lw, label='AUC = %0.3f' % auc10)
                # retrieving query mol2 info to report results
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 10 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n7 >> query/curve in green, with AUC = {} and best threshold = {}'.format(auc7, best_thresh7),
                                          '\n8 >> query/curve in brown, with AUC = {} and best threshold = {}'.format(auc8, best_thresh8),
                                          '\n9 >> query/curve in darkblue, with AUC = {} and best threshold = {}'.format(auc9, best_thresh9),
                                          '\n10 >> query/curve in purple, with AUC = {} and best threshold = {}'.format(auc10, best_thresh10),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                #cleaning
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("A7D7MergedActInactTCOutcome.csv")
                os.remove("A8D8MergedActInactTCOutcome.csv")
                os.remove("A9D9MergedActInactTCOutcome.csv")
                os.remove("A10D10MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") 
                os.remove("dfmA1A2A3A4A5A6.csv") 
                os.remove("dfmA1A2A3A4A5A6A7.csv") 
                os.remove("dfmA1A2A3A4A5A6A7A8.csv")
                os.remove("dfmA1A2A3A4A5A6A7A8A9.csv")
                os.remove("dfmA1A2A3A4A5A6A7A8A9A10.csv") # add a # at the beginning if u wanna keep these files
                
            elif existA1A2A3A4A5A6A7A8A9 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6A7A8A9.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                fpr7, tpr7, thresholds7 = roc_curve(rocdata['Outcome7'], rocdata['TanimotoCombo7'])
                fpr8, tpr8, thresholds8 = roc_curve(rocdata['Outcome8'], rocdata['TanimotoCombo8'])
                fpr9, tpr9, thresholds9 = roc_curve(rocdata['Outcome9'], rocdata['TanimotoCombo9'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                auc7 = auc(fpr7,tpr7)
                auc8 = auc(fpr8,tpr8)
                auc9 = auc(fpr9,tpr9)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                J7 = tpr7 - fpr7
                ix7 = argmax(J7)
                best_thresh7 = thresholds7[ix7]
                J8 = tpr8 - fpr8
                ix8 = argmax(J8)
                best_thresh8 = thresholds8[ix8]
                J9 = tpr9 - fpr9
                ix9 = argmax(J9)
                best_thresh9 = thresholds9[ix9]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                plt.plot(fpr7, tpr7, color='green',
             lw=lw, label='AUC = %0.3f' % auc7)
                plt.plot(fpr8, tpr8, color='saddlebrown',
             lw=lw, label='AUC = %0.3f' % auc8)
                plt.plot(fpr9, tpr9, color='darkblue',
             lw=lw, label='AUC = %0.3f' % auc9)
                
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 9 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n7 >> query/curve in green, with AUC = {} and best threshold = {}'.format(auc7, best_thresh7),
                                          '\n8 >> query/curve in brown, with AUC = {} and best threshold = {}'.format(auc8, best_thresh8),
                                          '\n9 >> query/curve in darkblue, with AUC = {} and best threshold = {}'.format(auc9, best_thresh9),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("A7D7MergedActInactTCOutcome.csv")
                os.remove("A8D8MergedActInactTCOutcome.csv")
                os.remove("A9D9MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") 
                os.remove("dfmA1A2A3A4A5A6.csv") 
                os.remove("dfmA1A2A3A4A5A6A7.csv") 
                os.remove("dfmA1A2A3A4A5A6A7A8.csv")
                os.remove("dfmA1A2A3A4A5A6A7A8A9.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4A5A6A7A8 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6A7A8.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                fpr7, tpr7, thresholds7 = roc_curve(rocdata['Outcome7'], rocdata['TanimotoCombo7'])
                fpr8, tpr8, thresholds8 = roc_curve(rocdata['Outcome8'], rocdata['TanimotoCombo8'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                auc7 = auc(fpr7,tpr7)
                auc8 = auc(fpr8,tpr8)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                J7 = tpr7 - fpr7
                ix7 = argmax(J7)
                best_thresh7 = thresholds7[ix7]
                J8 = tpr8 - fpr8
                ix8 = argmax(J8)
                best_thresh8 = thresholds8[ix8]
                
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                plt.plot(fpr7, tpr7, color='green',
             lw=lw, label='AUC = %0.3f' % auc7)
                plt.plot(fpr8, tpr8, color='saddlebrown',
             lw=lw, label='AUC = %0.3f' % auc8)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 8 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n7 >> query/curve in green, with AUC = {} and best threshold = {}'.format(auc7, best_thresh7),
                                          '\n8 >> query/curve in brown, with AUC = {} and best threshold = {}'.format(auc8, best_thresh8),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("A7D7MergedActInactTCOutcome.csv")
                os.remove("A8D8MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") 
                os.remove("dfmA1A2A3A4A5A6.csv") 
                os.remove("dfmA1A2A3A4A5A6A7.csv") 
                os.remove("dfmA1A2A3A4A5A6A7A8.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4A5A6A7 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6A7.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                fpr7, tpr7, thresholds7 = roc_curve(rocdata['Outcome7'], rocdata['TanimotoCombo7'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                auc7 = auc(fpr7,tpr7)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                J7 = tpr7 - fpr7
                ix7 = argmax(J7)
                best_thresh7 = thresholds7[ix7]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                plt.plot(fpr7, tpr7, color='green',
             lw=lw, label='AUC = %0.3f' % auc7)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 7 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n7 >> query/curve in green, with AUC = {} and best threshold = {}'.format(auc7, best_thresh7),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("A7D7MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") 
                os.remove("dfmA1A2A3A4A5A6.csv") 
                os.remove("dfmA1A2A3A4A5A6A7.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4A5A6 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5A6.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                fpr6, tpr6, thresholds6 = roc_curve(rocdata['Outcome6'], rocdata['TanimotoCombo6'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                auc6 = auc(fpr6,tpr6)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                J6 = tpr6 - fpr6
                ix6 = argmax(J6)
                best_thresh6 = thresholds6[ix6]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                plt.plot(fpr6, tpr6, color='red',
             lw=lw, label='AUC = %0.3f' % auc6)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 6 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n6 >> query/curve in red, with AUC = {} and best threshold = {}'.format(auc6, best_thresh6),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("A6D6MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv")
                os.remove("dfmA1A2A3A4A5A6.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4A5 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4A5.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                fpr5, tpr5, thresholds5 = roc_curve(rocdata['Outcome5'], rocdata['TanimotoCombo5'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                auc5 = auc(fpr5,tpr5)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                J5 = tpr5 - fpr5
                ix5 = argmax(J5)
                best_thresh5 = thresholds5[ix5]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                plt.plot(fpr5, tpr5, color='dodgerblue',
             lw=lw, label='AUC = %0.3f' % auc5)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 5 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n5 >> query/curve in blue, with AUC = {} and best threshold = {}'.format(auc5, best_thresh5),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("A5D5MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv")
                os.remove("dfmA1A2A3A4A5.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3A4 is True:
                rocdata = pd.read_csv("dfmA1A2A3A4.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                fpr4, tpr4, thresholds4 = roc_curve(rocdata['Outcome4'], rocdata['TanimotoCombo4'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                auc4 = auc(fpr4,tpr4)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                J4 = tpr4 - fpr4
                ix4 = argmax(J4)
                best_thresh4 = thresholds4[ix4]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                plt.plot(fpr4, tpr4, color='gold',
             lw=lw, label='AUC = %0.3f' % auc4)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 4 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n4 >> query/curve in yellow, with AUC = {} and best threshold = {}'.format(auc4, best_thresh4),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("A4D4MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv")
                os.remove("dfmA1A2A3A4.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2A3 is True:
                rocdata = pd.read_csv("dfmA1A2A3.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                fpr3, tpr3, thresholds3 = roc_curve(rocdata['Outcome3'], rocdata['TanimotoCombo3'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                auc3 = auc(fpr3,tpr3)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]                
                J3 = tpr3 - fpr3
                ix3 = argmax(J3)
                best_thresh3 = thresholds3[ix3]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                plt.plot(fpr3, tpr3, color='hotpink',
             lw=lw, label='AUC = %0.3f' % auc3)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 3 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n3 >> query/curve in pink, with AUC = {} and best threshold = {}'.format(auc3, best_thresh3),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("A3D3MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv")
                os.remove("dfmA1A2A3.csv") # add a # at the beginning if u wanna keep these files
            elif existA1A2 is True:
                rocdata = pd.read_csv("dfmA1A2.csv")
                plt.figure()
                lw = 2
                fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                fpr2, tpr2, thresholds2 = roc_curve(rocdata['Outcome2'], rocdata['TanimotoCombo2'])
                auc1 = auc(fpr1,tpr1)
                auc2 = auc(fpr2,tpr2)
                J1 = tpr1 - fpr1
                ix1 = argmax(J1)
                best_thresh1 = thresholds1[ix1]
                J2 = tpr2 - fpr2
                ix2 = argmax(J2)
                best_thresh2 = thresholds2[ix2]
                plt.plot(fpr1, tpr1, color='darkorange',
             lw=lw, label='AUC = %0.3f' % auc1)
                plt.plot(fpr2, tpr2, color='mediumturquoise',
             lw=lw, label='AUC = %0.3f' % auc2)
                parentdirectory = os.getcwd()
                for fileQ1 in os.listdir(parentdirectory):
                    if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                        nmfileQ1 = fileQ1.split('_')[2]
                        #print(nmfileQ1)
                        # saving results and writting it in txt file
                        archv = open('out_rocsvalidroc.txt', 'w+')
                        archv.writelines(['\nRocsValidRoc has finished running!',
                                          '\n',
                                          '\nHere are your obtained results:',
                                          '\n',
                                          '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 2 queries that will be summarized in the same order as they are listed in your query molecule file:',
                                          '\n',
                                          '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                          '\n2 >> query/curve in turquoise, with AUC = {} and best threshold = {}'.format(auc2, best_thresh2),
                                          '\n',
                                          '\nThank you for using RocsValidRoc :-)',
                                          '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                        archv.close()
                plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                plt.xlim([-0.02, 1.02])
                plt.ylim([-0.02, 1.02])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('ROC Curve')
                plt.legend(loc="lower right")
                plt.show()
                os.remove("A1D1MergedActInactTCOutcome.csv")
                os.remove("A2D2MergedActInactTCOutcome.csv")
                os.remove("dfmA1A2.csv") # add a # at the beginning if u wanna keep these files
            else:
                if existA1A2 is False:
                    rocdata = pd.read_csv("A1D1MergedActInactTCOutcome.csv")
                    plt.figure()
                    lw = 2
                    fpr1, tpr1, thresholds1 = roc_curve(rocdata['Outcome1'], rocdata['TanimotoCombo1'])
                    auc1 = auc(fpr1,tpr1)
                    J1 = tpr1 - fpr1
                    ix1 = argmax(J1)
                    best_thresh1 = thresholds1[ix1]
                    plt.plot(fpr1, tpr1, color='darkorange',
                 lw=lw, label='AUC = %0.3f' % auc1)
                    parentdirectory = os.getcwd()
                    for fileQ1 in os.listdir(parentdirectory):
                        if fnmatch(fileQ1, '*active*') and fileQ1.endswith("_1.rpt"):
                            nmfileQ1 = fileQ1.split('_')[2]
                            #print(nmfileQ1)
                            # saving results and writting it in txt file
                            archv = open('out_rocsvalidroc.txt', 'w+')
                            archv.writelines(['\nRocsValidRoc has finished running!',
                                              '\n',
                                              '\nHere are your obtained results:',
                                              '\n',
                                              '\nYou used the query molecule file {}'.format(nmfileQ1), '(.sdf or .mol2) with 1 query for which corresponding results are summarized next:',
                                              '\n',
                                              '\n1 >> query/curve in orange, with AUC = {} and best threshold = {}'.format(auc1, best_thresh1),
                                              '\n',
                                              '\nThank you for using RocsValidRoc :-)',
                                              '\nAuthor: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil'])
                            archv.close()
                    plt.plot([-0.02, 1.02], [-0.02, 1.02], color='silver', lw=lw, linestyle='--', label='random')
                    plt.xlim([-0.02, 1.02])
                    plt.ylim([-0.02, 1.02])
                    plt.xlabel('False Positive Rate')
                    plt.ylabel('True Positive Rate')
                    plt.title('ROC Curve')
                    plt.legend(loc="lower right")
                    plt.show()
                    os.remove("A1D1MergedActInactTCOutcome.csv")
                    #print(str('Oops.. something went wrong! Try again!'))

            print(str('Done! Thank you for using RocsValidRoc :-) '))
            print('Author: Guilherme M. Silva \nComputational Laboratory of Pharmaceutical Chemistry (LCQF) \nFCFRP, University of São Paulo, Brazil\n')
    

############################ Answer for option 7 ##############################
# Answer for option 7
    elif resposta == 7:
        cabeçalho('Now exiting RocsValidRoc... Seeya!')
        sleep(5)
        break

    else:
        print('Oooops! Pay attention, you typed something incorrect. Try again!')
        sleep(1)
