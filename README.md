# RocsValidRoc
Automatic and interactive way to perform retrospective validation of ROCS queries through ROC curve analysis.

# What does this script do? #
 RocsValidRoc fully automatize in an interactive way the validation procedure of using one (or various) molecule(s) as query(ies) in 3D shape overlaping (similarity search) using ROCS. For this, it automatically builds ROC curve(s) allowing you to check its metrics and to assess your query(ies) performance in ROCS retrospective exercise considering two databases of known actives and decoys.
<br />  It does the same as vROCS, but in command line (not in GUI). It also allows to set the desired number of threads in MPI flag (mpi_np). In particular, it automatically generates merged ROC curves to simultaneously validate multi-queries in ROCS.
<br />  Furthermore, it is very useful if you only wish to use several databases for a given query(ies) in a ROCS VS campaign. You will not have to type every single code line to run ROCS, for each database, as it is usually required in traditional procedure (for this, try ´Run only ROCS´).
              
# What do you need? #
* OpenEye applications installed and regularly working in your machine (including, obviously, ROCS), with appropriate OpenEye valid license (see https://docs.eyesopen.com/ for more info) 
* Python modules: subprocess, shutil, os, sys, smtplib, pandas,  matplotlib, numpy, time, fnmatch, datetime, scikit-learn. In general, they all can be easily installed by typing pip install.
             
# Place in folder: #
* RocsValidRoc.py script
* 1 ´actives´ and 1 ´decoys´ databases files previously processed by OMEGA conformer generation that you wish to run ROCS validation. We highly recommend using default names for databases: [omega]actives.oeb.gz and [omega]decoys.oeb.gz
* query file (in sdf or mol2 or sq). If you wish to use more than 1 query, you should merge them all into one single molecule file (in sdf or mol2 or sq) and inform it as the chosen query. Please note that maximum of 10 multi-queries (within 1 query file) is allowed in this current version.
* note: if you´re running only ´Build single OR multi query ROC curve from ROCS results´ you should place this script in folder containing ROCS results (in particular .rpt files) and then run it.
              
# Some keynotes: #
* Avoid running it in folders (path or directories) named with SPACES on it!
* RocsValidRoc will only list DATABASEs files containing the terms *actives* or *decoys* in corresponding filenames, within the current working directory (and sub-directories)! Moreover, it should contain term in brackets, e.g. [omega]
* If you´re using multi-QUERY (into one single molecule file, in sdf or mol2), results for each query will be in same order as molecules in such file, i.e. first molecule in query-file will be hits_1 and first curve, second hits_2 and second curve, so on... ROC curve output txt file will also help you to track your queries.
* RocsValidRoc will create a sub-directory that you set as your project name at your current working directory, and it´ll send all results there. In addition, it will generate time.txt file with start/end of ROCS run, plus output_rocsvalidroc.txt with all that you wanna know about your built ROC curves.
* Besides AUC values, we provide the optimal threshold in a given ROC curve, which is calculated following https://en.wikipedia.org/wiki/Youden%27s_J_statistic
* RocsValidRoc must be used with databases previously prepared by OMEGA conformer generation.
* Also, it is compatible to be used after a regular ROCS run (for this, try ´Build single OR multi query ROC curve from ROCS results´). However, RocsValidRoc already comes with ROCS function implemented. ** Just remember, if you use ROCS separetely/previous to RocsValidRoc, name your OMEGA processed databases with terms ´actives´ and ´decoys´ and then you may proceed to run RocsValidRoc **
* note: this homemade script was developed by me for personal purposes (especially to practice python skills) and it does not has a link to OpenEye.
              
> Author: Guilherme M. Silva (guimsilva@gmail.com)
