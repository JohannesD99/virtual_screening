import os
import gzip
import subprocess
import time
import requests

def backup(latest_config, scores):
    localtime = time.asctime(time.localtime(time.time()))

    data = latest_config
    data.insert(0, "Backup from: "+ localtime + "\n")
    data.insert(1, "\n")
    
    for i in scores.items():
        ID = i[0]
        delG = str(i[1])
        data.append(ID + ": " + delG + "\n")

    
    file = open("backup.txt", "w")
    file.writelines(data)
    file.close()

def download(ZINC):
    """Downloads the ZINC15 tranches with the urls given in the 'ZINC' file"""

    zinc_file = open(ZINC, "r")
    urls = zinc_file.readlines()
   
    for url in urls:

        url = url[:-1]
        file_name = url.split("/")
        r = requests.get(url)
        
        with open(os.path.join(screening_folder,file_name[-1]), "wb") as f:
            f.write(r.content) 
            f.close()

    zinc_file.close()

def unzip():
    """Extracts the downloaded 3d files in the screeningfolder and deletes .gz files"""

    ls = os.listdir(screening_folder)
    for f in ls:
        new_file = f[:-3]

        if f[-3:] == ".gz":

            with gzip.open(f, 'rb') as zipfile:
                file_content = zipfile.read()
                unzipped = open(new_file,"wb")
                unzipped.write(file_content)

                zipfile.close()
                unzipped.close()

            os.remove(f)

def pdbqt_split():
    """One .pdbqt file can hold multiple ZINC molecules from its tranch.
    This function splits them into single molecules and returns each as ZINC_ID.pdbqt. 
    The new files will be saved in a folder called 'ligands'"""
    os.mkdir(os.path.join(screening_folder,"ligands"))

    ls = os.listdir(screening_folder)
    for f in ls:

        if f[-6:] == ".pdbqt":
            mols = []

            with open(f, "r") as pdbqt_tranch:
                file_content = pdbqt_tranch.readlines()
                
                for idx, line in enumerate(file_content):

                    if "MODEL" in line:
                        ZINC_ID = file_content[idx + 1]
                        ZINC_ID = ZINC_ID.split()
                        ZINC_ID = ZINC_ID[-1]
                        mols.append((idx, ZINC_ID))
                
                for n, mol in enumerate(mols):

                    head = mol[0]

                    if n + 1 != len(mols):
                        tail = mols[n+1][0]
                    else:
                        tail = 0

                    file_name = mol[1] + ".pdbqt"
                    file_path = os.path.join(screening_folder,"ligands",file_name)

                    with open(file_path, "w") as pdbqt:
                        pdbqt.writelines(file_content[head+1:tail-1]) #removing the MODEL lines
                        pdbqt.close()
            
            pdbqt_tranch.close()
            os.remove(f)

def score(log_name, score_dict):
    """Takes result from logfile and compares it to previous ones"""

    delete = True
    log_file = open(log_name, "r")
    lines = log_file.readlines()
    log_file.close()
    best_score = lines[-10].split()

    ZINC_ID = log_name[5:-8]
    affinity = float(best_score[1])
    
    #add new score, if dict is under specified lenght or if new score is better than the worst one
    if len(score_dict) < topX:
        score_dict[ZINC_ID] = affinity
        score_dict = dict(sorted(score_dict.items(), key=lambda item: item[1]))
        delete = False
          
    else:
        worst = score_dict.popitem()

        if affinity < worst[1]:
            score_dict[ZINC_ID] = affinity
            score_dict = dict(sorted(score_dict.items(), key=lambda item: item[1]))
            delete = False
            os.remove("logs\\" + worst[0] + "_log.txt")
            os.remove("outputs\\" + worst[0] + "_out.pdbqt")

        else:
            score_dict[worst[0]] = worst[1]
    
    #delete reduntant output and score files
    if delete == True:
        os.remove(log_name)
        os.remove("outputs\\" + ZINC_ID + "_out.pdbqt")
    
    return score_dict

def exe_run(resume, vina):
    """takes the path to the vina.exe file and a logical variable as an input. 
    If resume == 'True' the backup.txt file will be used to continue the screening run"""
    
    if resume == True:

        backup = open("backup.txt", "r")
        content = backup.readlines()
        backup.close()

        score_dict = {}
        for i in content[14:]:
            line = i.split()
            ZINC_ID = line[0][:-1]
            score_dict[ZINC_ID] = float(line[1])

        #shortens the list of ligands to dock to the ones remaining
        ligands = os.listdir(os.path.join(screening_folder,"ligands"))
        last_lig = content[3][15:-1]
        index = ligands.index(last_lig)
        ligands = ligands[index+1:]
        
        

    
    else:
        arg_file = "config.txt"
        os.mkdir(os.path.join(screening_folder,"outputs"))
        os.mkdir(os.path.join(screening_folder,"logs"))

        score_dict = {} #for saving scores
        ligands = os.listdir(os.path.join(screening_folder,"ligands"))

    for lig in ligands:
      
        #adds ligandfile and outputfile to config.txt
        args = open(arg_file,'r')
        lines = args.readlines()
        ligand = "ligand=ligands\\" + lig + "\n"
        log = "log=logs\\" + lig[:-6] +"_log.txt\n"
        output = "out=outputs\\" + lig[:-6] + "_out.pdbqt\n"
        lines[1] = ligand
        lines[2] = log
        lines[3] = output

        args.close()
        args = open(arg_file,'w')
        args.writelines(lines)
        args.close()

        #running vina
        subprocess.Popen("cmd")
        command = vina + " --config config.txt"
        p = subprocess.Popen(command)
        p.wait()

        #extracts result from logfile
        score_dict = score(log[4:-1], score_dict)
        backup(lines, score_dict)


    return score_dict

# insert path to screeningfolder
# if a folder of .pdbqt ligands is provided by you, it should be named "ligands" and placed within the screening_folder

global screening_folder
screening_folder = r""

os.chdir(screening_folder)

# insert path to Zinc15 download file which is expected to contain urls
ZINC = r"" 

# path to vina.exe
vina = r""

# declare how many of the highest scoring molecules (based on affinity) you want to obtain as the docking result
global topX
topX = 5        # for example 5

###############################
#initialize config file

if os.path.exists("backup.txt") != True:
    #insert values for all empty arguments. Exhaustiveness ("ex") is optional. Refer to http://vina.scripps.edu/index.html

    rec = ""                                # insert the path to the protein to dock. Remember to correctly prepare it before docking.
    x_coord = ""
    y_coord = ""
    z_coord = ""
    x_size = ""
    y_size = ""
    z_size = ""
    ex = ""
    config = open(os.path.join(screening_folder,"config.txt"),"w")
    content = ["receptor=" + rec,
                "ligand=",
                "log=",
                "out=",
                "center_x=" + x_coord,
                "center_y=" + y_coord,
                "center_z=" + z_coord,
                "size_x=" + x_size,
                "size_y=" + y_size,
                "size_z=" + z_size,
                "exhaustiveness=" + ex,
                "\n"
                ]

    content = "\n".join(content)
    config.writelines(content)
    config.close()

###############################

# download and preparation           
download(ZINC)
unzip()
pdbqt_split()

# run the virtual screening
# set "resume" True if you want to continue from a breakpoint ("backup.txt" file is expected in screening_folder)
# for a new run set False
resume = False
score = exe_run(resume, vina)

# write your own lines to process the dictionary with the molecules and respective scores in it
# or simply open "backup.txt" it will contain the best to worst in descending order
