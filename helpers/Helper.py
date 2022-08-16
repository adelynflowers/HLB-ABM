
import json

#Open default config file
config_file = None 

with open("configs/bioConfig.json", "r") as read_file:
    config_file = json.load(read_file)

#Create a list of bio configs
import datetime
import copy 
filePrefix = f'{datetime.datetime.now()}'.replace(' ', '-').replace(':','.')
configList = []
for modality in range(1,7):
    for grove_i in range(0,3):
        for grove_j in range(0,3):
            config = json.loads(json.dumps(config_file))
            config["invasionModality"] = modality
            groveID = (grove_i * 10) + grove_j 
            config["invasionGrove"] = groveID 
            config["fileName"] = f'{filePrefix}_m{modality}_g{groveID}_bio.csv'
            configList.append(config)

#generate two econ configs

#Everybody set to group spray
# groupEcon = None
# with open("configs/econConfig.json", "r") as read_file:
#     groupEcon = json.load(read_file)

#everybody set to no action
# noneEcon = copy.deepcopy(groupEcon)
# for i in range(0,3):
#     for j in range(0,3):
#         noneEcon[f"g{i}{j}_behavior"] = 0

econConfig = None 
with open("configs/econConfig.json", "r") as read_file:
    econConfig = json.load(read_file)




#List to hold all of our triples
masterConfigList = []

def appendMasterTriples(econ, biocons, folder, master):
    for bio in biocons:
        econ_copy = copy.deepcopy(econ)
        bio_copy = copy.deepcopy(bio)
        econ_copy["outputFilename"] = bio_copy["fileName"].replace("bio", "econ")
        master.append([econ_copy, bio_copy, folder])



#Lambda 
lambda_vals = [0.25, 0.5, 0.75, 1]
alpha_vals = [0.25, 0.5, 0.75, 1]
for lambdaVal in lambda_vals:
    for alphaVal in alpha_vals:
        gEcon = copy.deepcopy(econConfig)
        for i in range(0,3):
            for j in range(0,3):
                gEcon[f"g{i}{j}_lambda"] = lambdaVal
                gEcon[f"g{i}{j}_alpha"] = alphaVal 
        appendMasterTriples(gEcon, configList, f"l{int(lambdaVal*100)}_a{int(alphaVal*100)}", masterConfigList)

#No action, no cooperation
# appendMasterTriples(noneEcon, configList, "noaction_nocooperation", masterConfigList)

# #No action, all cooperation
# for i in range(0,3):
#     for j in range(0,3):
#         gEcon = copy.deepcopy(groupEcon)
#         gEcon[f"g{i}{j}_behavior"] = 0
#         appendMasterTriples(gEcon, configList, f"noaction_allcooperation/g{i}{j}", masterConfigList)

# #Individual action, no cooperation
# for i in range(0,3):
#     for j in range(0,3):
#         gEcon = copy.deepcopy(noneEcon)
#         gEcon[f"g{i}{j}_behavior"] = 1
#         appendMasterTriples(gEcon, configList, f"individualaction_nocooperation/g{i}{j}", masterConfigList)

# #Individual action, all cooperation
# for i in range(0,3):
#     for j in range(0,3):
#         gEcon = copy.deepcopy(groupEcon)
#         gEcon[f"g{i}{j}_behavior"] = 1
#         appendMasterTriples(gEcon, configList, f"individualaction_allcooperation/g{i}{j}", masterConfigList)




# #Group action, no cooperation for spray and windows
# for i in range(0,3):
#     for j in range(0,3):
#         gEcon = copy.deepcopy(noneEcon)
#         gEcon[f"g{i}{j}_behavior"] = 2
#         appendMasterTriples(gEcon, configList, f"groupaction_nocooperation/g{i}{j}", masterConfigList)

#Group action, all cooperation
#appendMasterTriples(groupEcon, configList, "groupaction_allcooperation", masterConfigList)
#Window spray combos
# for window in sprayWindows:
#     for efficacy in sprayEfficacy:
#         gEcon = copy.deepcopy(groupEcon)
#         gEcon["groupWindow"] = window 
#         gEcon["sprayingPopEff"] = efficacy 
#         appendMasterTriples(gEcon, configList, f"group_w{window}_e{int(efficacy*100)}", masterConfigList)


#And now we run each of them

def runInstance(config):
    biofileName = config[1]["fileName"][:-4]
    econFileName = config[0]["outputFilename"][:-4]
    os.makedirs(f'configs/initial_info/{config[2]}', exist_ok=True)
    os.makedirs(f'output/initial_info/{config[2]}', exist_ok=True)
    config[0]["outputFilename"] = f"output/initial_info/{config[2]}/{econFileName}.csv"
    config[1]["fileName"] = f"output/initial_info/{config[2]}/{biofileName}.csv"
    with open(f"configs/initial_info/{config[2]}/{biofileName}.json", "w") as write_file:
        json.dump(config[1], write_file)
    with open(f"configs/initial_info/{config[2]}/{econFileName}.json", "w") as write_file:
        json.dump(config[0], write_file)   
    os.system(f"./serial.out \"configs/initial_info/{config[2]}/{econFileName}.json\" \"configs/initial_info/{config[2]}/{biofileName}.json\"")



import os
#for config in masterConfigList:
#    runInstance(config)


## RUN SIMULATIONS

import multiprocessing as mp
import tqdm
pool = mp.Pool(mp.cpu_count())

for _ in tqdm.tqdm(pool.imap_unordered(runInstance, masterConfigList), total=len(masterConfigList)):
    pass 

pool.close()


## CREATE OUTPUT FILES
# import pandas as pd
# combined_ev = None
# evFiles = []
# csv_files = os.listdir("output")
# for file in csv_files:
#     outputFilename = file[:-4] + "-ev.csv" 
#     os.system(f"Rscript Analysis/CreateEVFile.R output/{file} output/{outputFilename}")
#     evFiles.append(f"output/{outputFilename}")

# #make mega file
# empty = True
# for file in os.listdir("output"):
#     if file.endswith("-ev.csv"):
#         df = pd.read_csv(f"output/{file}")
#         if empty:
#             combined_ev = df
#             empty = False 
#         else:
#             combined_ev = combined_ev.append(df)

# combined_ev.to_csv("Analysis/combined-ev.csv", index=False)


import smtplib,ssl

from email.message import EmailMessage
port = 465
password = "CitrusABM21"
context = ssl.create_default_context()
message = EmailMessage()
message.set_content("The tasks have finished running, please check the output folder for results.")
message['Subject'] = "Your HPC tasks have finished running"
message['From'] = "citrusabmalerts@gmail.com"
message['To'] = "shaynes@csus.edu"

with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
    server.login("citrusabmalerts@gmail.com", password)
    server.send_message(message)

