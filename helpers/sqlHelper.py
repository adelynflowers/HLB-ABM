import json
import os
import psycopg2
import multiprocessing as mp
import tqdm
import random


add_cluster = ("INSERT INTO simulation_clusters "
               "(cluster_id, description) "
               "VALUES (DEFAULT, %(descrip)s) RETURNING cluster_id;")
def createCluster(cnx, cursor, description):
    cluster_data = { 'descrip': description}
    cursor.execute(add_cluster, cluster_data)
    cnx.commit()
    return cursor.fetchone()[0]


### CREATE SET
add_set = ("INSERT INTO simulation_sets"
           "(set_id, description, cluster_id)"
           "VALUES (DEFAULT, %(description)s, %(cluster_id)s) RETURNING set_id;")
def createSet(cnx, cursor, description, clusterid):
    set_data = {'description': description, 'cluster_id': clusterid}
    cursor.execute(add_set, set_data)
    cnx.commit()
    return cursor.fetchone()[0]

## CREATE EXPERIMENT
add_experiment = ("INSERT INTO experiments VALUES"
                  "(DEFAULT, %(spray_efficacy)s, %(spray_window_individual)s, %(spray_window_group)s, %(invasion_modality)s, %(invasion_grove_id)s,"
                  "%(fresh_yield)s, %(juice_yield)s, %(fruit_price)s, %(juice_price)s, %(planning_length)s, %(fixed_costs)s, %(projection_length)s, %(spray_cost)s, %(bio_param_id)s, %(set_id)s, %(focus_grove)s)"
                  "RETURNING experiment_id;"
)

def createExperiment(cnx, cursor, paramList):
    bio_id = 1 #default bio (haven't changed it yet)
    experiment_data = {
        "spray_efficacy": paramList['spray_efficacy'],
        'spray_window_individual': paramList['spray_window_individual'],
        'spray_window_group': paramList['spray_window_group'],
        'invasion_modality': paramList['invasion_modality'],
        'invasion_grove_id': paramList['invasion_grove_id'],
        'fresh_yield': paramList['fresh_yield'],
        'juice_yield': paramList['juice_yield'],
        'fruit_price': paramList['fruit_price'],
        'juice_price': paramList['juice_price'],
        'planning_length': paramList['planning_length'],
        'fixed_costs': paramList['fixed_costs'],
        'projection_length': paramList['projection_length'],
        'spray_cost': paramList['spray_cost'],
        'bio_param_id': paramList['bio_param_id'],
        'set_id': paramList['set_id'],
        'focus_grove': paramList['focus_grove']
    }
    cursor.execute(add_experiment, experiment_data)
    cnx.commit()
    return cursor.fetchone()[0]

def loadDataFile(cnx, path, tbl):
    fQuery = f"copy {tbl} FROM \'{path}\' DELIMITER \',\' csv header;"
    cursor = cnx.cursor()
    cursor.execute(fQuery)
    cnx.commit()

cnx_a = psycopg2.connect(host="localhost", user="postgres", password="CitrusABM21", database="citrus")
cursor = cnx_a.cursor()

cluster_name = "randomized alpha lambda premium"
cluster_id = createCluster(cnx_a,cursor,cluster_name)
global_bioid = 1


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


econConfig = None 
with open("configs/econConfig.json", "r") as read_file:
    econConfig = json.load(read_file)




#List to hold all of our triples
masterConfigList = []

def getExperimentID(cnx, setid, econConfig, bioConfig, bioid):
    paramList = {
        "spray_efficacy": econConfig['sprayingPopEff'],
        'spray_window_individual': econConfig['individualWindow'],
        'spray_window_group': econConfig['groupWindow'],
        'invasion_modality': bioConfig['invasionModality'],
        'invasion_grove_id': bioConfig['invasionGrove'],
        'fresh_yield': econConfig['freshYield'],
        'juice_yield': econConfig['juiceYield'],
        'fruit_price': econConfig['freshPrice'],
        'juice_price': econConfig['juicePrice'],
        'planning_length': econConfig['planningLength'],
        'fixed_costs': econConfig['costs'],
        'projection_length': econConfig['projectionLength'],
        'spray_cost': econConfig['sprayCost'],
        'bio_param_id': bioid,
        'set_id': setid,
        'focus_grove': 'N/A'
    }
    cursor = cnx.cursor()
    eID = createExperiment(cnx,cursor, paramList)
    cnx.commit()
    return eID

setid_dict = {}
def appendMasterTriples(cnx, econ, biocons, folder, master, single=False, nameSuffix=""):
    if single:
        rndBio = random.randint(0, len(biocons)-1)
        bio_copy = copy.deepcopy(biocons[rndBio])
        econ_copy = copy.deepcopy(econ) 
        econ_copy["outputFilename"] = bio_copy["fileName"].replace("bio", f"econ_{nameSuffix}") 
        bio_copy["fileName"] = econ_copy["outputFilename"].replace("econ", "bio")
        master.append([econ_copy, bio_copy, folder])
    else: 
        cursor = cnx.cursor()
        setid = createSet(cnx, cursor, folder, cluster_id)
        setid_dict[folder] = setid
        for bio in biocons:
            econ_copy = copy.deepcopy(econ)
            bio_copy = copy.deepcopy(bio)
            econ_copy["outputFilename"] = bio_copy["fileName"].replace("bio", "econ")
            master.append([econ_copy, bio_copy, folder])


#One set
# setid = createSet(cnx_a, cursor, "all", cluster_id)
# setid_dict["all"] = setid
# #Randomize
# for _ in range(0,2000):
#     econ_copy = copy.deepcopy(econConfig)
#     for i in range(0,3):
#         for j in range(0,3):
#             econ_copy[f"g{i}{j}_lambda"] = random.random()
#             econ_copy[f"g{i}{j}_alpha"] = random.random() 
#     rndm_bio_idx = random.randint(0, len(configList) - 1)
#     bio_copy = copy.deepcopy(configList[rndm_bio_idx])
#     masterConfigList.append([econ_copy, bio_copy, "all"])

#REDOING ALL
# for lmbda in [0,0.25,0.5,0.75,1]:
#     for alpha in [0,0.25,0.5,0.75,1]:
#         econ_copy = copy.deepcopy(econConfig)
#         for i in range(0,3):
#             for j in range(0,3):
#                 econ_copy[f"g{i}{j}_alpha"] = alpha
#                 econ_copy[f"g{i}{j}_lambda"] = lmbda
#         appendMasterTriples(cnx_a,econ_copy, configList, f"l{int(lmbda*100)}_a{int(alpha*100)}", masterConfigList)
    
#RANDOMIZE
setid = createSet(cnx_a, cursor, "allSet", cluster_id)
setid_dict["all set"] = setid
numSims = 5000
for k in range(0,numSims):
    econ_copy = copy.deepcopy(econConfig)
    for i in range(0,3):
        for j in range(0,3):
            rndA = random.uniform(0,1) 
            rndL = random.uniform(0,1)
            rndP = random.uniform(0, 56)
            econ_copy[f"g{i}{j}_alpha"] = round(rndA,2)
            econ_copy[f"g{i}{j}_lambda"] = round(rndL,2) 
            econ_copy[f"g{i}{j}_premium"] = round(rndP,2)
    appendMasterTriples(cnx_a,econ_copy, configList, "allSet", masterConfigList,True,f"{k}")
setid_dict["allSet"] = setid

def runInstance(config):
    biofileName = config[1]["fileName"][:-4]
    econFileName = config[0]["outputFilename"][:-4]
    os.makedirs(f'configs/{cluster_name}/{config[2]}', exist_ok=True)
    os.makedirs(f'output/{cluster_name}/{config[2]}', exist_ok=True)
    expid = config[3]
    config[0]["outputFilename"] = f"output/{cluster_name}/{config[2]}/{econFileName}.csv"
    config[1]["fileName"] = f"output/{cluster_name}/{config[2]}/{biofileName}.csv"
    config[0]["experimentID"] = expid
    with open(f"configs/{cluster_name}/{config[2]}/{biofileName}.json", "w") as write_file:
        json.dump(config[1], write_file)
    with open(f"configs/{cluster_name}/{config[2]}/{econFileName}.json", "w") as write_file:
        json.dump(config[0], write_file)   
    os.system(f"./serial.out \"configs/{cluster_name}/{config[2]}/{econFileName}.json\" \"configs/{cluster_name}/{config[2]}/{biofileName}.json\"")


revisedMaster = []
for config in masterConfigList:
    expid = getExperimentID(cnx_a, setid_dict[config[2]], config[0], config[1], global_bioid)
    newConfig = copy.deepcopy(config)
    newConfig.append(expid)
    revisedMaster.append(newConfig)



pool = mp.Pool(mp.cpu_count())

for _ in tqdm.tqdm(pool.imap_unordered(runInstance, revisedMaster), total=len(revisedMaster)):
    pass 

pool.close()

basePath = f"/home/instr1/repo/EconABM_accounting/EconABM/output/{cluster_name}"
print("Starting data loading...")
for config in tqdm.tqdm(revisedMaster):
    loadDataFile(cnx_a, f"{basePath}/{config[2]}/{config[0]['outputFilename']}", "econ")




cnx_a.close()



