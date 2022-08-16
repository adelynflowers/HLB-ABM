from mysql.connector import (connection)
import pandas as pd
print("importing data")

#Add biological params
# add_bio = ("INSERT INTO bio_parameters"
#            "VALUES"
#            "(NULL, 30, 20, 17, 13, 0.4, 0.95, 0.05, 4, 13, 40, 40, 10, 0.083, 0.3, 15, 0.8614, 0.9847, 17, 6, 1825, 0.18, 300, 80, 40000, 0.01, 80, 140, 180, 195, 250, 280, 75, 33);"


#
# 
def getBehavior(row):
    if row['behavior'] == "No_Action":
        return 1
    elif row['behavior'] == "Individual_Spray":
        return 2
    elif row['behavior'] == "Group_Spray":
        return 3
### CREATE CLUSTER
add_cluster = ("INSERT INTO simulation_clusters "
               "(cluster_id, description) "
               "VALUES (NULL, %(descrip)s);")
def createCluster(cursor, description):
    cluster_data = { 'descrip': description}
    cursor.execute(add_cluster, cluster_data)
    return cursor.lastrowid


### CREATE SET
add_set = ("INSERT INTO simulation_sets"
           "(set_id, description, cluster_id)"
           "VALUES (NULL, %(description)s, %(cluster_id)s);")
def createSet(cursor, description, clusterid):
    set_data = {'description': description, 'cluster_id': clusterid}
    cursor.execute(add_set, set_data)
    return cursor.lastrowid

## CREATE EXPERIMENT
add_experiment = ("INSERT INTO experiments VALUES"
                  "(NULL, %(efficacy)s, %(individual)s, %(group)s, %(modality)s, %(invasion_grove)s,"
                  "0.892, 0.770, 14.00, 2.70, 91, 12716, 1825, 15, %(bio_param)s, %(set_id)s, %(focus_grove)s);"
)

def createExperiment(cursor, efficacy, individual, group, modality, grove, set, focus):
    bio_id = 1 #default bio (haven't changed it yet)
    experiment_data = {
        "efficacy": efficacy,
        'individual': individual,
        'group': group,
        'modality': modality,
        'invasion_grove': grove,
        'bio_param': bio_id,
        'set_id': set,
        'focus_grove': focus
    }
    cursor.execute(add_experiment, experiment_data)
    return cursor.lastrowid

## CREATE BIO
add_bio = ("INSERT INTO bio VALUES (%s, %s, %s, %s, %s, %s, %s);")

def createBio(cursor, df, experiment_id):
    df.t = df.t.astype(int)
    df.i = df.i.astype(int)
    df.j = df.j.astype(int)
    df.numPsyllids = df.numPsyllids.astype(int)
    df.numInfectedPsyllids = df.numInfectedPsyllids.astype(int)
    df.hlbSeverity = df.hlbSeverity.astype(float)
    df['experiment_id'] = experiment_id
    cursor.executemany(add_bio, df.values.tolist())

add_econ = ("INSERT INTO econ VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);")
def createEcon(cursor,df_og, experiment_id):
    df = df_og[["t","i","j"]]
    df.t = df.t.astype(int)
    df.i = df.i.astype(int)
    df.j = df.j.astype(int)
    df['behavior'] = df_og.apply(getBehavior, axis=1)
    df.behavior = df.behavior.astype(int)
    df["costs"] = int(0)
    df["returns"] = int(0)
    df["profit"] = int(0)
    df["lambda"] = int(0)
    df["last_NA_ev"] = int(0)
    df["last_IS_ev"] = int(0)
    df["last_GS_ev"] = int(0)
    df["risk_perception"] = int(0)
    df["alpha_perception"] = int(0)
    df["experiment_id"] = experiment_id
    cursor.executemany(add_econ, df.values.tolist())
    
    

#Establish connection
cnx = connection.MySQLConnection(user='sam', password='CitrusABM21',
                                 host='127.0.0.1',
                                 database='citrus')
cursor = cnx.cursor()

basepath = "../output"
#FIRST, BIO CALIBRATES
import os
def importCalibrateFiles(path, focus, set_id):
    print(f"Checking {path}")
    # Create list of bio files
    print(f"Focus is {focus}")
    experiment_ids = {}
    for file in os.listdir(f"{path}"):
        if "bio" in file:
            modality = int(file.split("_")[1][1])
            invasion_grove = int(file.split("_")[2][1:])
            print(f"{file}: m: {modality} g: {invasion_grove} bio")
            df = pd.read_csv(f"{path}/{file}")
            experiment_id = 0
            if f"{modality}_{invasion_grove}" in experiment_ids:
                experiment_id = experiment_ids[f"{modality}_{invasion_grove}"]
            else:
                experiment_id = createExperiment(cursor, (efficacy/100), 60, 21, modality, invasion_grove, set_id, focus)
                experiment_ids[f"{modality}_{invasion_grove}"] = experiment_id
            createBio(cursor, df, experiment_id)
        elif "econ" in file:
            modality = int(file.split("_")[1][1])
            invasion_grove = int(file.split("_")[2][1:])
            print(f"{file}: m: {modality} g: {invasion_grove} econ")
            df = pd.read_csv(f"{path}/{file}")
            experiment_id = 0
            if f"{modality}_{invasion_grove}" in experiment_ids:
                experiment_id = experiment_ids[f"{modality}_{invasion_grove}"]
            else:
                experiment_id = createExperiment(cursor, (efficacy/100), 60, 21, modality, invasion_grove, set_id, focus)
                experiment_ids[f"{modality}_{invasion_grove}"] = experiment_id
            createEcon(cursor, df, experiment_id)



def importBioCalibrate(efficacy):
    #First, create a cluster
    cluster_id = createCluster(cursor, f"Bio calibrate {efficacy}")
    #Each set corresponds to a folder
    top_folder = f"bio_calibrate_{efficacy}"
    for folder in os.listdir(f"{basepath}/{top_folder}"):
        #Create set, name of folder is strategy
        set_id = createSet(cursor, folder, cluster_id)
        # Seperate by focuses
        if folder in ["noaction_nocooperation", "groupaction_allcooperation"]:
            print(f"Starting {top_folder}/{folder}...")
            importCalibrateFiles(f"{basepath}/{top_folder}/{folder}", "N/A", set_id)
        else:
            for grove in os.listdir(f"{basepath}/{top_folder}/{folder}"):
                print(f"Starting {top_folder}/{folder}/{grove}...")
                importCalibrateFiles(f"{basepath}/{top_folder}/{folder}/{grove}", grove, set_id)
            

for efficacy in [65,75,85]:
    importBioCalibrate(efficacy)





cnx.commit()
cursor.close()
cnx.close()