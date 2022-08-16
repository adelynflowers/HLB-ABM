INSERT INTO simulation_clusters VALUES 
(
    NULL,
    'Test cluster'
);

INSERT INTO simulation_sets VALUES
(
    NULL,
    'Test set',
    (SELECT cluster_id 
        FROM simulation_clusters 
        WHERE description='Test cluster')
);

INSERT INTO experiments VALUES
(
    NULL,
    0.90,
    21,
    60,
    1,
    0,
    0.86,
    0.75,
    14.56,
    12.33,
    91,
    12756,
    1825,
    25.50,
    (SELECT set_id 
        FROM simulation_sets
        WHERE description='Test set')
);

INSERT INTO bio VALUES
(
    1,
    21,
    33,
    40000,
    91425,
    0.97,
    (SELECT experiment_id
        FROM experiments
        WHERE invasion_modality=1
        AND invasion_grove_id=0
        AND set_id=
            (SELECT set_id 
            from simulation_sets 
            WHERE description='Test set')
    )
);