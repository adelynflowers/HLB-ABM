CREATE TABLE simulation_clusters
(
    cluster_id smallserial,
    description VARCHAR(255) UNIQUE,
    PRIMARY KEY(cluster_id)
);

CREATE TABLE simulation_sets
(
    set_id smallserial,
    description VARCHAR(255),
    cluster_id smallint,
    PRIMARY KEY(set_id),
    FOREIGN KEY(cluster_id)
        REFERENCES simulation_clusters(cluster_id)
);

CREATE TABLE bio_parameters 
(
    param_id serial,
    max_flush_age smallint,
    flush_emerging smallint,
    egg_adult_transition smallint,
    duration_young_flush smallint,
    proportion_migrating smallint,
    within_row_p real,
    between_row_p real,
    egg_duration smallint,
    nymph_duration smallint,
    shoot_capacity smallint,
    shoot_egg_capacity smallint,
    eggs_per_female_adult smallint,
    transmission_flush_nymph real,
    transmission_adult_flush real,
    latent_period smallint,
    egg_survival_p real,
    adult_survival_p real,
    nymph_min_age_to_infect smallint,
    nymph_min_age_to_be_infected smallint,
    model_duration integer,
    initial_infected_portion real,
    initial_num_psyllids integer,
    invasion_day smallint,
    carrying_capacity integer,
    border_crossing_p real,
    spring_flush_start smallint,
    spring_flush_end smallint,
    summer_flush_start smallint,
    summer_flush_end smallint,
    fall_flush_start smallint,
    fall_flush_end smallint,
    row_length smallint,
    num_rows smallint,
    PRIMARY KEY (param_id)
);

CREATE TABLE experiments
(
    experiment_id serial,
    invasion_modality smallint,
    invasion_grove_id smallint,
    fresh_yield real,
    juice_yield real,
    fruit_price real,
    juice_price real,
    fixed_costs real,
    bio_param_id integer,
    set_id integer,
    PRIMARY KEY(experiment_id),
    FOREIGN KEY(set_id)
        REFERENCES simulation_sets(set_id),
    FOREIGN KEY(bio_param_id)
        REFERENCES bio_parameters(param_id)
);

CREATE TABLE econ
(
    t smallint,
    id smallint,
    costs real,
    returns real,
    profit real,
    hlb_severity real,
    strategy_names text,
    strategy_params text,
    experiment_id integer,
    FOREIGN KEY(experiment_id)
        REFERENCES experiments(experiment_id)
);

CREATE TABLE bio 
(
    t smallint,
    i smallint,
    j smallint,
    num_psyllids integer,
    num_infected_psyllids integer,
    hlb_severity real,
    experiment_id integer,
    FOREIGN KEY(experiment_id)
        REFERENCES experiments(experiment_id)
);


INSERT INTO bio_parameters VALUES(
DEFAULT, 30, 20, 17, 13, 0.4, 0.05, 0.95, 
4, 13, 40, 40, 10, 0.083, 0.3, 15, 0.8614, 0.9847, 17, 6, 
1825, 0.18, 300, 80, 40000, 0.01, 80, 150, 180, 195, 250, 280, 75, 33);