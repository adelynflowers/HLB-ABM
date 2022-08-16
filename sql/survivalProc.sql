/*
Step 1: Get all bio data with the requested experiment ID
OPTIONAL STEP FOR FOCUS GROVES
----
Step 2: filter I and J to be within the values of that grove
    - the corresponding bio params entry has row_length and num_rows for calculating bounds
----
Step 3: Assign grove IDs
Step 4: If a grove ever has hlb_severity > 0, pick the minimum day that happens 
Step 5: If not, pick maximum t
Step 6: Assign E accordingly
Step 7: Insert into survival

*/
CREATE PROCEDURE Create_Survival_Experiment @ID INT UNSIGNED,
AS 



GO;



