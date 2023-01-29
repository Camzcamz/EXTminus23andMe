# Example for using Rscript 

### Usage ##
# Rscript EXT_RG.R 1_INPUT 2_OUTPUT_PREFIX 3_NUMBER_OF_THREADS(optional)

        # - 1_INPUT: delimited text file with
                # column 1: path to munged sumstat files
                # column 2: trait names
                # column 3: sample prevalence (optional)
                # column 4: populaton prevalence (optional) 

### Example ###
Rscript EXT_RG.R Example_input_list.txt OUT_EXAMPLE 3

# Example_input_list.txt: 
        # ./example_data/age_start_smoke.sumstats.gz  age_started_smoking
        # ./example_data/cig_per_day.sumstats.gz  cigarettes_per_day
        # ./example_data/family_relationship_satisfaction.sumstats.gz family_relationship_satisfaction

