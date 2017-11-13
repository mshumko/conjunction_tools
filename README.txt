Conjunction tool kit was written by Mykhaylo Shumko. Please direct all 
comments, questions, or bugs to msshumko at gmail.com. 

To run the conjunction calculator, read the __init__ function for all of the 
required and optional inputs. 

The basic workflow for MagneticConjunctionCalc class is as follow.

1. IF the data is in JSON headed ASCII format. 
        THEN: Call the class with paths to magnetic ephemeris data. Specify 
              conjunction criteria.
        ELSE: LOAD in the data yourself into a dictionary with time, L and MLT.
              Specify the new key's for that data if necessary. Specify 
              conjunction criteria.

2. Sync the two data sets using:
   A: the periodic_data_indicies(), which assumes both datasets are continous 
      and have the same cadence. <<Implement an indexing array similar to the
      non-periodic case>>
   B: the non_periodic_data_indicies(), which does not assume the following,
      and will match up times, with a little bit of a time threshold, tThresh.

3. Filter the data to either the periodic indicies, or non period indicies.

4. Use the filtered data to find the change in L and MLT via 
   calc_magnetic_seperation().

5. Find out which indicies satisfy the conjunction criteria via
   calc_indicies_L_MLT_conjunction(), and lower_L_bound().

6. Find the duration of the conjunctions via calc_conjunction_duration(). The 
   plotting/output data will be created here.

8. Plot/save the data to a JSON headed ASCII via plotLMLTandDLDMLT() or 
   save_to_file().

