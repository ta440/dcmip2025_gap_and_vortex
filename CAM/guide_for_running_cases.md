## How to run test cases in CAM.

Here are instructions for running the test cases in CAM.


Add discussion of making the new ic file, and defining it as a new option within the namelist, etc. .


1. Move into a local CAM clone.
2. Create a new case, using
   cime/scripts/create_newcase --case [your_casename] --res [dycore_res] --compset FHS94 --run-unsupported --project [project_number] .

   For the dycore_res option:
   SE: ne60_ne60_mg16
   FV3: C192_C192_mg17
   MPAS mpasa60_mpasa60

3. Move into the directory specified by [your_casename]
4. Perform the xmlchanges outlined in xmlchanges.txt
5. Edit the user_nl_cam files based on the template provided for your dycore (SE, FV3, MPAS).
6. Build the case
7. Run the case!

Note, running FV3 with the small Earth modification may require some other changes. Refer to the [] text file for more information.
