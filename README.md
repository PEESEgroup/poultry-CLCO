# Combined Techno-Economic and Life-Cycle Optimization of Poultry Manure Management Pathways

This readme document is intended to aid others in using the model and reproducing results from the paper "Consequential life cycle optimization of poultry manure management technologies in a food-energy-water-waste nexus" under consideration at Environmental Science & Technology.

# Running the model
## Installation
Download the github repository.  To run any CLCO scenario, open the CLCO.py file in your IDE of choice, and follow the steps in the generating results section below.  To run the FLP problem, run the CLCO_FLP.py file.  

It is unnecessary to run any other python file in this repository.

CLCO_Data.py holds parameter information for the program, and sources for this data can be found in CLCO_TEA_LCA_Data.xlsx in the data analysis folder, values which are also repeated throughout the SI.  

To run this model with your own parameters, simply change the relevant parameters in the CLCO_Data.py file.

## Necessary software

 - Python 3.8 or greater
 - Gurobi Optimizer (version 9.5.2 build v9.5.2rc0 (win64) was used to generate the results)
 - Python dependencies 

## Model runtime

A Dell XPS 15 7590 with Windows 10 Home and Intel(R) Core(TM) i-7-9750H CPU @2.60GHz, 2592 Mhz, 6 Core(s), 12 Logical Processor(s) can run an optimization for 1 county in 10 seconds or less.  The calculation time necessary for each Pareto front is about 1-2 hours.

## Generating results

To run the model, simply choose a scenario from the list below and edit the line of code in the main method of CLCO.py.  This line of code is located toward the bottom of the python file.  

To run scenario 3:
```    S = [3]```

To run scenario 3, 428, and 1503:
```    S = [3, 428, 1503]```

### NPV and LCA figures

Run scenarios 2, 5, 6, 7, 8 for direct land application, pyrolysis, HTL, HTC, and AD for figures 2, 3, and 5.
Run scenarios 2, 5, 8, 10, 11, 12, and 10103 and then extract the entries for Onondaga county to create figures S1 and S2.  The normalization technology was chosen based on what technology would have the lowest range of LCA impacts.

### FLP

Figure 4 was created by running the CLCO_FLP.py file.  This program does a full search over all maximum transportation distances and number of facilities.  Looking at the reported objective function values (likely in the SI) it was easy to determine the minimum transportation distance for each number of facilities.  Then, the solutions were manually adjusted as a few counties with very little manure were mislabeled due to a high mipgap to decrease computation time.

### Pareto Fronts

To find the 3-D Pareto fronts, run scenarios 4501, 4502, 4503, 4511, 4512, and 4513.  

To find the 2-D Pareto fronts in the SI, run scenarios 1501, 1502, 1503, 1511, 1512, 1513 for the CLCA GWP/NPV Pareto Front and scenarios 2501, 2502, 2503, 2511, 2512, 2513 for the ALCA GWP/NPV Pareto Front.  Pareto fronts are automatically plotted using matplotlib and can be found in the folder of the scenario number.
Pareto fronts found in the figures in the manuscript are manually adjusted at the anchor points.

### Sensitivity analysis

To replicate the sensitivity analysis, run scenarios 10003, 10103, 10203 for NPV max, tradeoff, and GWP min respectively.  Modify CLCO_Data.py to change key parameters.

### Changing model location

To change the model location, update FEEDSTOCK_SUPPLY and INTRA_COUNTY_TRANSPORT_DISTANCE in CLCO_Data.py.  The data is stored in a list, so the nth county in the feedstock supply list corresponds to the nth county in the transportation distance list

**It should be noted that the model has not been tested for counties with more than 6,000 tons of poultry manure.  Extremely large inputs may require modification to the piecewise linear approximtion list lpa_xvals and bounded variables**

## Excel plotting

For any scenario that does not have an existing S[Scenario Number].xlsx file, one can be created from the STemplate.xlsx.

- First, delete the existing TEA and LCA tabs, as they are being replaced by references to new directories.
- Then go to the data tab in the ribbon
- Click on "Get Data" on the left part of the ribbon, then select "from file" and then "from folder".
- Next, navigate in the data file directory to the folder with the scenario number you want
- Enter that folder, and then select TEA.
- In the pop-up window, on the combine drop-down, select combine and load, and then click "OK".
- Rename the sheet to TEA from "TEA (2)" or whatever it loads in as.
- Repeat the above steps, except select the LCA folder instead.
- Some simple figures should already be populated, although they may need additional modification (i.e. to include all 62 counties in the figure and the data table), though this shouldn't be an issue for any figure that appears in the text or the SI.

## Scenario List
Not every scenario in the scenario list below was used to create figures for this model or validate the model.
```
    SCENARIO LIST:
    1: AD + CHP w/ disposal of digestate, NPV max, county level
    2: Direct Land Application, NPV max, county level
    3: Optimal County Level Results, NPV max, county level
    4: FLP, NPV max, county level
    5: Pyrolysis + CHP, NPV max, county level
    6: HTL + CHP, NPV max, county level
    7: HTC + CHP, NPV max, county level
    8: AD + CHP, NPV max, county level
    9: AD + Pyrolysis + CHP, NPV max, county level
    10: HTC w/ selling hydrochar, NPV max, county level
    11: HTC w/ DLA, NPV max, county level
    12: HTL w/ combustion of hydrochar, NPV max, county level
    50: Pyrolysis + CHP GWP min, county level
    51: Pyroylsis + CHP Onondaga county Pareto Front min GWP max NPV
    52: Pyroylsis + CHP Jefferson county Pareto Front min GWP max NPV
    421: AD + CHP w/ disposal of digestate, NPV max, two optimal facilities
    422: Direct Land Application, NPV max, two optimal facilities
    423: Optimal County Level Results, NPV max, two optimal facilities
    424: FLP, NPV max, two optimal facilities
    425: Pyrolysis + CHP, NPV max, two optimal facilities
    426: HTL + CHP, NPV max, two optimal facilities
    427: HTC + CHP, NPV max, two optimal facilities
    428: AD + CHP, NPV max, two optimal facilities
    429: AD + Pyrolysis, NPV max, two optimal facilities
    431: AD + CHP w/ disposal of digestate, NPV max, three optimal facilities
    432: Direct Land Application, NPV max, three optimal facilities
    433: Optimal County Level Results, NPV max, three optimal facilities
    434: FLP, NPV max, three optimal facilities
    435: Pyrolysis + CHP, NPV max, three optimal facilities
    436: HTL + CHP, NPV max, three optimal facilities
    437: HTC + CHP, NPV max, three optimal facilities
    438: AD + CHP, NPV max, three optimal facilities
    439: AD + Pyrolysis, NPV max, three optimal facilities
    441: AD + CHP w/ disposal of digestate, NPV max, four optimal facilities
    442: Direct Land Application, NPV max, four optimal facilities
    443: Optimal County Level Results, NPV max, four optimal facilities
    444: FLP, NPV max, four optimal facilities
    445: Pyrolysis + CHP, NPV max, four optimal facilities
    446: HTL + CHP, NPV max, four optimal facilities
    447: HTC + CHP, NPV max, four optimal facilities
    448: AD + CHP, NPV max, four optimal facilities
    449: AD + Pyrolysis, NPV max, four optimal facilities
    1001:, CLCA,  Direct Land Application, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1002:, CLCA,  Direct Land Application, Pareto Front NPV max GWP min, Onondaga county
    1003:, CLCA,  Direct Land Application, Pareto Front NPV max GWP min, Jefferson county
    1011:, CLCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1012:, CLCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1013:, CLCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1101:, CLCA,  Pyrolysis, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1102:, CLCA,  Pyrolysis, Pareto Front NPV max GWP min, Onondaga county
    1103:, CLCA,  Pyrolysis, Pareto Front NPV max GWP min, Jefferson county
    1111:, CLCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1112:, CLCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1113:, CLCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1201:, CLCA,  HTL, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1202:, CLCA,  HTL, Pareto Front NPV max GWP min, Onondaga county
    1203:, CLCA,  HTL, Pareto Front NPV max GWP min, Jefferson county
    1211:, CLCA,  HTL, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1212:, CLCA,  HTL, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1213:, CLCA,  HTL, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1301:, CLCA,  HTC, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1302:, CLCA,  HTC, Pareto Front NPV max GWP min, Onondaga county
    1303:, CLCA,  HTC, Pareto Front NPV max GWP min, Jefferson county
    1311:, CLCA,  HTC, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1312:, CLCA,  HTC, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1313:, CLCA,  HTC, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1401:, CLCA,  AD, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1402:, CLCA,  AD, Pareto Front NPV max GWP min, Onondaga county
    1403:, CLCA,  AD, Pareto Front NPV max GWP min, Jefferson county
    1411:, CLCA,  AD, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1412:, CLCA,  AD, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1413:, CLCA,  AD, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1501:, CLCA,  All, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1502:, CLCA,  All, Pareto Front NPV max GWP min, Onondaga county
    1503:, CLCA,  All, Pareto Front NPV max GWP min, Jefferson county
    1511:, CLCA,  All, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1512:, CLCA,  All, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1513:, CLCA,  All, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2001:, ALCA,  Direct Land Application, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2002:, ALCA,  Direct Land Application, Pareto Front NPV max GWP min, Onondaga county
    2003:, ALCA,  Direct Land Application, Pareto Front NPV max GWP min, Jefferson county
    2011:, ALCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2012:, ALCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2013:, ALCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2101:, ALCA,  Pyrolysis, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2102:, ALCA,  Pyrolysis, Pareto Front NPV max GWP min, Onondaga county
    2103:, ALCA,  Pyrolysis, Pareto Front NPV max GWP min, Jefferson county
    2111:, ALCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2112:, ALCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2113:, ALCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2201:, ALCA,  HTL, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2202:, ALCA,  HTL, Pareto Front NPV max GWP min, Onondaga county
    2203:, ALCA,  HTL, Pareto Front NPV max GWP min, Jefferson county
    2211:, ALCA,  HTL, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2212:, ALCA,  HTL, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2213:, ALCA,  HTL, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2301:, ALCA,  HTC, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2302:, ALCA,  HTC, Pareto Front NPV max GWP min, Onondaga county
    2303:, ALCA,  HTC, Pareto Front NPV max GWP min, Jefferson county
    2311:, ALCA,  HTC, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2312:, ALCA,  HTC, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2313:, ALCA,  HTC, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2401:, ALCA,  AD, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2402:, ALCA,  AD, Pareto Front NPV max GWP min, Onondaga county
    2403:, ALCA,  AD, Pareto Front NPV max GWP min, Jefferson county
    2411:, ALCA,  AD, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2412:, ALCA,  AD, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2413:, ALCA,  AD, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2501:, ALCA,  All, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2502:, ALCA,  All, Pareto Front NPV max GWP min, Onondaga county
    2503:, ALCA,  All, Pareto Front NPV max GWP min, Jefferson county
    2511:, ALCA,  All, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2512:, ALCA,  All, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2513:, ALCA,  All, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    3000: 1 pyrolysis plant for the whole state - max NPV
    3001: 2 pyrolysis plants for the whole state - max NPV
    3002: pyrolysis plant in every county, calculate biochar break-even prices - max NPV
    3003: plant in every county, calculate biochar break-even prices - max NPV
    3100: 1 pyrolysis plant for the whole state - max NPV and min GWP
    3101: 2 pyrolysis plants for the whole state - max NPV and min GWP
    3102: pyrolysis plant in every county, calculate biochar break-even prices - max NPV and min GWP
    3103: plant in every county, calculate biochar break-even prices - max NPV and min GWP
    3200: 1 pyrolysis plant for the whole state - min GWP
    3201: 2 pyrolysis plants for the whole state - min GWP
    3202: pyrolysis plant in every county, calculate biochar break-even prices - min GWP
    3203: plant in every county, calculate biochar break-even prices - min GWP
    4501:, CLCA,  All, Pareto Front FE min GWP min, largest facility from 2 facility FLP
    4502:, CLCA,  All, Pareto Front FE min GWP min, Onondaga county
    4503:, CLCA,  All, Pareto Front FE min GWP min, Jefferson county
    4511:, ALCA,  All, Pareto Front FE min GWP min, largest facility from 2 facility FLP
    4512:, ALCA,  All, Pareto Front FE min GWP min, Onondaga county
    4513:, ALCA,  All, Pareto Front FE min GWP min, Jefferson county
    10003: first calculated optimal plants at NPV max in every county, then implemented constraints on those plants for sensitivity analysis
    10103: first calculated optimal plants at a tradeoff in every county, then implemented constraints on those plants for sensitivity analysis
    10203: first calculated optimal plants at GWP min in every county, then implemented constraints on those plants for sensitivity analysis
```

