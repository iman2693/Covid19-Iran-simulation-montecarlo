import sys
from numpy import random as rn
import random
import matplotlib 
import matplotlib.pyplot as plt 
def Simulation(days=300 , nd=30 , Rt=None , muT=4 , sizeV=1 , limit=1000000 , pp=0.001 , n0=1 ):
    # days: observation period
    # nd: simulation period
    # Rt = rr  # infection rate pattern
    #  muT  is the mean time an infected person will transmit the virus to (i.e., infect) another person.
    # We assume that the independence among those ones being infected.  The default value is set as muT = 4 (days).
    # sizeV: the dispersion parameter so that variance = mu + mu^2/size. The default value is set as sizeV =1.
    # limit: the target/study population size
    # pp: the proportion of people with immunity in the population
    # n0: the initial number of infectious persons.  
    # The default setting assumes one virus carrier/infectious person in the beginning, i.e., n0=1.
    kk = [0 for i in range(days)] # kk: daily new cases
    atrisk =[0 for i in range(days)] # atrisk: number of active cases each day; simulation period of nn days
    tt = 0   # the cumulative total number of confirmed cases. 
    if nd > len(Rt):
        print("The length of Rt should not be smaller than nd.")
        sys.exit(0)
    stoplimit = limit*(1-pp)
    nk = n0   # The initial number of existing infectious persons.  
    # there must be a first patient to kick off the transmission process! 
    #------ First Day Of Simulation ------
    for k in range(nk):
        if tt>stoplimit:
            Rt[0]=0.001
        ni = rn.poisson(Rt[0],1)[0]    # how many people will be infected by this existing virus carrier person.
        imuind = rn.choice(2,1,True,[1-pp,pp])[0] # if people with immunity ni=0 
        if(imuind==1):
            ni=0
        tt=tt+ni
        if(ni > 0):
            tk=[0 for i in range(ni)]   
            for i in range(ni):
                tk[i]= rn.negative_binomial(1,sizeV/(sizeV + muT),size=round(sizeV))[0]+1  # this is the nth day on which a new case occurs
                kk[tk[i]-1] = kk[tk[i]-1] + 1
            pastevent =[1 for i in range(max(tk)-1)]+[0 for i in range((days-max(tk)+1))] 
            atrisk = [sum(i) for i in zip(atrisk, pastevent)]   #atrisk = atrisk + pastevent
    #----------------------------------------

    #------ Day 2 to nd ---------
    for j in range(1,nd):
        nk = kk[j-1]    # this is the number of people newly infected (i.e., new cases) on (j-1)th day
        if(nk > 0):
            for k in range(nk):
                if(tt>stoplimit): 
                    Rt[j]=0.001
                ni = rn.poisson(Rt[j],1)[0]     # how many people will be infected by this existing virus carrier person.
                imuind = rn.choice(2,1,True,[1-pp,pp])[0]  # This Person is immunity or not  1= immunity        0= not immunity
                if(imuind==1):                  # if This person is immunity , it can not transmit the disease.
                    ni=0
                tt=tt+ni
                if(ni > 0):
                    tk=[0 for i in range(ni)]
                    for i in range(ni):
                        tk[i] = rn.negative_binomial(1,sizeV/(sizeV + muT),size=round(sizeV))[0]+1+j  # this is the nth day on which a new case occurs
                        kk[tk[i]-1] = kk[tk[i]-1] + 1
                    pastevent = [0 for l in range(j-1)]+[1 for l in range(max(tk)+1-j)]+[0 for l in range(days-max(tk))]
                    atrisk = [sum(i) for i in zip(atrisk, pastevent)]
    return [atrisk,kk,tt] # riskpopu = atrisk, dailynew = kk, total=tt
# --------------------------------------
# observed number of confirmed infection cases in Iran over the period of 1 March to 29 April 2020  
# ir1=[978 , 1501 , 2336 , 2922 , 3513 , 4747 , 5823 , 6566 , 7161 , 8042 , 9000 , 10075 , 11364 ,
#  12729 , 13938 , 14991 , 16169 , 17361 , 18407 , 19644 , 20610 , 21638 , 23049 , 24811 , 27017 ,
#   29406 , 32332 , 35408 , 38309 , 41495 , 44605 , 47593 , 50468 , 53183 , 55743 , 58226 , 60500 ,
#    62589 , 64586 , 66220 , 68192 , 70029 , 71686 , 73303 , 74877 , 76389 , 77995 , 79494 , 80868 ,
#    82211 , 83505 , 84802 , 85996 , 87026 , 88194 , 89328 , 90481 , 91472 , 92584 , 93657]
# Number of new cases in Iran over the period of 1 March to 29 April 2020  
ir1_newcases=[385 , 523 , 835 , 586 , 591 , 1234 , 1076 , 743 , 595 , 881 , 958 , 1075 , 1289 ,
 1365 , 1209 , 1053 , 1178 , 1192 , 1046 , 1237 , 966 , 1028 , 1411 , 1762 , 2206 , 2389 , 2926 ,
  3076 , 2901 , 3186 , 3110 , 2988 , 2875 , 2715 , 2560 , 2483 , 2274 , 2089 , 1997 , 1634 , 1972 ,
   1837 , 1657 , 1617 , 1574 , 1512 , 1606 , 1499 , 1374 , 1343 , 1294 , 1297 , 1194 , 1030 , 1168 ,
    1134 , 1153 , 991 , 1112 , 1073]

# observed number of confirmed infection cases in Iran over the period of 30 April to 28 June 2020
# ir2=[94640 , 95646 , 96448 , 97424 , 98647 , 99970 , 101650 , 103135 , 104691 , 106220 , 107603 ,
#  109286 , 110767 , 112725 , 114533 , 116635 , 118392 , 120198 , 122492 , 124603 , 126949 , 129341 ,
#   131652 , 133521 , 135701 , 137724 , 139511 , 141591 , 143849 , 146668 , 148950 , 151466 , 154445 ,
#    157562 , 160696 , 164270 , 167156 , 169425 , 171789 , 173832 , 175927 , 177938 , 180156 , 182525 , 
#    184955 , 187427 , 189876 , 192439 , 195051 , 197647 , 200262 , 202584 , 204952 , 207525 , 209970 ,
#     212501 , 215096 , 217724 , 220180 , 222669]
# Number of new cases in Iran over the period of 30 April to 28 June 2020 
ir2_newcases=[983 , 1006 , 802 , 976 , 1223 , 1323 , 1680 , 1485 , 1556 , 1529 , 1383 , 1683 , 1481 ,
 1958 , 1808 , 2102 , 1757 , 1806 , 2294 , 2111 , 2346 , 2392 , 2311 , 1869 , 2180 , 2023 , 1787 , 2080 ,
  2258 , 2819 , 2282 , 2516 , 2979 , 3117 , 3134 , 3574 , 2886 , 2269 , 2364 , 2043 , 2095 , 2011 , 2218 ,
   2369 , 2430 , 2472 , 2449 , 2563 , 2612 , 2596 , 2615 , 2322 , 2368 , 2573 , 2445 , 2531 , 2595 , 2628 , 2456 , 2489]

#observed number of confirmed infection cases in Iran over the period of 29 June to 27 Aguest
# ir3=[225205 , 227662 , 230211 , 232863 , 235429 , 237878 , 240438 , 243051 , 245688 , 248379 , 250458 ,
#  252720 , 255117 , 257303 , 259652 , 262173 , 264561 , 267061 , 269440 , 271606 , 273788 , 276202 , 278827 ,
#   281413 , 284034 , 286523 , 288839 , 291172 , 293606 , 296273 , 298909 , 301530 , 304204 , 306752 , 309437 ,
#    312035 , 314786 , 317483 , 320117 , 322567 , 324692 , 326712 , 328844 , 331189 , 333699 , 336324 , 338825 ,
#     341070 , 343203 , 345450 , 347835 , 350279 , 352558 , 354764 , 356792 , 358905 , 361150 , 363363 , 365606 , 367796]
# Number of new cases in Iran over the period of 29 June to 27 Aguest
ir3_newcases=[2536 , 2457 , 2549 , 2652 , 2566 , 2449 , 2560 , 2613 , 2637 , 2691 , 2079 , 2262 , 2397 , 2186 , 2349 ,
 2521 , 2388 , 2500 , 2379 , 2166 , 2182 , 2414 , 2625 , 2586 , 2621 , 2489 , 2316 , 2333 , 2434 , 2667 , 2636 , 2621 ,
  2674 , 2548 , 2685 , 2598 , 2751 , 2697 , 2634 , 2450 , 2125 , 2020 , 2132 , 2345 , 2510 , 2625 , 2501 , 2245 , 2133 ,
   2247 , 2385 , 2444 , 2279 , 2206 , 2028 , 2113 , 2245 , 2213 , 2243 , 2190]

#observed number of confirmed infection cases in Iran over the period of 28 Aguest to 26 October
# ir4=[369911 , 371816 , 373570 , 375212 , 376894 , 378752 , 380746 , 382772 , 384666 , 386658 , 388810 , 391112 , 393425 ,
#  395488 , 397801 , 399940 , 402029 , 404648 , 407353 , 410334 , 413149 , 416198 , 419043 , 422140 , 425481 , 429193 ,
#   432798 , 436319 , 439882 , 443086 , 446448 , 449960 , 453637 , 457219 , 461044 , 464596 , 468119 , 471772 , 475674 ,
#    479825 , 483844 , 488236 , 492378 , 496253 , 500075 , 504281 , 508389 , 513219 , 517835 , 522387 , 526490 , 530380 ,
#     534631 , 539670 , 545286 , 550757 , 556891 , 562705 , 568896 , 574856]
# Number of new cases in Iran over the period of 28 Aguest to 26 October
ir4_newcases=[2115 , 1905 , 1754 , 1642 , 1682 , 1858 , 1994 , 2026 , 1894 , 1992 , 2152 , 2302 , 2313 , 2063 , 2313 , 2139 ,
 2089 , 2619 , 2705 , 2981 , 2815 , 3049 , 2845 , 3097 , 3341 , 3712 , 3605 , 3521 , 3563 , 3204 , 3362 , 3512 , 3677 , 3582 ,
  3825 , 3552 , 3523 , 3653 , 3902 , 4151 , 4019 , 4392 , 4142 , 3875 , 3822 , 4206 , 4108 , 4830 , 4616 , 4552 , 4103 , 3890 ,
   4251 , 5039 , 5616 , 5471 , 6134 , 5814 , 6191 , 5960]

#observed number of confirmed infection cases in Iran over the period of 26 October to 26 December 2020
# ir5=[581824 , 588648 , 596941 , 604952 , 612772 , 620491 , 628780 , 637712 , 646164 , 654936 , 663800 , 673250 , 682486 ,
#  692949 , 703288 , 715068 , 726585 , 738322 , 749525 , 762068 , 775121 , 788473 , 801894 , 815117 , 828377 , 841308 , 854361 ,
#   866821 , 880542 , 894385 , 908346 , 922397 , 935799 , 948749 , 962070 , 975951 , 989572 , 1003494 , 1016835 , 1028986 ,
#    1040547 , 1051374 , 1062397 , 1072620 , 1083023 , 1092407 , 1100818 , 1108269 , 1115770 , 1123474 , 1131077 , 1138530 ,
#     1145651 , 1152072 , 1158384 , 1164535 , 1170743 , 1177004 , 1183182 , 1189203]
# Number of new cases in Iran over the period of 26 October to 26 December 2020
ir5_newcases=[6968 , 6824 , 8293 , 8011 , 7820 , 7719 , 8289 , 8932 , 8452 , 8772 , 8864 , 9450 , 9236 , 10463 , 10339 , 11780 ,
 11517 , 11737 , 11203 , 12543 , 13053 , 13352 , 13421 , 13223 , 13260 , 12931 , 13053 , 12460 , 13721 , 13843 , 13961 , 14051 ,
  13402 , 12950 , 13321 , 13881 , 13621 , 13922 , 13341 , 12151 , 11561 , 10827 , 11023 , 10223 , 10403 , 9384 , 8411 , 7451 ,
   7501 , 7704 , 7603 , 7453 , 7121 , 6421 , 6312 , 6151 , 6208 , 6261 , 6178 , 6021]
#--------------------------------------

#Initialize Rt for each Simulation (infection rate patterns)
rr_ir1 = [5.95,2.25,3.2,0.03,0.86,4.15,0.75,0,0.3,2.48,1.3,1.3,1.92,1.2,0.55,0.58,1.31,1.2,0.25,1.8,0.33,1.3,2,1.6,2,1.24,
        1.78,1.15,0.85,1.3,0.92,0.8,0.87,0.8,1,0.76,0.7,0.84,0.68,0.6,1.65,0.65,0.7,0.98,0.8,0.9,1.1,0.85,0.8,0.83,0.77,1.2,0.7,
        0.585,1.35,0.95,1.1,0.65,1.35,0.92]

rr_ir2 = [2.8,1.15,0.25,1.915,1.966,1.15,1.68,0.72,1.213,0.94,0.687,1.81,0.69,1.93,0.71,1.5,0.5,1.07,1.81,0.83,1.38,1.07,0.81,
        0.45,1.52,0.687,0.716,1.49,1.25,1.83,0.488,1.284,1.52,1.022,1.04,1.43,0.49,0.3,1.1,0.58,1.03,0.78,1.44,1.12,1.3,0.94,1.03,
        1.18,0.91,1.02,1.03,0.6,1.13,1.31,0.7,1.1,1.08,1.1,0.82,0.96]

rr_ir3 = [3.3, 0.9, 1.2, 1.09, 0.92, 0.875, 1.13, 0.98, 1, 1.1, 0.39, 1.32, 1.16, 0.63, 1.3, 1.205, 0.854, 1.13, 0.775,
          0.64, 1.1, 1.35, 1.235, 0.973, 1.01, 0.8, 0.95, 0.91, 1.26, 1.23, 0.95, 0.99, 1.07, 0.9, 1.1, 0.99, 1.205,
          0.9, 0.94, 0.75, 0.55, 0.77, 1.23, 1.4, 1.25, 1.13, 0.66, 0.7, 0.83, 1.23, 1.28, 1.08, 0.81, 0.83, 0.72, 1.2,
          1.12, 1.02, 1.05, 0.87]

rr_ir4=[3.1,0.70,0.65,0.88,1.118,1.25,1.27,1,0.8,1.28,1.154,1.137,1.1,0.77,1.28,0.80,0.855,1.84,1.05,1.43,0.8,1.18,0.75,
        1.33,1.28,1.33,0.95,0.8,1.12,0.67,1.22,1.1,1.12,1,1.25,0.7,1,1.1,1.25,1.2,0.79,1.35,0.77,0.85,1.02,1.3,0.86,
        1.57,0.95,0.92,0.75,0.83,1.3,1.65,1.25,1.4,1,1.3]

rr_ir5=[3.705,0.964,1.665,0.87,0.95,1.01,1.2,1.205,0.79,1.18,1.11,1.17,0.94,1.34,0.955,1.456,0.92,1.087,0.9,1.332,1.11,
    1.072,1.03,0.945,1.004,0.922,1.033,0.837,1.33,1.035,1.037,1.032,0.848,0.883,1.082,1.17,0.95,1.038,0.883,0.809,0.79,0.795,
    1.014,0.801,1.042,0.66,0.68,0.68,1.01,1.124,0.9,1.02,0.79,0.66,0.98,0.895,1.079,1.014,0.92,0.91]
#--------------------------------------
# A bootstrap procedure for the point estimates and interval estimates:
random.seed(10) # Set Seed
rn.seed(10)     # Set Seed
outMA = []
newMA = []
gtotalA = []
nd=60
iteration=1  # number of iteration of Simulation

# Main Code
for m in range(iteration):
    simu=Simulation(nd=60,Rt=rr_ir1, muT=2.8,sizeV=1.3,limit=80000000,n0=205)    # ir1
    #simu=Simulation(nd=60,Rt=rr_ir2, muT=2.8,sizeV=1.3,limit=80000000,n0=1073)   # ir2
    #simu=Simulation(nd=60,Rt=rr_ir3, muT=2.8,sizeV=1.3,limit=80000000,n0=2489)   # ir3
    #simu=Simulation(nd=60,Rt=rr_ir4, muT=2.8,sizeV=1.3,limit=80000000,n0=2190)   # ir4
    #simu=Simulation(nd=60,Rt=rr_ir5, muT=2.8,sizeV=1.3,limit=80000000,n0=5960)   # ir5
    this=simu[0][:nd]
    outMA.append(this)
    thisnew=simu[1][:nd]
    newMA.append(thisnew)
    gtotalA=gtotalA+[simu[2]]
print(outMA)
print(newMA)
print(gtotalA)
plt.bar([i for i in range(1,61)],ir1_newcases)
plt.bar([i for i in range(1,61)],newMA[0],0.4)
plt.show()
