# MMSB2022_HW5
110.2 BEBI5009_生物系統模擬


# Stochastic Simulation Algorithm（Gillespie Algorithm）

- bistable：http://www.best.org.tw/upload/downloads_upload/2020_12_BEST%E5%AD%A3%E8%A8%8A-yswu.pdf p.16, 1.2 基因調控中的雙穩態系統.

- [Result]()


# Agent-based Modeling

## Question 1：

Please simulate a combination of 
- beta_max = 0.05, 0.04, 0.03, 0.02, 0.01 (representing personal hygiene: e.g. Wearing a mask) and
- isolated = 0.0, 0.5, 0.7, 0.8, 0.9. (Representing social distancing and/or lockdown measures) 
And thus there will be 25 simulations in total. 

Plot the number of infected individuals over 5000 steps with 5 beta_max parametegrs, each isolated parameter on a different plot.

(That is, plot #1 is isolated = 0.0 with 5 time series: beta_max = 0.05, 0.04, 0.03, 0.02, 0.01, plot #2 is isolated = 0.5 with 5 time series: beta_max = 0.05, 0.04, 0.03, 0.02, 0.01, and so on.)

### SIR Model, isolated = 0.0
![SIR Model, isolated = 0.0](./PNG/SIR_Model_isolated_000.png)
```
Q2, maxima (peaks), isolated = 0.0：

- βmax = 0.05  ==> 1000
- βmax = 0.04  ==> 1000
- βmax = 0.03  ==> 1000
- βmax = 0.02  ==> 1000
- βmax = 0.01  ==> 1000
```

### SIR Model, isolated = 0.5
![SIR Model, isolated = 0.5](./PNG/SIR_Model_isolated_050.png)
```
Q2, maxima (peaks), isolated = 0.5：

- βmax = 0.05  ==> 983
- βmax = 0.04  ==> 983
- βmax = 0.03  ==> 983
- βmax = 0.02  ==> 983
- βmax = 0.01  ==> 983
```

### SIR Model, isolated = 0.7
![SIR Model, isolated = 0.7](./PNG/SIR_Model_isolated_070.png)
```
Q2, maxima (peaks), isolated = 0.7：

- βmax = 0.05  ==> 857
- βmax = 0.04  ==> 878
- βmax = 0.03  ==> 873
- βmax = 0.02  ==> 876
- βmax = 0.01  ==> 882
```

### SIR Model, isolated = 0.8
![SIR Model, isolated = 0.8](./PNG/SIR_Model_isolated_080.png)
```
Q2, maxima (peaks), isolated = 0.8：

- βmax = 0.05  ==> 754
- βmax = 0.04  ==> 769
- βmax = 0.03  ==> 757
- βmax = 0.02  ==> 771
- βmax = 0.01  ==> 776
```

### SIR Model, isolated = 0.9
![SIR Model, isolated = 0.9](./PNG/SIR_Model_isolated_090.png)
```
Q2, maxima (peaks), isolated = 0.9：

- βmax = 0.05  ==> 376
- βmax = 0.04  ==> 418
- βmax = 0.03  ==> 406
- βmax = 0.02  ==> 398
- βmax = 0.01  ==> 415
```

## Question 2：

- List the maxima (peaks) of infected individuals in the 25 simulations. Which parameter set is the most effective in "flattening the curve" (having the lowest peak infected individuals)?

        isolated 0.9, βmax = 0.05
 
- Compared to decreasing the beta_max parameter, is increasing the "isolated" parameter more effective in flattening the curve?

        YES
