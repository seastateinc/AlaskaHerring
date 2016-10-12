## —————————————————————————————————————————————————————————————————————————— ##
##                  DESIGN MATRIX FOR PARAMETER CONTROLS                      ##
##  Prior descriptions   Parameter values                                     ##                                                                            ##
##  -0 uniform           (0,0)                                                ##
##  -1 normal            (p1=mu,p2=sig)                                       ##
##  -2 lognormal         (p1=log(mu),p2=sig)                                  ##
##  -3 beta              (p1=alpha,p2=beta)                                   ##
##  -4 gamma             (p1=alpha,p2=beta)                                   ##
## —————————————————————————————————————————————————————————————————————————— ##
##  init  lower   upper    est  prior
## value  bound   bound    phz   type     p1    p2   # PARAMETER              ##
## —————————————————————————————————————————————————————————————————————————— ##
   -1.05   -6.79    1.00      2      1  -0.69  0.05   # log_natural_mortality
    7.50   -6.00   12.00      1      0      0     0   # log_rinit
    7.50   -6.00   12.00      1      0      0     0   # log_rbar
    7.50   -6.00   12.00      2      0      0     0   # log_ro
    1.50    0.00   12.00      2      0      0     0   # log_reck
    0.91   -9.00    1.00     -2      0      0     0   # log_sigma_r
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##                CONTROLS FOR TIME-VARYING MATURITY                          ##
## —————————————————————————————————————————————————————————————————————————— ##
## nMatBlocks
  1
## a50    a95     phz   terminalBlockYear
   4.5    7.0      -2           2015
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##            CONTROLS FOR TIME-VARYING LN(NATURAL MORTALITY DEVS)            ##
## KEY:
##  Type: 1 = constant M
##  Type: 2 = interpolated using cubic spline.
## —————————————————————————————————————————————————————————————————————————— ##
## Type
  2
## Phase for estimation if nMortBlocks > 1
  2
## nMortBlocks, or Nodes in the case of cubic spline interpolation
  4
## The terminal year of each block
 1977 1998 2002 2007
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##                    CONTROLS FOR SELECTIVITY PARAMETERS                     ##
## —————————————————————————————————————————————————————————————————————————— ##
## - Each selectivity block can have different functional forms based on selType
## - LEGEND:
##   - SelType = 1, logistic selectivity, 2 parameters.
##   - SelType = 2, logistic with 50% & 95% parameters
##  nSelexblocks
    3
## —————————————————————————————————————————————————————————————————————————— ##
##  Gear  Sel     sel   sel   age   year  phz    | start end
##  Index Type    mu    sd    nodes nodes mirror | block block
## —————————————————————————————————————————————————————————————————————————— ##
    1     1       3.0   0.5   0     0     -1        1971  1980
    1     1       5.0   0.3   0     0      2        1981  2000
    1     1       5.0   0.3   0     0      2        2001  2015
## —————————————————————————————————————————————————————————————————————————— ##


## —————————————————————————————————————————————————————————————————————————— ##
##                        OTHER MISCELLANEOUS CONTROLS                        ##
## —————————————————————————————————————————————————————————————————————————— ##
## number of controls to read in.
5
## Value    # - Description
0.90718     # - Catch Scaler (convert from short tons to metric tons)
1           # - Condition on Catch = 0, Condition of Ft = 1
25000       # - harvest threshold
0.10        # - target harvest rate
20000       # - threshold denominator
## EOF
999