# NNLOJET combine settings:

Overriding the default stettings of the combine procedure is done within the `[Options]` section of the `combine.ini` file. The outlier rejection is controlled by two numbers: 
1. The theshold (something like a deviation from the mean in units of sigma; but computed in terms of the inter-quantile distances). *Lowering* this value will flag *more* points as outliers (default = `4`).
2. The highest fraction of points that are allowed to be removed. If step 1 flags too many points in a bin, adjust it dynamically to always stay below this fraction. This is a safety feature to avoid rejecting too many points (default = `0.05 = 5%`).

## default:
```
trim = (4, 0.05)
```

##  more aggresive outlier identification 
* reduce the threshold => flag more outliers
```
# instead of ~4 sigma, start flagging at ~3 sigma
trim = (3, 0.05)
```

## relax constraint on the maximal fraction
* increase max fraction => allow more outliers to be removed
```
# instead of 5%, allow up to 8% of points to be dropped 
trim = (4, 0.08)
```

## combine both settings
```
trim = (3, 0.08)
```
