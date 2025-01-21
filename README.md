# Cell line annotation

This parses cellosaurus data into a table.

## Requirements

### Python
- pyarrow

### R
- io
- arrow

## Procedure

1. Download the data
```
( cd data && ./get.sh )
```

2. Parse the data
```
( cd feather && python parse.py )
```

3. Pre-process the data
```
cd rds
Rscript preprocess.R
```

Then, the output `rds` files will be in the `rds` directory.

4. Filter and obtain human cell lines by running filter.R
```
cd filter
Rscript filter.R
```

5. Fetch GEO information from embryonic stem cells
``` 
( cd geo && ./run_ESC.sh )

```

6. Download cancer cell line data from CCLE and GDSC
```
( cd cell-lines && ./get.sh )
```
This will create folders ccle and gdsc and download data into them.


