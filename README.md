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

