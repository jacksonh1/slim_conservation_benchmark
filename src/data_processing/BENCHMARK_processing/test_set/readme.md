calculate pairwise scores with specified parameters and associated keys



inputs:
- input table (csv file). contains a reference index for each entry
- parameters json file (created by the user):

```json
{
    "column name": {
        "scoring function": "function name",
        "parameters": {
            "parameter name": "parameter value"
        },
        "output folder": "folder name"
    }
}
```
- import table
- get the list info files
- for each column name in the parameters file:
    - execute with multiprocessing:
        - for each info file
            - calculate the score using the "scoring function" and the "parameters"
            - export the score results to a folder
                - file={output folder}/{reference index}_{column_name_no_spaces}.json
        - save the score results file name in a multiprocessing results file (reference index, score results file name)
    - save a copy of the parameters file in the "output folder"
    - import the multiprocessing results file
    - use the reference index from the table and results file to add a new column to the input table with the score results file name
- save the input table with the new columns






setting up folders and reindexing table file
calculating alignment scores
calculating pairwise matrices
calculating embedding pairwise matrices
13
no idr
no idr
no idr
no idr
calculating kmer scores from pairwise matrices
done
total time elapsed: 7.261046449343364 minutes (435.6627869606018 seconds)




