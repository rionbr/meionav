# Merge results from gene screening

Code in this folder is used to merge data from the selected core meiotic genes to the data from the experimental screening of DM genes.

Files `data/core_DM_screened.csv` and `data/core_DM_control.csv` contain the results from the experimental screening.

Instructions:
1. Run `merge_screening_data.py` to merge the core DM genes with the screened data.

Other files:
- `plot_screened.py` generates the plot for the experimental screening data.


## Set of genes that were screened based on a previous pipeline but do not appear in the new pipeline

```python
['FBgn0036557', 'FBgn0038960', 'FBgn0037880', 'FBgn0035707', 'FBgn0037440', 'FBgn0261610', 'FBgn0038053', 'FBgn0037359', 'FBgn0032471', 'FBgn0023512', 'FBgn0038486', 'FBgn0037009', 'FBgn0036536', 'FBgn0040650', 'FBgn0036316', 'FBgn0015031', 'FBgn0037008', 'FBgn0036354', 'FBgn0000303', 'FBgn0025186', 'FBgn0033644', 'FBgn0050035']
 ```