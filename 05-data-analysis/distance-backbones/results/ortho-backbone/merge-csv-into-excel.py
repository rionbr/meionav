
# code to merge files into excel, because excel can't read columns as strings.
import pandas as pd

dfHS = pd.read_csv('net-ortho-backbone-HS.csv.gz', index_col=0, usecols=[0, 2, 3])
dfMM = pd.read_csv('net-ortho-backbone-MM.csv.gz', index_col=0, usecols=[0, 2, 3])
dfDM = pd.read_csv('net-ortho-backbone-DM.csv.gz', index_col=0, usecols=[0, 2, 3])

dfHS.index.name = 'id'
dfHS.rename(columns={'label': 'gene'}, inplace=True)
dfHS.reset_index()

dfMM.index.name = 'id'
dfMM.rename(columns={'label': 'gene'}, inplace=True)

dfDM.index.name = 'id'
dfDM.rename(columns={'label': 'gene'}, inplace=True)


with pd.ExcelWriter('net-ortho-backbone-HS-MM-DM.xlsx') as writer:
    dfHS.to_excel(writer, sheet_name='Human')
    dfMM.to_excel(writer, sheet_name='Mouse')
    dfDM.to_excel(writer, sheet_name='Fruit fly')
