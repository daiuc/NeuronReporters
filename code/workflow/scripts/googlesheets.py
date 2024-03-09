import pandas as pd


sheet_id = '1dOU-b__H9pAhNUjS4vXvkZhgo1tP_8QMpxWTRiFXi-Y'

sheet_name = 'diffbind'
url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}'

sheet = pd.read_csv(url)
sheet.to_csv('config/diffBind_samplesheet.csv', index=False)

