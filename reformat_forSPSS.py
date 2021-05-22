import pandas as pd
import collections
df = pd.read_excel('exportedforSPSS.xls')
data = {}
for i, row in df.iterrows():
    tup = (row['subjectID'], row['group'])
    if tup not in data:
        data[tup] = {}
    session = row['session'].lstrip('time')
    source = str(row['source'])
    detector = str(row['detector'])
    typ = row['type']
    cond = row['cond']
    item_name = 'Session' + session + '_' + 'Source' + source + '_' + 'D' + detector + '_' + typ + '_' + cond + 'Beta'
    data[tup][item_name] = [row['beta']]
data_out = {}
for tup in data:
    data_out[tup] = collections.OrderedDict(sorted(data[tup].items()))
df_out = []
for k, v in data_out.items():
    df_tmp = pd.DataFrame.from_dict(v)
    df_tmp.insert(0,'subjectId',[k[0]])
    df_tmp.insert(1,'group',[k[1]])
    df_out.append(df_tmp)
pd.concat(df_out).to_excel('result_rearranged.xls',index=False)