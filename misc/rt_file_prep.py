#!/usr/bin/env python3

import re
import argparse
import pandas as pd
import dask.dataframe as dd
'''
This converts ionbot output to input for deeplc - for retention time predictions
'''
def fixed_to_mod(seq,aa="C",mod_name="Carbamidomethyl"):
    ret_str = []
    matches = re.finditer(aa,seq)
    for match in matches:
        ret_str.append("{}|{}".format(match.start(0)+1,mod_name))
    return "|".join(ret_str)

def overall_mod(dd):
    dd=dd[dd['ri_126.1277']>0].fillna('')
    dd['fixed']=dd['modifications'].apply(fixed_to_mod,meta=pd.Series(dtype='str',name='fixed'))
    dd['j']=dd['fixed'].apply(lambda x: '' if x=='' else '|',meta=pd.Series(dtype='str',name='j'))
    dd['j2']=dd['unexpected_modification'].apply(lambda x: '' if x=='' else '|',meta=pd.Series(dtype='str',name='j2'))
    dd['modifications']=dd['modifications']+dd['j']+dd['fixed']+dd['j2']+dd['unexpected_modification']
    return(dd.compute(scheduler='processes',num_workers=30))
    
def main():
    parser = argparse.ArgumentParser(description='Ionbot output analysis')
    parser.add_argument('--ib', help='Directory of ionbot output files', required=True)
    parser.add_argument('--out', help='Output prefix', required=True)
    args = vars(parser.parse_args())

    output_all=args['out']+'.csv'
    output_calibration=args['out']+'_calibration.csv'
    ionbotdf=import_data(dd.read_csv(f"{args['ib'].strip('/')}/*.mgf.ionbot.csv"))
    ionbotdf=ionbotdf.sort_values('percolator_psm_score',ascending=False)
    ionbotdf=ionbotdf[['matched_peptide','modifications','rt']].drop_duplicates(subset=['matched_peptide','modifications']).rename({'matched_peptide':'seq','rt':'tr'},axis=1)
    ionbotdf.to_csv(output_all,index=False)
    ionbotdf.head(1000).to_csv(output_calibration,index=False)
    

if __name__ == "__main__":
    main()
