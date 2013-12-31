# -*- coding: utf-8 -*-

import pandas as pd

cosmic_df = pd.read_csv('CosmicV67-Recep.txt', sep='\t')
cosmic_driver_df = cosmic_df[cosmic_df['mType'] == 'Driver']
uniprots_to_precalculate = set(cosmic_driver_df['uniprot'])
