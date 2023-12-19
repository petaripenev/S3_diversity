import os
from collections import defaultdict
import pandas as pd


df_sample = pd.read_csv(config["samples"], sep="\t")
df_sample = df_sample.fillna("")
df_sample[df_sample.columns] = df_sample.apply(lambda x: x.str.strip())

name_to_path = defaultdict(lambda: defaultdict(dict))

for header_col in df_sample:
    if header_col.endswith("_path"):
        name_header = "_".join([header_col.split("_path")[0], "name"])

        df_sample[name_header] = df_sample[header_col].apply(os.path.basename)

sample_to_info = df_sample.set_index("sample_name").to_dict("index")
