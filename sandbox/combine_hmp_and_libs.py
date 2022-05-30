import pandas as pd
import feather

output = "outputs/hash_tables/all_unnormalized_abund_hashes_wide.feather"

libs = pd.read_feather("outputs/hash_tables/libs_unnormalized_abund_hashes_wide.feather")
hmp =  pd.read_feather("outputs/hash_tables/hmp_unnormalized_abund_hashes_wide.feather")
all = pd.concat([libs, hmp], axis = 0, ignore_index = True, sort=False)
all = all.fillna(0)
all.to_feather(str(output))