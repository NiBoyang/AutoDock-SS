def create_subsample(lismol):
    
    import gzip
    import random
    from rdkit import Chem
    
    def slice_per(source, step):
        return [source[i::step] for i in range(step)]
    
    with gzip.open(lismol) as temp:
        with Chem.ForwardSDMolSupplier(temp) as gzsuppl:
            mols = [x for x in gzsuppl if x is not None]
        names = [mol.GetProp("_Name") for mol in mols]
        filtered_name = list(dict.fromkeys(names))
        sliced_lis = slice_per(filtered_name, len(filtered_name)//10)
        random.seed(38420)
        return set(filtered_name) - set([random.choice(k) for k in sliced_lis]), filtered_name