from datasets import load_dataset

### GreenHyperSpectra: unlabeled ###
ds_un = load_dataset("Avatarr05/GreenHyperSpectra", "unlabeled")
GreenHyperSpectra = ds_un['train'].to_pandas().drop(['Unnamed: 0'], axis=1)

display(GreenHyperSpectra.head())

### Labeled data: labeled_all ###
ds = load_dataset("Avatarr05/GreenHyperSpectra", "labeled_all")
df = ds['train'].to_pandas().drop(['Unnamed: 0'], axis=1)

display(df.head())

annotated_ds_train = load_dataset("Avatarr05/GreenHyperSpectra", 'labeled_splits', split="train")
annotated_ds_train = annotated_ds_train['train'].to_pandas().drop(['Unnamed: 0'], axis=1)

annotated_ds_test = load_dataset("Avatarr05/GreenHyperSpectra", 'labeled_splits', split="test")
annotated_ds_test = annotated_ds_test['train'].to_pandas().drop(['Unnamed: 0'], axis=1)

display(annotated_ds_train.head())
display(annotated_ds_test.head())