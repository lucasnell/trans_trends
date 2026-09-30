# `data`



Folder contents:

```
.
├── gadm41_ISL_0.geojson
├── myv_arth.csv
├── Myvatn_WSGUTM28.geojson
├── pitfall_locations.csv
├── raw
│   ├── infalls-08-17.csv
│   ├── infalls-13-19.csv
│   ├── metadata-infalls.csv
│   ├── metadata-pitfalls.csv
│   ├── pitfalls-08-17.csv
│   └── pitfalls-11-19.csv
└── README.md
```




File descriptions:

- `gadm41_ISL_0.geojson`: outline of Iceland
- `myv_arth.csv`: cleaned dataset; see below for column descriptions
- `Myvatn_WSGUTM28.geojson`: outline of Lake Myvatn, Iceland
- `pitfall_locations.csv`: UTM (zone 28, WGS84) locations of all pitfall traps 
  by site and distance from the lake
- `raw/infalls-08-17.csv`: raw data for infall traps catching midge emergences
  (2008-06-03 to 2017-07-29);
  many rows are duplicated with `raw/infalls-13-19.csv`
- `raw/infalls-13-19.csv`: raw data for infall traps catching midge emergences
  (2013-06-18 to 2019-08-08);
  many rows are duplicated with `raw/infalls-08-17.csv`
- `raw/metadata-infalls.csv`: column names and descriptions for
  `raw/infalls*.csv` files
- `raw/metadata-pitfalls.csv`: column names and descriptions for
  `raw/pitfalls*.csv` files
- `raw/pitfalls-08-17.csv`: raw data for pitfall traps catching ground-dwelling 
  predatory arthropods (2008-06-25 to 2017-07-29);
  many rows are duplicated with `raw/pitfalls-11-19.csv`
- `raw/pitfalls-11-19.csv`: raw data for pitfall traps catching ground-dwelling 
  predatory arthropods (2011-07-17 to 2019-08-08);
  many rows are duplicated with `raw/pitfalls-08-17.csv`


Cleaned data column descriptions:

- `trans`: Transect location at Myvatn; unique values are
  `"btl"` (site Fellshóll), `"hag"` (Haganes), `"kal"` (Kálfaströnd),
  `"non"` (Nóntangi), and `"vin"` (Vindbelgur)
- `dist`: Distance from the lake; unique values are 5, 50, 150, 500
- `year`: Year sampled (2008--2019)
- `midges`: Number of total midges in sample on the sampled date
- `small_midges`: Number of small midges in sample on the sampled date
- `big_midges`: Number of large midges in sample on the sampled date
- `taxon`: Taxon that count refers to; unique values are
  `"acar"` (mites), `"cara"` (ground beetles), `"stap"` (rove beetles), 
  `"opil"` (harvestmen), `"gnap"` (ground spiders), `"sheet"` (sheet weavers), 
  and `"lyco"` (wolf spider).
  Mites were not included in our analyses because they are not particularly
  mobile, and many are not predators.
- `count`: Number of individual arthropods of the specified taxon caught 
  in the sample
- `first_day`: First day that (non-midge) arthropod sample was out
- `last_day`: Last day that (non-midge) arthropod sample was out
- `season_days`: Number of days that (non-midge) arthropod sample was out
- `midge_first_day`: First day that midge sample was out
- `midge_last_day`: Last day that midge sample was out
- `midge_days`: Number of days that midge sample was out
