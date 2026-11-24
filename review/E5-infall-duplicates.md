# Question about duplicated infall samples at Mývatn (2008 and 2012)

*Note: dates are written year/month/day (e.g., 2012/06/05).*

While cleaning the Mývatn infall (midge) data, I found four samples that
appear twice in the 2008–2017 infall file (`infalls-08-17.csv`, formerly
`myvatn_infalls.csv`), with the same transect, distance, set date, and
collection date but different values. Could you check them against the
field sheets, or tell me who could?

None of these dates is in the 2013–2019 file, so the duplicates are in the
original data, not created by combining the two files. Earlier curation notes
in the file's `comments.2` column already flag all four pairs; I've quoted
them below.

## Why it matters

We model midge catch per trap-day as a predictor. At the moment both copies
of each pair are kept, so their midges are added together and their days
counted twice. Before choosing which row to keep, I'd like to know which one
is correct.

## The four pairs

Line numbers are lines in `infalls-08-17.csv`, counting the header as
line 1. "Big" and "small" are `bgch` and `smch` (raw counts, before dividing
by the subsample fraction `fract_bgch` or `fract_smch`).

### 1. Haganes 5 m, 2012/06/05 to 2012/06/19

| Line | Days | Big | Small | `fract_smch` | `simu` |
|---|---|---|---|---|---|
| 863 | 14 | 4 | 291 | 1 | NA |
| 864 | 14 | 7 | 468 | 1 | 5 |

Curation note: *"There are two entries with the same transect, distance, and
dates. One may be BTL 50 from the next sampling event, which is missing. No
way to tell."*

The file does fit that note. Fellshóll 50 m has no 2012/06/19 to 2012/07/03
sample (it has 5 samples in 2012, where most plots have 6), and Haganes 5 m
has 7.

### 2. Nóntangi 50 m, 2012/07/17 to 2012/07/31

| Line | Days | Big | Small | `fract_smch` |
|---|---|---|---|---|
| 942 | 14 | 0 | 360 | 1 |
| 944 | 14 | 0 | 601 | 0.229 |

Curation note: *"One of these is possibly kal 5. No way to tell."*

This also fits: Kálfaströnd 5 m has no 2012/07/17 to 2012/07/31 sample
(5 samples in 2012), and Nóntangi 50 m has 7.

### 3. Nóntangi 150 m, 2008/06/04 to 2008/06/10

| Line | `daysout` | Big | Small | `comments` |
|---|---|---|---|---|
| 26 | 6 | NA | NA | "Missing data??" |
| 27 | 7 | 26 | 173 | |

According to the curation note, line 27's collection date was changed from
2008/06/11 to 2008/06/10 because of the "missing data" entry. Its `daysout`
is still 7, while the dates now give 6. The next sample at this plot is set
on 2008/06/10.

### 4. Vindbelgur 5 m, 2008/05/30 to 2008/06/09

| Line | Days | Big | Small | `comments` |
|---|---|---|---|---|
| 17 | 10 | 1,842.5 | 28,482.4 | |
| 18 | 10 | NA | NA | "may have fallen (rubber band failure)??" |

According to the curation note, line 18 originally had zeros, which were
changed to NA because the other entry has values.

## Questions

1. Haganes 5 m (lines 863–864): is one of these actually Fellshóll 50 m,
   2012/06/19 to 2012/07/03? If so, which one?
2. Nóntangi 50 m (lines 942 and 944): is one of these actually Kálfaströnd
   5 m, 2012/07/17 to 2012/07/31? If so, which one? Is the 0.229 subsample
   fraction on line 944 correct?
3. Nóntangi 150 m (lines 26–27): is the collection date for line 27 2008/06/10
   or 2008/06/11? Can line 26 be dropped as a placeholder for a missing
   sample?
4. Vindbelgur 5 m (lines 17–18): was line 18 the same trap (fallen, so no
   usable sample), or a second trap? If it was the same trap, I'll drop
   line 18 and keep line 17.
5. More generally: was there ever more than one infall trap per plot? That
   would change how I treat any remaining duplicates.

Thanks for your help.
