# Question about overlapping 2012 pitfall dates at Mývatn

*Note: dates are written year/month/day (e.g., 2012/06/05), or month/day when the year is clear from context (e.g., 06/05).*

While cleaning the Mývatn pitfall data, I found 2012 samples whose
deployment dates overlap with another sample from the same trap. Could you
check them against the field sheets, or tell me who could?

## What I'm seeing

All of these rows are in the 2008–2017 pitfall file
(`pitfalls-08-17.csv`, formerly `myvatn_pitfalls.csv`). 54 samples fall into
23 groups in which one sample's set date comes before the previous sample's
collection date at the same plot. The overlaps follow three patterns:

1. **June, all five sites:** a sample collected 06/22 or 06/23, followed by
   one set 06/19 (a 3–4 day overlap).
2. **July, Fellshóll (btl), Haganes (hag), and Nóntangi (non):** two samples
   with the same 07/17 set date, one collected 07/31 and the other 08/09.
3. **Nóntangi 150 m and Vindbelgur (vin) 5, 50, and 150 m:** a 44-day sample
   (06/26 to 08/09) that spans both the 06/19–07/03 and 07/17–07/31 samples.

Counts differ within every group, so these look like separate catches rather
than duplicated rows. Two other things stood out:

- Most of the 2012 half counts (0.5, 1.5, 6.5, and so on) are in the
  06/19–07/03 samples.
- In the July pairs, the sample ending 08/09 usually has fewer animals than
  the one ending 07/31, despite being out 9 days longer.

## Why it matters

We model catch per trap-day, so I need the correct number of days each trap
was open. If a set date was entered wrongly, every sample stays and only its
days change. If traps really did overlap, the days are right as they stand.
I'd rather not guess which.

## Questions

1. In June, should the 06/19 set dates be 06/22 or 06/23 (the day the
   previous sample was collected)?
2. In July, should the second 07/17 set date be 07/31?
3. Are the 06/26 to 08/09 dates at Nóntangi 150 m and Vindbelgur 5, 50, and
   150 m correct? A 44-day deployment is about three times the usual length.
4. Are the half counts in the 06/19–07/03 samples averages of two traps (or of
   two counts)? If they're averages of two traps, I'd need to either sum the
   counts or double the trap-days.
5. Are there other 2012 trap-level notes (lost, flooded, or extra traps) I
   should know about?

Thanks for your help.

## Appendix: all overlapping samples

Dates are 2012 set date–collection date (month/day). Days = collection date −
set date. Columns are raw-file column names; "aran + liny" is
`aran_other + liny` (sheet weavers).

| Plot | Dates | Days | cara | stap | opil | gnap | aran + liny | lyco_total |
|---|---|---|---|---|---|---|---|---|
| btl 5 m | 07/17–07/31 | 14 | 1 | 0 | 1 | 0 | 6 | 0 |
|  | 07/17–08/09 | 23 | 1 | 0 | 0 | 0 | 9 | 1 |
| btl 50 m | 06/05–06/22 | 17 | 7 | 5 | 6 | 18 | 18 | 19 |
|  | 06/19–07/03 | 14 | 2 | 0 | 0 | 0 | 9 | 0 |
| btl 50 m | 07/17–07/31 | 14 | 1 | 0 | 0 | 1 | 1 | 2 |
|  | 07/17–08/09 | 23 | 1 | 0 | 3 | 0 | 0 | 4 |
| btl 150 m | 06/05–06/22 | 17 | 10 | 6 | 14 | 0 | 47 | 31 |
|  | 06/19–07/03 | 14 | 16 | 3 | 4 | 5 | 14 | 10 |
| btl 150 m | 07/17–07/31 | 14 | 1 | 0 | 2 | 0 | 7 | 0 |
|  | 07/17–08/09 | 23 | 0 | 0 | 1 | 0 | 0 | 0 |
| hag 5 m | 06/05–06/22 | 17 | 3 | 0 | 8 | 9 | 10 | 22 |
|  | 06/19–07/03 | 14 | 3.5 | 3 | 1.5 | 0 | 4.5 | 19.5 |
| hag 5 m | 07/17–07/31 | 14 | 7 | 9 | 2 | 2 | 48 | 22 |
|  | 07/17–08/09 | 23 | 3 | 0 | 1 | 2 | 4 | 11 |
| hag 50 m | 06/05–06/22 | 17 | 0 | 1 | 12 | 0 | 9 | 12 |
|  | 06/19–07/03 | 14 | 1.5 | 0.5 | 1 | 1 | 4 | 15.5 |
| hag 50 m | 07/17–07/31 | 14 | 1 | 2 | 0 | 5 | 2 | 11 |
|  | 07/17–08/09 | 23 | 0 | 0 | 0 | 1 | 0 | 9 |
| hag 150 m | 06/05–06/22 | 17 | 0 | 0 | 2 | 8 | 8 | 17 |
|  | 06/19–07/03 | 14 | 0.5 | 0 | 0.5 | 1 | 6.5 | 7 |
| hag 150 m | 07/17–07/31 | 14 | 0 | 0 | 0 | 1 | 1 | 7 |
|  | 07/17–08/09 | 23 | 0 | 0 | 0 | 1 | 0 | 2 |
| hag 500 m | 06/05–06/22 | 17 | 0 | 1 | 2 | 7 | 15 | 22 |
|  | 06/19–07/03 | 14 | 1 | 0 | 1 | 1 | 6.5 | 18.5 |
| hag 500 m | 07/17–07/31 | 14 | 0 | 0 | 3 | 1 | 1 | 5 |
|  | 07/17–08/09 | 23 | 0 | 0 | 0 | 0 | 0 | 2 |
| kal 5 m | 06/05–06/22 | 17 | 2 | 0 | 8 | 4 | 23 | 10 |
|  | 06/19–07/03 | 14 | 1 | 0 | 4 | 2 | 3.5 | 2.5 |
| kal 50 m | 06/05–06/22 | 17 | 2 | 1 | 1 | 10 | 5 | 17 |
|  | 06/19–07/03 | 14 | 0 | 0 | 0.5 | 4 | 4.5 | 20.5 |
| kal 150 m | 06/05–06/22 | 17 | 0 | 0 | 1 | 9 | 7 | 15 |
|  | 06/19–07/03 | 14 | 0 | 0 | 0 | 4 | 7.5 | 6.5 |
| non 5 m | 06/05–06/23 | 18 | 13 | 18 | 3 | 1 | 23 | 0 |
|  | 06/19–07/03 | 14 | 47 | 8.5 | 2.5 | 0 | 2.5 | 0 |
| non 5 m | 07/17–07/31 | 14 | 21 | 3 | 0 | 2 | 0 | 0 |
|  | 07/17–08/09 | 23 | 8 | 3 | 0 | 5 | 0 | 0 |
| non 50 m | 06/05–06/23 | 18 | 1 | 1 | 0 | 0 | 3 | 4 |
|  | 06/19–07/03 | 14 | 1 | 1 | 0 | 1.5 | 8 | 2.5 |
| non 150 m | 06/05–06/23 | 18 | 0 | 1 | 3 | 0 | 12 | 15 |
|  | 06/19–07/03 | 14 | 0 | 0 | 0 | 1 | 3.5 | 1.5 |
|  | 06/26–08/09 | 44 | 0 | 0 | 1 | 7 | 3 | 5 |
|  | 07/17–07/31 | 14 | 13 | 0 | 1 | 1 | 2 | 2 |
| vin 5 m | 06/05–06/23 | 18 | 2 | 1 | 4 | 19 | 40 | 80 |
|  | 06/19–07/03 | 14 | 3 | 5.5 | 8 | 3.5 | 9 | 10 |
|  | 06/26–08/09 | 44 | 0 | 7 | 6 | 2 | 2 | 3 |
|  | 07/17–07/31 | 14 | 1 | 18 | 9 | 1 | 5 | 2 |
| vin 50 m | 06/05–06/23 | 18 | 1 | 0 | 5 | 72 | 31 | 90 |
|  | 06/19–07/03 | 14 | 4 | 0.5 | 1 | 6.5 | 5 | 39.5 |
|  | 06/26–08/09 | 44 | 0 | 1 | 1 | 1 | 1 | 10 |
|  | 07/17–07/31 | 14 | 0 | 9 | 0 | 1 | 2 | 22 |
| vin 150 m | 06/05–06/23 | 18 | 3 | 1 | 1 | 16 | 15 | 32 |
|  | 06/19–07/03 | 14 | 1.5 | 0.5 | 0.5 | 2.5 | 8.5 | 16.5 |
|  | 06/26–08/09 | 44 | 0 | 0 | 0 | 1 | 4 | 6 |
|  | 07/17–07/31 | 14 | 0 | 0 | 0 | 0 | 18 | 8 |
