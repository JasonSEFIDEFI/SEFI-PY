# SEFI:ELEMENTS working research package

Start with [MANUSCRIPT.md](MANUSCRIPT.md) and [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md).

This is a local exploratory reconstruction of the lost original element mapping. It includes real AME2016/AME2020/NUBASE2020 data, 52 mass and 52 lifetime method/split results plus 40 earlier-evaluation mass tests, held-out predictions, cluster bootstrap uncertainty, interval coverage, and tests. No new element is established.

Reproduce with Python 3.14.7 and requirements.lock.txt:

```
python download_data.py
python analyze.py
python historical.py
python uncertainty.py
python plot_paper.py
python -m pytest -q -p no:cacheprovider test_analysis.py
python build_paper.py
```

Download only accesses the three public IAEA URLs. Existing raw data are reused. Analysis overwrites generated result files, so copy the package if preserving a prior run. Geometric and matched conventional kernel predictions coincide. All reported target uncertainty intervals are empirical and fail under some extrapolation splits; read the manuscript before interpreting them.
