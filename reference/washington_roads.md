# Washington Road Crashes

Crashes on Washington primary roads from 2016, 2017, and 2018. Data
acquired from Washington Department of Transportation through the
Highway Safety Information System (HSIS).

## Usage

``` r
washington_roads
```

## Format

Data frame compiled from roadway, traffic, and police-reported crash
data that has 1,501 rows and 13 columns:

- ID:

  Anonymized road ID. Factor.

- Year:

  Year. Integer.

- AADT:

  Annual Average Daily Traffic (AADT). Double.

- Length:

  Segment length in miles. Double.

- Total_crashes:

  Total crashes. Integer.

- lnaadt:

  Natural logarithm of AADT. Double.

- lnlength:

  Natural logarithm of length in miles. Double.

- speed50:

  Indicator of whether the speed limit is 50 mph or greater. Binary.

- ShouldWidth04:

  Indicator of whether the shoulder is 4 feet or wider. Binary.

- Fatal_crashes:

  Total number of non-intersection fatal crashes for the road segment

- Injury_crashes:

  Total number of non-intersection Injury crashes for the road segment

- Animal:

  Total number of non-intersection animal-related crashes for the road
  segment

- Rollover:

  Total number of non-intersection rollover crashes for the road segment

## Source

<https://highways.dot.gov/research/safety/hsis>
