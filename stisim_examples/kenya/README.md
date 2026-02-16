# Kenya HIV Model Data

Location-specific data for HIV modeling in Kenya.

## Status

⚠️ **Data files not yet populated**

This location is a placeholder for future Kenya-specific HIV modeling data.

## Required Data Files

To use this location, you need to provide the following files in this directory:

### STI-Specific Data (in stisim_examples/kenya/)
- `init_prev_hiv.csv` - Initial HIV prevalence by risk group/sex/sw
- `condom_use.csv` - Condom use by partnership type over time
- `art_coverage.csv` - ART coverage over time
- `vmmc_coverage.csv` - VMMC coverage over time

### Demographic Data (in stisim/data/files/)
- `kenya_asfr.csv` - Age-specific fertility rates
- `kenya_deaths.csv` - Death rates by age/sex/year
- `kenya_age_YYYY.csv` - Age distribution for year YYYY

## File Format References

See `stisim_examples/zimbabwe/README.md` and `stisim_examples/demo/README.md` for data file format specifications.

## External Data

If you have Kenya data in an external repository, you can use it directly:

```python
import hivsim as hs

sim = hs.Sim(
    location='kenya',
    datafolder='/path/to/kenya/data',
    n_agents=5000,
    dur=40
)
sim.run()
```

## Contributing

To add Kenya data to stisim_examples:

1. Obtain or calibrate HIV epidemic data for Kenya
2. Format data files according to specifications in zimbabwe/demo READMEs
3. Copy data files to appropriate locations
4. Create `kenya/sim.py` configuration helpers
5. Test with `stx.Sim(demographics='kenya', diseases='hiv')`

## Data Sources

Potential data sources for Kenya HIV modeling:
- UNAIDS Kenya Country Profile
- Kenya Demographic and Health Surveys (KDHS)
- UN World Population Prospects
- Kenya National AIDS Control Council (NACC) reports
