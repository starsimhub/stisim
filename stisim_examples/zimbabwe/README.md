# Zimbabwe HIV Model Data

Location-specific data for HIV modeling in Zimbabwe.

## Data Files

### init_prev_hiv.csv
Initial HIV prevalence by risk group, sex, and sex worker status.

**Columns:**
- `risk_group`: Risk group (0=low, 1=mid, 2=high)
- `sex`: Sex (female, male)
- `sw`: Sex worker status (0=no, 1=yes)
- `init_prev`: Initial HIV prevalence (0-1 scale)

**Source:** Calibrated to Zimbabwe epidemiological data

### condom_use.csv
Condom use by partnership type over time.

**Columns:**
- `partnership`: Partnership type (e.g., "(0,0)", "(fsw,client)")
- Year columns: Condom use proportion (0-1 scale)

**Source:** Zimbabwe Demographic and Health Surveys

### art_coverage.csv
ART coverage (absolute numbers on ART) over time.

**Columns:**
- `year`: Year
- `n_art`: Number of people on ART

**Source:** UNAIDS estimates for Zimbabwe

### vmmc_coverage.csv
Voluntary medical male circumcision coverage over time.

**Columns:**
- `year`: Year
- `coverage`: VMMC coverage (0-1 scale)

**Source:** Zimbabwe Ministry of Health reports

## Demographic Data

Demographic data (births, deaths, age distribution) is stored in `stisim/data/files/`:
- `zimbabwe_asfr.csv` - Age-specific fertility rates
- `zimbabwe_deaths.csv` - Death rates by age/sex/year
- `zimbabwe_age_1980.csv` - Age distribution for 1980

**Source:** UN World Population Prospects

## Usage

### Via stisim_examples
```python
import stisim_examples as stx

# Zimbabwe HIV model with pre-configured data
sim = stx.Sim(demographics='zimbabwe', diseases='hiv')
sim.run()

# Or use HIV-specific convenience function
sim = stx.HIVSim(location='zimbabwe')
sim.run()
```

### Via hivsim_examples
```python
import hivsim_examples as hx

sim = hx.Sim(location='zimbabwe')
sim.run()
```

### Via base stisim/hivsim
```python
import stisim as sti
import hivsim as hs

# Load only demographics, manual disease configuration
sim = sti.Sim(demographics='zimbabwe', diseases='hiv', n_agents=5000, dur=40)
sim.run()

# Or with hivsim
sim = hs.Sim(location='zimbabwe', n_agents=5000, dur=40)
sim.run()
```

### With custom data folder
```python
# Use data from external repository
sim = hs.Sim(
    location='zimbabwe',
    datafolder='/Users/robynstuart/gf/hiv_zim/data',
    n_agents=5000,
    dur=40
)
sim.run()
```

## Notes

- Data has been calibrated to match Zimbabwe's HIV epidemic trajectory
- ART coverage is provided as absolute numbers, not proportions
- Condom use varies by partnership type and has increased over time
- VMMC programs began scaling up in 2010

## References

- UNAIDS Zimbabwe Country Profile
- Zimbabwe Demographic and Health Surveys
- UN World Population Prospects
