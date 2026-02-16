# Demo Location Data

This is a generic/demo location with minimal working defaults suitable for testing and tutorials.

## Data Files

### init_prev_hiv.csv
Initial HIV prevalence by risk group, sex, and sex worker status.

**Columns:**
- `risk_group`: Risk group (0=low, 1=mid, 2=high)
- `sex`: Sex (female, male)
- `sw`: Sex worker status (0=no, 1=yes)
- `init_prev`: Initial HIV prevalence (0-1 scale)

### condom_use.csv
Condom use by partnership type over time.

**Columns:**
- `partnership`: Partnership type (e.g., "(0,0)", "(fsw,client)")
- Year columns (1980, 1990, 2000, etc.): Condom use proportion (0-1 scale)

### art_coverage.csv
ART coverage scale-up over time.

**Columns:**
- `year`: Year
- `coverage`: ART coverage (0-1 scale)

### vmmc_coverage.csv
Voluntary medical male circumcision coverage over time.

**Columns:**
- `year`: Year
- `coverage`: VMMC coverage (0-1 scale)

## Usage

```python
import stisim_examples as stx

# Create demo HIV simulation
sim = stx.Sim(demographics='demo', diseases='hiv')
sim.run()

# Or use HIV-specific convenience function
sim = stx.HIVSim(location='demo')
sim.run()
```

## Notes

- Data values are illustrative and not based on any specific real-world location
- Useful for testing, tutorials, and understanding data formats
- For real modeling work, use location-specific data (zimbabwe, kenya, etc.)
