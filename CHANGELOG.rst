==========
What's new
==========

.. currentmodule:: stisim

All notable changes to the codebase are documented in this file.


Version 1.4 (2025-08-12)
--------------------------
- Add location arg and Sim class
- Update to work with Starsim v3.
- *GitHub info*: PR `148 <https://github.com/starsimhub/stisim/pull/148>`_


Version 1.3 (2025-06-27)
--------------------------
- Fixes to the pair-matching algorithm within the sexual network to better align partner ages 
- Improvements to networks, including analyzers for debut age and partner age differences
- *GitHub info*: PR `143 <https://github.com/starsimhub/stisim/pull/143>`_


Version 1.2 (2025-06-10)
--------------------------
- Improvements to networks, including analyzers for relationship duration and network degree
- Adds a `PriorPartners` network for recalling past relationships - for use in partner notification
- *GitHub info*: PR `135 <https://github.com/starsimhub/stisim/pull/135>`_


Version 1.1.2 (2025-06-03)
--------------------------
- Bugfix to calibration class for multisims
- *GitHub info*: PR `138 <https://github.com/starsimhub/stisim/pull/138>`_


Version 1.1.1 (2025-05-23)
--------------------------
- Bugfixes to calibration class and BV connector
- Replaces the `match_pairs` method of the `StructuredSexual` network with the faster option that was previously in the `FastStructuredSexual` network (now removed).
- *GitHub info*: PR `https://github.com/starsimhub/stisim/pull/124`_


Version 1.1.0 (2025-05-13)
--------------------------
- Improvements to the Calibration class: this class now inherits directly from the Starsim calibration class, so users will have access to easier parameter constraints, plotting, flexible fit functions, etc
- Generalization of the coinfection class to handle any two diseases
- Addition of a more detailed BV model
- *GitHub info*: PR `119 <https://github.com/starsimhub/stisim/pull/119>`_


Version 1.0.5 (2025-05-08)
--------------------------
- Adds results for syphilis transmission by disease stage
- *GitHub info*: PR `112 <https://github.com/starsimhub/stisim/pull/112>`_

Version 1.0.4 (2025-05-07)
--------------------------
- Adds results for overtreatment among pregnant women
- *GitHub info*: PR `111 <https://github.com/starsimhub/stisim/pull/111>`_

Version 1.0.3 (2025-04-30)
--------------------------
- Bugfixes for congenital syphilis
- *GitHub info*: PR `85 <https://github.com/starsimhub/stisim/pull/85>`_

Version 1.0.2 (2025-04-14)
--------------------------
- Bugfixes for syphilis and GUD
- *GitHub info*: PR `83 <https://github.com/starsimhub/stisim/pull/83>`_

Version 1.0.1 (2025-03-31)
--------------------------
- Track HIV prevalence for 15-49 year olds
- *GitHub info*: PR `79 <https://github.com/starsimhub/stisim/pull/79>`_

Version 1.0.0 (2024-12-11)
--------------------------
- Updates to work with Starsim v2.2.1
- *GitHub info*: PR `63 <https://github.com/starsimhub/stisim/pull/63>`_


Version 0.2.0 (2024-11-01)
--------------------------
- Updates to work with Starsim v2.0
- *GitHub info*: PR `62 <https://github.com/starsimhub/stisim/pull/62>`_


Version 0.1.0 (2024-10-02)
--------------------------
- Collection of updates related to the NG/CT/TV work
- *GitHub info*: PR `59 <https://github.com/starsimhub/stisim/pull/59>`_


Version 0.0.2 (2024-06-07)
--------------------------
- Initial version of STIsim with structured sexual networks, models of HIV and syphilis, worksflows for model calibration, and interventions for testing and treatment.
- *GitHub info*: PR `22 <https://github.com/starsimhub/stisim/pull/22>`_


Version 0.0.1 (2024-05-15)
--------------------------
- Pre-release version

