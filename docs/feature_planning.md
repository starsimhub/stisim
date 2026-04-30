# Feature map / roadmap / planning doc


## What are we doing here?

The file "docs/HIVsim projects.xlsx" contains a list of upcoming HIVsim research projects that the STIsim model will need to support. Our task now is to create a feature map for the highest-priority projects in this list. We will then do a gap analysis to determine which features STIsim currently has, and which ones need to be developed. We will aim to create a feature map as well as a release roadmap. In terms of dates and timing, we're currently at the end of 2026-Q1 and we're planning the next 4 quarters. Sprint length is 2 weeks. We should look at the existing list of issues (https://github.com/starsimhub/stisim/issues) and we're aiming to create a project board view into these. I'm almost inclined to create a projects subfolder within STIsim itself for leadership tracking(?), or brainstorm suggestions on alternatives.

## Iterating on the sparkling-crunching-hickey.md plan

I've now reviewed the first version of the plan. It's good but we need to refine, and we shouldn't be overindexing on the GH issues that are already there, as we'll likely need to make a whole bunch more and some of what's there will just sit in backlog and can be addressed if priorities change or if someone has spare time.

Let me start with Adam's Uganda project. I'm wondering if it would be helpful to stub out what this will look like. We'll need some low-level initial tasks (make a repo, add project plan too the repo, gather necessary data, set up baseline model). Then in terms of features: we will need a robust PrEP implementation, which we already have several GH issues around. We'll also need all the maternal transmission issues dealt with; these are currently gathered under the https://github.com/starsimhub/stisim/milestone/18 milestone but we can reoragnize these as needed. We DO have a mechanism to route STI-positive contacts to PrEP/HIV testing, via the eligibility argument that all interventions have. However, we could think about adding stubs for this in the stubbed-out model implementation that we could construct as part of the initial tasks. Ditto, we DO have an API for one intervention to trigger enrolment in another - either (1) IntvA sets something like self.positive[uids] = True and then IntvB uses an eligibility check on this, OR (2) IntvA can directly set the eligibility for B. 

Next, for the sake of project trakcing and reporting, let's include P1.2 as a completed item and track the project over 2026-Q1 to 2027-Q2 (18m). Please revew syph_dx_zim and make a note of any features that this analysis used that STIsim supports. 

For P1.3, the main feature that will be needed is not the multi-country parameterization but rahter the VOI pipelines described in that repo. The timeline for this project is to deliver by 2026-Q3.

For P1.4, let's list a new project: priority 1, pathogens = NG+CT+TV+syph, project name = something about how novel POC diagnostics might help generate demand for testing either via partner notification or otherwise. As well as tracking the resulting improvements in health outcomes on bacterial STIs and HIV, we are also interested in how the avoidance of unncessary partner notification might reduce gender-based violence. Partners on this project are WHI and LSHTM. Timeslines: deliver by end of Q2 2026. This should actually be the highest-priority project. However, in terms of feature requirements, I believe the only things missing are the GBV linkages. 

Next: P1.5 will be HIVsim Eswatini model validation. This is the project whose core aim is to show that HIVsim can replicate any analysis that an EMOD-HIV would want to do. Unsure yet what features are needed for this. Timeline: 2027-Q1 for project completion but we'll want some intermediary goals.

And P1.6 will be Emilia's project on DoxyPep. I think that the core features needed are DoxyPep. 

Then we need the overall "model hardening" goal. This is cross-cutting and less related to individual projects. Here we need to identify gaps in documentation, test coverage, multi-country deployment, model validation in various dimensions including network validation... we need to replace the abstract issues that just say something vague like "model validation" with something that we're actually testing. We shouldn't leave it all to the end but rather integrate it throughout, as diverse team members will focus on different things. 
