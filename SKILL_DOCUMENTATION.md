---
name: vibrio-natriegens-modeling
description: "Test genomic variant hypotheses in Vibrio natriegens using GSMM (iLC858_v1.1) and RBA models. Gene knockouts, mutations, pathway engineering, metabolic flux analysis for fastest-growing organism."
---

# Vibrio natriegens Metabolic Modeling

## When to Use This Skill

Use this skill when you need to:
- **Test genomic variant hypotheses** - Predict how gene deletions, mutations, or additions affect V. natriegens metabolism
- **Analyze metabolic phenotypes** - Understand growth characteristics, substrate utilization, or metabolic capabilities
- **Design metabolic engineering strategies** - Evaluate pathway modifications, heterologous gene expression, or enzyme improvements
- **Compare modeling approaches** - Contrast constraint-based GSMM vs. proteome-constrained RBA predictions
- **Study fast-growth mechanisms** - Investigate metabolic strategies of the fastest-growing organism

**Prerequisites**: This skill assumes you have access to the `cobrapy` skill for basic constraint-based modeling operations.

## Quick Start

### Testing a Gene Knockout Hypothesis

```python
import cobra

# Load V. natriegens GSMM
model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')

# Get baseline
wt_growth = model.optimize().objective_value
print(f"Wild-type growth: {wt_growth:.4f} h⁻¹")

# Test cytochrome bo3 knockout
with model:
    model.genes.get_by_id("PN96_21415").knock_out()
    ko_growth = model.optimize().objective_value
    print(f"Knockout growth: {ko_growth:.4f} h⁻¹")
    print(f"Growth impact: {(ko_growth/wt_growth - 1)*100:+.1f}%")
```

### Running RBA Simulation

```python
import sys
sys.path.append("/path/to/RBApy")  # Update with your path
import rba

# Load RBA model (requires prior generation)
model = rba.RbaModel.from_xml("RBA/model")

# Solve for optimal growth
result = model.solve(bissection_tol=0.01)
print(f"RBA growth rate: {result.mu_opt:.4f} h⁻¹")

# Export fluxes
result.write_fluxes('rba_results.csv', file_type="csv")
```

## Overview

This skill provides specialized capabilities for working with *Vibrio natriegens* genome-scale models to test metabolic hypotheses about genomic variants. V. natriegens is the fastest-growing organism known (doubling time <10 minutes), making it an important model for studying rapid growth metabolism.

**Available Models:**
- **iLC858_v1.1 (GSMM)** - 858 genes, ~1000 reactions, optimized for constraint-based modeling
- **RBA Model** - Proteome-constrained model accounting for ribosome capacity, protein allocation, and macromolecular composition

**Model Features:**
- Corrected electron transport chain (multiple terminal oxidases)
- Accurate quinone metabolism (ubiquinone/menaquinone)
- PHB (polyhydroxybutyrate) metabolism
- ATP maintenance requirements
- Validated growth on multiple carbon sources

## Core Capabilities

### 1. Genomic Variant Testing Scenarios

#### Gene Deletion (Knockout)

Test the effect of removing one or more genes:

```python
# Single gene knockout
with model:
    gene = model.genes.get_by_id("PN96_21415")  # Cytochrome bo3 subunit
    gene.knock_out()
    solution = model.optimize()
    print(f"Growth: {solution.objective_value:.4f} h⁻¹")

# Operon deletion (multiple genes)
operon_genes = ["PN96_21415", "PN96_21420", "PN96_21425", "PN96_21430"]
with model:
    for gene_id in operon_genes:
        model.genes.get_by_id(gene_id).knock_out()
    solution = model.optimize()
    print(f"Operon KO growth: {solution.objective_value:.4f} h⁻¹")

# Systematic knockout screen (use cobrapy skill)
from cobra.flux_analysis import single_gene_deletion
results = single_gene_deletion(model)
essential = results[results['growth'] < 0.01]
print(f"Essential genes: {len(essential)}")
```

#### Enzyme Efficiency Changes (Mutations)

Model beneficial or deleterious enzyme mutations:

```python
# Simulate 30% efficiency improvement in cytochrome bd
rxn = model.reactions.get_by_id("rxn10806_c")
original_ub = rxn.upper_bound

rxn.upper_bound = original_ub * 1.3  # 30% increase
solution = model.optimize()
print(f"Growth with improved enzyme: {solution.objective_value:.4f} h⁻¹")

# For RBA: modify kapp values
import pandas as pd
kapp_data = pd.read_csv('RBA/data/kapp_resp.tsv', sep='\t')
kapp_data.loc[kapp_data['reaction'] == 'R_rxn10806_c', 'kapp'] *= 1.3
kapp_data.to_csv('RBA/data/kapp_mutant.tsv', sep='\t', index=False)

rba_model = rba.RbaModel.from_xml("RBA/model")
rba_model.set_enzyme_efficiencies('RBA/data/kapp_mutant.tsv')
result = rba_model.solve()
print(f"RBA growth: {result.mu_opt:.4f} h⁻¹")
```

#### Pathway Addition (Horizontal Gene Transfer)

Add heterologous metabolic capabilities:

```python
from cobra import Reaction, Metabolite

# Example: Add lactose utilization
lactose_c = Metabolite('cpd00192_c', formula='C12H22O11',
                       name='Lactose', compartment='c')

# Beta-galactosidase reaction
lac_hydrolysis = Reaction('LAC_hydrolysis')
lac_hydrolysis.name = 'Lactose hydrolysis'
lac_hydrolysis.lower_bound = 0
lac_hydrolysis.upper_bound = 1000
lac_hydrolysis.gene_reaction_rule = "lacZ"  # Heterologous gene

lac_hydrolysis.add_metabolites({
    lactose_c: -1,
    model.metabolites.get_by_id('cpd00001_c'): -1,  # H2O
    model.metabolites.get_by_id('cpd00027_c'): 1,   # Glucose
    model.metabolites.get_by_id('cpd00559_c'): 1    # Galactose
})

# Add transport and exchange reactions
model.add_reactions([lac_hydrolysis])
model.add_boundary(lactose_c, type="exchange", lb=-10, ub=1000)

# Test growth on lactose
model.reactions.EX_cpd00027_e.lower_bound = 0  # Block glucose
solution = model.optimize()
print(f"Growth on lactose: {solution.objective_value:.4f} h⁻¹")
```

#### Gene Overexpression

Simulate pathway upregulation:

```python
# Overexpress PHB synthesis genes
phb_reactions = ['rxn01453_c', 'rxn27735_c']

for rxn_id in phb_reactions:
    rxn = model.reactions.get_by_id(rxn_id)
    rxn.upper_bound *= 2.0  # 2x overexpression

# Add PHB production objective
phb_demand = Reaction('DM_PHB')
phb_demand.add_metabolites({
    model.metabolites.get_by_id('cpd15380_c'): -1
})
phb_demand.bounds = (0, 1000)
model.add_reactions([phb_demand])

# Optimize for PHB production
model.objective = 'DM_PHB'
solution = model.optimize()
print(f"PHB production: {solution.objective_value:.4f} mmol/gDW/h")
```

#### Metabolic Burden Analysis

Quantify the cost of heterologous expression:

```python
# Add synthetic protein production drain
synthetic_protein = Reaction('SYNTH_PROTEIN')
synthetic_protein.name = 'Synthetic protein burden'
synthetic_protein.lower_bound = 1.0  # Constitutive expression
synthetic_protein.upper_bound = 1.0

# Simplified amino acid and ATP cost
synthetic_protein.add_metabolites({
    model.metabolites.get_by_id('cpd00002_c'): -100,  # ATP
    model.metabolites.get_by_id('cpd00062_c'): -100,  # GTP
    model.metabolites.get_by_id('cpd00035_c'): -10,   # Alanine
    model.metabolites.get_by_id('cpd00041_c'): -10,   # Aspartate
    # Add more amino acids as needed
    model.metabolites.get_by_id('cpd00008_c'): 100,   # ADP
    model.metabolites.get_by_id('cpd00015_c'): 100,   # GDP
})

model.add_reactions([synthetic_protein])
solution = model.optimize()

burden = (wt_growth - solution.objective_value) / wt_growth * 100
print(f"Metabolic burden: {burden:.2f}%")
```

### 2. Model-Specific Operations

#### Working with V. natriegens Gene IDs

Gene locus tags follow the **PN96_xxxxx** format:

```python
# Find reactions for a specific gene
gene = model.genes.get_by_id("PN96_09285")
print(f"Gene {gene.id} participates in:")
for rxn in gene.reactions:
    print(f"  {rxn.id}: {rxn.name}")
    print(f"    GPR: {rxn.gene_reaction_rule}")
```

#### Key V. natriegens Metabolic Reactions

Important reactions for hypothesis testing:

```python
# Electron transport chain terminal oxidases
etc_reactions = {
    'rxn10113_c': 'Cytochrome bo3 oxidase (high O2)',
    'rxn10806_c': 'Cytochrome bd oxidase (low O2)',
    'rxn14426_c': 'Cytochrome cbb3 oxidase',
    'rxn19357_c': 'Cytochrome aa3 oxidase',
    'rxn35348_c': 'Ubiquinol-cytochrome c reductase (bc1 complex)'
}

# PHB metabolism
phb_reactions = {
    'rxn01453_c': '3-hydroxybutanoyl-CoA dehydrogenase',
    'rxn27735_c': 'PHB synthesis'
}

# Key metabolic nodes
central_reactions = {
    'rxn00288_c': 'Glucose-6-phosphate isomerase',
    'rxn01512_c': 'Citrate synthase',
    'rxn10042_c': 'ATP synthase',
    'rxn00062_c': 'ATP maintenance'
}

# Check fluxes through these reactions
solution = model.optimize()
for rxn_id, description in etc_reactions.items():
    flux = solution.fluxes.get(rxn_id, 0)
    print(f"{rxn_id} ({description}): {flux:.4f}")
```

#### Growth Medium Configuration

V. natriegens can grow on various carbon sources:

```python
# Standard glucose medium
def set_glucose_medium(model, uptake_rate=10):
    """Configure minimal glucose medium."""
    # Block all carbon sources first
    for rxn in model.exchanges:
        if 'carbon' in rxn.name.lower():
            rxn.lower_bound = 0

    model.reactions.EX_cpd00027_e.lower_bound = -uptake_rate  # Glucose
    model.reactions.EX_cpd00007_e.lower_bound = -1000  # O2
    model.reactions.EX_cpd00009_e.lower_bound = -1000  # Phosphate
    model.reactions.EX_cpd00013_e.lower_bound = -1000  # NH3

# Acetate medium
def set_acetate_medium(model, uptake_rate=10):
    """Configure acetate medium."""
    set_glucose_medium(model, 0)  # Clear glucose
    model.reactions.EX_cpd00029_e.lower_bound = -uptake_rate  # Acetate

# Anaerobic conditions
def set_anaerobic(model):
    """Remove oxygen availability."""
    model.reactions.EX_cpd00007_e.lower_bound = 0

# Test under multiple conditions
conditions = {
    'Aerobic glucose': lambda m: set_glucose_medium(m, 10),
    'Aerobic acetate': lambda m: set_acetate_medium(m, 10),
    'Anaerobic glucose': lambda m: (set_glucose_medium(m, 10), set_anaerobic(m))
}

for cond_name, setup_func in conditions.items():
    with model:
        setup_func(model)
        sol = model.optimize()
        print(f"{cond_name}: {sol.objective_value:.4f} h⁻¹")
```

### 3. RBA Model Operations

RBA models require additional setup and provide more realistic proteome-constrained predictions.

#### RBA Model Generation

```python
# From RBA/ directory
# IMPORTANT: Update paths in generate_model.py first!

# Edit generate_model.py:
# sys.path.append("/your/path/to/RBApy")

# Then generate model:
# python generate_model.py

# This creates model files in RBA/model/:
# - metabolism.xml (metabolic network from iLC858_v1.1)
# - enzymes.xml (catalytic parameters)
# - proteins.xml (protein composition)
# - rnas.xml (RNA transcripts)
# - processes.xml (macromolecular processes)
# - parameters.xml (model parameters)
```

#### Solving RBA Models

```python
import sys
sys.path.append("/path/to/RBApy")
import rba

# Load model
model = rba.RbaModel.from_xml("RBA/model")

# Solve (finds maximum growth rate)
result = model.solve(bissection_tol=0.01)
print(f"Optimal growth rate: {result.mu_opt:.4f} h⁻¹")

# Access results
result.write_fluxes('results/fluxes.csv', file_type="csv")

# Note: RBA optimization is slower than GSMM FBA
# Typical solve time: 10-60 seconds vs <1 second for FBA
```

#### Modifying RBA Enzyme Efficiencies

```python
import pandas as pd

# Load existing efficiencies
kapp_df = pd.read_csv('RBA/data/kapp_resp.tsv', sep='\t')
print("Current enzyme efficiencies:")
print(kapp_df.head())

# Modify specific enzyme (simulate mutation)
mask = kapp_df['reaction'] == 'R_rxn10806_c'
kapp_df.loc[mask, 'kapp'] *= 1.5  # 50% improvement

# Save modified file
kapp_df.to_csv('RBA/data/kapp_variant.tsv', sep='\t', index=False)

# Load model with modified efficiencies
model = rba.RbaModel.from_xml("RBA/model")
model.set_enzyme_efficiencies('RBA/data/kapp_variant.tsv')
result = model.solve()
print(f"Growth with modified enzyme: {result.mu_opt:.4f} h⁻¹")

# To knockout a reaction in RBA: set kapp to 0
kapp_df.loc[mask, 'kapp'] = 0  # Complete knockout
```

#### RBA Medium Configuration

```python
# Medium is defined in RBA/data/medium.tsv
# Format: exchange_reaction_id <tab> uptake_rate

# Modify medium programmatically:
import pandas as pd

medium_df = pd.read_csv('RBA/data/medium.tsv', sep='\t',
                        names=['reaction', 'uptake'])

# Change glucose uptake
medium_df.loc[medium_df['reaction'] == 'EX_cpd00027_e', 'uptake'] = 15

# Save modified medium
medium_df.to_csv('RBA/data/medium_modified.tsv', sep='\t',
                 index=False, header=False)

# Load model with new medium
model = rba.RbaModel.from_xml("RBA/model")
model.set_medium('RBA/data/medium_modified.tsv')
result = model.solve()
```

### 4. Comparative Analysis

#### GSMM vs RBA Comparison

```python
import cobra
import rba
import pandas as pd

# Test same variant in both models
variant_name = "cytochrome_bo3_knockout"

# GSMM simulation
gsmm = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')
gsmm_wt = gsmm.optimize().objective_value

with gsmm:
    gsmm.genes.get_by_id("PN96_21415").knock_out()
    gsmm_ko = gsmm.optimize().objective_value

# RBA simulation
kapp_df = pd.read_csv('RBA/data/kapp_resp.tsv', sep='\t')
kapp_df.loc[kapp_df['reaction'] == 'R_rxn10113_c', 'kapp'] = 0
kapp_df.to_csv('RBA/data/kapp_ko.tsv', sep='\t', index=False)

rba_model = rba.RbaModel.from_xml("RBA/model")
rba_wt = rba_model.solve().mu_opt

rba_model.set_enzyme_efficiencies('RBA/data/kapp_ko.tsv')
rba_ko = rba_model.solve().mu_opt

# Compare results
comparison = pd.DataFrame({
    'Model': ['GSMM', 'RBA'],
    'WT_Growth': [gsmm_wt, rba_wt],
    'KO_Growth': [gsmm_ko, rba_ko],
    'Growth_Ratio': [gsmm_ko/gsmm_wt, rba_ko/rba_wt]
})

print(comparison)

# Interpretation:
# - GSMM often predicts higher growth rates (no proteome constraints)
# - RBA provides more realistic predictions accounting for protein costs
# - Growth ratios (KO/WT) may differ if variant affects proteome allocation
```

## Common Workflows

### Workflow 1: Test Gene Knockout Hypothesis

```python
import cobra
import pandas as pd

# Hypothesis: "PN96_21415 deletion reduces aerobic growth"

# Load model
model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')

# Establish baseline
wt_solution = model.optimize()
wt_growth = wt_solution.objective_value
print(f"WT growth: {wt_growth:.4f} h⁻¹")

# Test knockout
with model:
    gene = model.genes.get_by_id("PN96_21415")
    print(f"Gene affects {len(gene.reactions)} reactions:")
    for rxn in gene.reactions:
        print(f"  - {rxn.id}: {rxn.name}")

    gene.knock_out()
    ko_solution = model.optimize()
    ko_growth = ko_solution.objective_value

# Analyze results
growth_ratio = ko_growth / wt_growth
print(f"\nKnockout growth: {ko_growth:.4f} h⁻¹")
print(f"Growth ratio: {growth_ratio:.3f}")
print(f"Growth change: {(growth_ratio - 1)*100:+.1f}%")

# Interpret
if ko_growth < 0.01:
    print("✓ Essential gene - lethal knockout")
elif growth_ratio < 0.5:
    print("✓ Severely deleterious - major growth defect")
elif growth_ratio < 0.9:
    print("✓ Deleterious - significant growth impact")
else:
    print("⚠ Minimal impact - gene may be non-essential")

# Check flux redistribution
flux_changes = ko_solution.fluxes - wt_solution.fluxes
major_changes = flux_changes[flux_changes.abs() > 1.0].sort_values()
print(f"\nMajor flux changes (>1.0 mmol/gDW/h):")
print(major_changes)
```

### Workflow 2: Multi-Condition Variant Testing

```python
import cobra
import pandas as pd
import matplotlib.pyplot as plt

# Test variant across multiple growth conditions

model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')
gene_to_test = "PN96_21415"

# Define conditions
test_conditions = [
    ('Aerobic_Glucose', {'EX_cpd00027_e': -10, 'EX_cpd00007_e': -1000}),
    ('Aerobic_Acetate', {'EX_cpd00029_e': -10, 'EX_cpd00007_e': -1000}),
    ('Microaerobic_Glucose', {'EX_cpd00027_e': -10, 'EX_cpd00007_e': -5}),
    ('Anaerobic_Glucose', {'EX_cpd00027_e': -10, 'EX_cpd00007_e': 0}),
]

results = []

for condition_name, settings in test_conditions:
    # Configure medium
    with model:
        # Block all carbon
        for rxn in model.exchanges:
            if 'cpd00027' in rxn.id or 'cpd00029' in rxn.id:
                rxn.lower_bound = 0

        # Apply condition settings
        for rxn_id, bound in settings.items():
            if rxn_id in model.reactions:
                model.reactions.get_by_id(rxn_id).lower_bound = bound

        # WT growth
        wt_growth = model.optimize().objective_value

        # Knockout growth
        model.genes.get_by_id(gene_to_test).knock_out()
        ko_growth = model.optimize().objective_value

        results.append({
            'Condition': condition_name,
            'WT_Growth': wt_growth,
            'KO_Growth': ko_growth,
            'Ratio': ko_growth / wt_growth if wt_growth > 0 else 0
        })

# Create results dataframe
df = pd.DataFrame(results)
print(df)

# Visualize
fig, ax = plt.subplots(figsize=(10, 6))
x = range(len(df))
width = 0.35

ax.bar([i - width/2 for i in x], df['WT_Growth'], width,
       label='Wild-type', color='#2E86AB')
ax.bar([i + width/2 for i in x], df['KO_Growth'], width,
       label=f'{gene_to_test} KO', color='#F24236')

ax.set_xlabel('Growth Condition', fontsize=12)
ax.set_ylabel('Growth Rate (h⁻¹)', fontsize=12)
ax.set_title(f'Growth Impact of {gene_to_test} Deletion', fontsize=14)
ax.set_xticks(x)
ax.set_xticklabels(df['Condition'], rotation=45, ha='right')
ax.legend()
plt.tight_layout()
plt.savefig(f'results/{gene_to_test}_knockout.png', dpi=300)
plt.show()
```

### Workflow 3: Compensatory Mutation Analysis

```python
import cobra

# Hypothesis: "Improved cytochrome bd can compensate for bo3 deletion"

model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')

# Baseline
wt_growth = model.optimize().objective_value
wt_bd_flux = model.optimize().fluxes['rxn10806_c']

# Knockout bo3
ko_model = model.copy()
ko_model.genes.get_by_id("PN96_21415").knock_out()
ko_growth = ko_model.optimize().objective_value
ko_bd_flux = ko_model.optimize().fluxes['rxn10806_c']

# Compensatory mutation (increase bd efficiency)
comp_model = ko_model.copy()
bd_rxn = comp_model.reactions.get_by_id('rxn10806_c')
bd_rxn.upper_bound *= 1.5  # 50% improvement
comp_growth = comp_model.optimize().objective_value
comp_bd_flux = comp_model.optimize().fluxes['rxn10806_c']

# Analysis
print("=== Compensatory Mutation Analysis ===")
print(f"WT growth: {wt_growth:.4f} h⁻¹, bd flux: {wt_bd_flux:.4f}")
print(f"KO growth: {ko_growth:.4f} h⁻¹, bd flux: {ko_bd_flux:.4f}")
print(f"Compensated: {comp_growth:.4f} h⁻¹, bd flux: {comp_bd_flux:.4f}")

growth_lost = wt_growth - ko_growth
growth_recovered = comp_growth - ko_growth
compensation = (growth_recovered / growth_lost) * 100 if growth_lost > 0 else 0

print(f"\nGrowth lost by knockout: {growth_lost:.4f} h⁻¹")
print(f"Growth recovered by compensation: {growth_recovered:.4f} h⁻¹")
print(f"Compensation level: {compensation:.1f}%")

if compensation > 90:
    print("✓ Strong compensation - mutation nearly fully restores growth")
elif compensation > 50:
    print("⚠ Partial compensation - mutation helps but doesn't fully restore")
else:
    print("✗ Weak compensation - minimal rescue effect")
```

### Workflow 4: Pathway Engineering Design

```python
import cobra
from cobra import Reaction, Metabolite

# Design pathway for novel product synthesis

model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')

# Example: Engineer 1,3-propanediol (1,3-PDO) production
# Pathway: Glycerol → 3-hydroxypropionaldehyde → 1,3-propanediol

# Add 1,3-PDO metabolite
pdo = Metabolite('cpd00141_c', formula='C3H8O2',
                 name='1,3-Propanediol', compartment='c')

# Add glycerol dehydratase
gly_dehydratase = Reaction('GLY_DEHYDRATASE')
gly_dehydratase.name = 'Glycerol dehydratase'
gly_dehydratase.bounds = (0, 1000)
gly_dehydratase.gene_reaction_rule = "dhaB1 and dhaB2"  # Heterologous

gly_dehydratase.add_metabolites({
    model.metabolites.get_by_id('cpd00100_c'): -1,  # Glycerol
    Metabolite('cpd15441_c', name='3-HPA', compartment='c'): 1,  # Intermediate
    model.metabolites.get_by_id('cpd00001_c'): 1   # H2O
})

# Add 1,3-propanediol oxidoreductase
pdo_reductase = Reaction('PDO_REDUCTASE')
pdo_reductase.name = '1,3-Propanediol reductase'
pdo_reductase.bounds = (0, 1000)
pdo_reductase.gene_reaction_rule = "dhaT"

pdo_reductase.add_metabolites({
    model.metabolites.get_by_id('cpd15441_c'): -1,  # 3-HPA
    model.metabolites.get_by_id('cpd00004_c'): -1,  # NADH
    pdo: 1,  # 1,3-PDO
    model.metabolites.get_by_id('cpd00003_c'): 1   # NAD+
})

# Add to model
model.add_reactions([gly_dehydratase, pdo_reductase])

# Add exchange reactions
model.add_boundary(pdo, type='demand', reaction_id='DM_1_3_PDO')
model.add_boundary(model.metabolites.get_by_id('cpd00100_c'),
                   type='exchange', lb=-10, ub=1000)

# Test production
model.objective = 'DM_1_3_PDO'

# On glycerol
model.reactions.EX_cpd00027_e.lower_bound = 0  # Block glucose
model.reactions.EX_cpd00100_e.lower_bound = -10  # Allow glycerol

solution = model.optimize()
print(f"1,3-PDO production: {solution.objective_value:.4f} mmol/gDW/h")

# Check growth impact
model.objective = 'bio1'
growth = model.optimize().objective_value
print(f"Growth rate: {growth:.4f} h⁻¹")

# Production envelope (use cobrapy skill)
from cobra.flux_analysis import production_envelope
envelope = production_envelope(
    model,
    reactions=['DM_1_3_PDO'],
    objective='bio1'
)
print("\nProduction vs growth trade-off:")
print(envelope)
```

### Workflow 5: Essential Gene Identification

```python
import cobra
from cobra.flux_analysis import single_gene_deletion
import pandas as pd

# Identify essential genes under specific conditions

model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')

# Set condition (e.g., minimal glucose medium)
model.reactions.EX_cpd00027_e.lower_bound = -10
model.reactions.EX_cpd00007_e.lower_bound = -1000

# Perform systematic knockout
print("Running systematic gene deletion analysis...")
results = single_gene_deletion(model)

# Classify genes
essential = results[results['growth'] < 0.01]
severe = results[(results['growth'] >= 0.01) & (results['growth'] < 0.5)]
moderate = results[(results['growth'] >= 0.5) & (results['growth'] < 0.9)]

print(f"\nGene classification:")
print(f"Essential (growth < 1%): {len(essential)} genes")
print(f"Severe defect (1-50%): {len(severe)} genes")
print(f"Moderate defect (50-90%): {len(moderate)} genes")

# Export results
results.to_csv('results/gene_essentiality.csv')

# Identify essential genes in specific pathways
print("\nEssential respiratory genes:")
resp_genes = [g for g in essential.index
              if any('PN96_21' in g or 'PN96_08' in g or 'PN96_05' in g
                     for g in essential.index)]
for gene_id in resp_genes:
    gene = model.genes.get_by_id(gene_id)
    print(f"  {gene_id}: {[r.name for r in gene.reactions]}")
```

## Key Concepts

### Gene Nomenclature

V. natriegens genes use the **PN96_** prefix:
- Format: `PN96_xxxxx` (5-digit locus tag)
- Example: `PN96_21415` (cytochrome bo3 subunit)
- Gene-reaction rules use these locus tags
- Multiple genes can catalyze the same reaction (isozymes)

### Model Versions

**iLC858_v1.1** includes critical updates over iLC858:
- Corrected electron transport chain stoichiometry
- Fixed quinone metabolism (ubiquinone/menaquinone)
- Updated PHB synthesis genes
- Added ATP maintenance requirement (3.15 mmol/gDW/h)
- Thermodynamic corrections (irreversible reactions)
- Always use v1.1 for new analyses

### Metabolite Identifiers

Models use ModelSEED compound IDs:
- Format: `cpd#####_compartment`
- Common metabolites:
  - `cpd00002_c`: ATP (cytoplasm)
  - `cpd00008_c`: ADP
  - `cpd00006_c`: NADH
  - `cpd00027_c`: Glucose
  - `cpd00029_c`: Acetate
  - `cpd00007_e`: O2 (extracellular)
- Compartments: `c` (cytoplasm), `e` (extracellular)

### Reaction Identifiers

Reactions use ModelSEED IDs:
- Format: `rxn#####_compartment`
- Exchange reactions: `EX_cpd#####_e`
- Positive flux = product formation
- Negative flux = substrate consumption

### Growth Rate Interpretation

V. natriegens growth rates (aerobic, rich medium):
- Experimental: ~2.0 h⁻¹ (doubling time ~20 min)
- Model predictions: 1.5-2.5 h⁻¹ (depending on conditions)
- GSMM typically predicts higher than RBA
- Essential genes: growth < 0.01 h⁻¹

### GSMM vs RBA Trade-offs

| Feature | GSMM (iLC858_v1.1) | RBA Model |
|---------|-------------------|-----------|
| Speed | Fast (<1s) | Slow (10-60s) |
| Realism | Metabolic constraints only | + Proteome constraints |
| Use case | Rapid screening | Detailed mechanistic |
| Growth prediction | Often high | More realistic |
| Gene expression | Not modeled | Implicitly modeled |
| Setup complexity | Simple | Requires generation |

**When to use each:**
- **GSMM**: Initial screening, large-scale knockouts, pathway design, quick iteration
- **RBA**: Detailed mechanistic understanding, protein burden effects, growth rate optimization

## Best Practices

### 1. Always Establish Wild-Type Baseline

```python
# Before testing variants
wt_growth = model.optimize().objective_value
wt_fluxes = model.optimize().fluxes.copy()

# Store for comparison
baseline = {'growth': wt_growth, 'fluxes': wt_fluxes}
```

### 2. Use Context Managers for Testing

```python
# Ensures model isn't permanently modified
with model:
    # Make temporary changes
    model.genes.get_by_id("PN96_21415").knock_out()
    test_result = model.optimize()
# Model automatically reverted
```

### 3. Test Under Multiple Conditions

Don't rely on single-condition tests:
```python
conditions = ['aerobic_glucose', 'aerobic_acetate',
              'microaerobic', 'anaerobic']
# Test variant under all conditions
```

### 4. Validate Biological Plausibility

```python
# Check for unrealistic results
if solution.objective_value > 3.0:
    print("⚠ Warning: Unrealistically high growth rate")
    print("Check for thermodynamic loops or model errors")

# Use loopless FBA for validation
from cobra.flux_analysis.loopless import loopless_solution
loopless_sol = loopless_solution(model)
```

### 5. Document All Modifications

```python
# Keep modification log
modifications = []

def log_change(change_type, description, genes=None, reactions=None):
    modifications.append({
        'type': change_type,
        'description': description,
        'genes': genes or [],
        'reactions': reactions or []
    })

# Example
log_change('knockout', 'Deleted cytochrome bo3',
           genes=['PN96_21415'], reactions=['rxn10113_c'])
```

### 6. Check Flux Variability

Use FVA to identify alternative pathways:
```python
from cobra.flux_analysis import flux_variability_analysis

fva = flux_variability_analysis(model, fraction_of_optimum=0.95)

# Identify flexible reactions
flexible = fva[(fva['maximum'] - fva['minimum']) > 1.0]
print(f"Reactions with alternative flux states: {len(flexible)}")
```

### 7. Compare GSMM and RBA Predictions

For important hypotheses, validate with both models:
```python
# GSMM prediction
gsmm_growth = gsmm_model.optimize().objective_value

# RBA prediction
rba_growth = rba_model.solve().mu_opt

# If predictions differ significantly, investigate why
if abs(gsmm_growth - rba_growth) / gsmm_growth > 0.2:
    print("⚠ GSMM and RBA predictions differ by >20%")
    print("Consider proteome allocation constraints")
```

### 8. Visualize Results

Always visualize multi-condition comparisons:
```python
import matplotlib.pyplot as plt

# Bar plots for growth comparisons
# Heatmaps for flux changes
# Line plots for production envelopes
```

## Troubleshooting

### Infeasible Solution

**Symptom**: `model.optimize().status == 'infeasible'`

**Causes & Solutions**:
```python
# 1. Check medium constraints
print("Current medium:")
for rxn in model.exchanges:
    if rxn.lower_bound < 0:
        print(f"  {rxn.id}: {rxn.lower_bound}")

# 2. Check if essential reaction was blocked
with model:
    # Try relaxing all bounds temporarily
    for rxn in model.reactions:
        rxn.bounds = (-1000, 1000)
    if model.optimize().status == 'optimal':
        print("Issue is in reaction bounds")

# 3. Check biomass reaction
biomass = model.reactions.get_by_id('bio1')
print("Biomass reaction metabolites:")
for met, coef in biomass.metabolites.items():
    print(f"  {met.id}: {coef}")
```

### Unrealistically High Growth Rate

**Symptom**: Growth rate > 3.0 h⁻¹

**Solution**: Check for thermodynamic loops
```python
from cobra.flux_analysis.loopless import loopless_solution

# Remove thermodynamically infeasible cycles
loopless = loopless_solution(model)
print(f"Loopless growth: {loopless.objective_value:.4f} h⁻¹")
```

### RBA Path Errors

**Symptom**: `ModuleNotFoundError: No module named 'rba'`

**Solution**: Update sys.path in RBA scripts
```python
# In solve_model.py and generate_model.py
import sys
import os

# Update this path to your RBApy installation
sys.path.append("/path/to/your/RBApy")
```

### Gene ID Not Found

**Symptom**: `KeyError: 'gene_id'`

**Solution**: Verify gene ID format
```python
# V. natriegens genes use PN96_ prefix
# Check available genes:
gene_ids = [g.id for g in model.genes if 'PN96_21' in g.id]
print(f"Found {len(gene_ids)} genes with PN96_21 prefix")

# Search for genes by pattern
matching_genes = model.genes.query("PN96_21", attribute='id')
```

### Slow Optimization

**Symptom**: `model.optimize()` takes >5 seconds

**Solutions**:
```python
# 1. Try different solver
model.solver = 'glpk'  # or 'cplex', 'gurobi'

# 2. Use slim_optimize for objective value only
growth = model.slim_optimize()  # Faster

# 3. Reduce model complexity
# Check number of reactions
print(f"Reactions: {len(model.reactions)}")
```

## Example: Complete Hypothesis Testing Pipeline

```python
import cobra
import pandas as pd
import matplotlib.pyplot as plt
from cobra.flux_analysis import flux_variability_analysis

# ==============================================
# HYPOTHESIS: "Cytochrome bd upregulation
# compensates for cytochrome bo3 deletion"
# ==============================================

# Step 1: Load model and establish baseline
model = cobra.io.read_sbml_model('GSMM/iLC858_v1.1.sbml')

wt_solution = model.optimize()
wt_growth = wt_solution.objective_value
wt_bo3_flux = wt_solution.fluxes['rxn10113_c']
wt_bd_flux = wt_solution.fluxes['rxn10806_c']

print("="*50)
print("WILD-TYPE BASELINE")
print("="*50)
print(f"Growth: {wt_growth:.4f} h⁻¹")
print(f"Cytochrome bo3 flux: {wt_bo3_flux:.4f}")
print(f"Cytochrome bd flux: {wt_bd_flux:.4f}")

# Step 2: Test bo3 knockout
ko_model = model.copy()
ko_model.genes.get_by_id("PN96_21415").knock_out()
ko_solution = ko_model.optimize()
ko_growth = ko_solution.objective_value
ko_bd_flux = ko_solution.fluxes['rxn10806_c']

print("\n" + "="*50)
print("CYTOCHROME BO3 KNOCKOUT")
print("="*50)
print(f"Growth: {ko_growth:.4f} h⁻¹ ({(ko_growth/wt_growth-1)*100:+.1f}%)")
print(f"Cytochrome bd flux: {ko_bd_flux:.4f} (change: {(ko_bd_flux/wt_bd_flux-1)*100:+.1f}%)")

# Step 3: Check if bd is already maxed out
fva_ko = flux_variability_analysis(
    ko_model,
    reaction_list=['rxn10806_c'],
    fraction_of_optimum=1.0
)
bd_max = fva_ko.loc['rxn10806_c', 'maximum']
print(f"Cytochrome bd maximum capacity: {bd_max:.4f}")

if abs(ko_bd_flux - bd_max) < 0.01:
    print("✓ Cytochrome bd is already at maximum capacity")
    print("  Compensation requires increased bd expression/efficiency")
else:
    print("⚠ Cytochrome bd not at maximum - pathway may be limited elsewhere")

# Step 4: Test compensatory mutation
comp_model = ko_model.copy()
bd_rxn = comp_model.reactions.get_by_id('rxn10806_c')
bd_rxn.upper_bound *= 1.5  # 50% efficiency improvement
comp_solution = comp_model.optimize()
comp_growth = comp_solution.objective_value
comp_bd_flux = comp_solution.fluxes['rxn10806_c']

print("\n" + "="*50)
print("KNOCKOUT + COMPENSATORY MUTATION")
print("="*50)
print(f"Growth: {comp_growth:.4f} h⁻¹ ({(comp_growth/wt_growth-1)*100:+.1f}%)")
print(f"Cytochrome bd flux: {comp_bd_flux:.4f}")

# Step 5: Calculate compensation level
growth_lost = wt_growth - ko_growth
growth_recovered = comp_growth - ko_growth
compensation_pct = (growth_recovered / growth_lost * 100) if growth_lost > 0 else 0

print("\n" + "="*50)
print("COMPENSATION ANALYSIS")
print("="*50)
print(f"Growth lost by knockout: {growth_lost:.4f} h⁻¹")
print(f"Growth recovered by mutation: {growth_recovered:.4f} h⁻¹")
print(f"Compensation level: {compensation_pct:.1f}%")

if compensation_pct > 90:
    print("\n✓ CONCLUSION: Strong compensation")
    print("  Cytochrome bd upregulation nearly fully restores growth")
elif compensation_pct > 50:
    print("\n⚠ CONCLUSION: Partial compensation")
    print("  Cytochrome bd helps but doesn't fully compensate")
else:
    print("\n✗ CONCLUSION: Weak compensation")
    print("  Other factors limit compensation")

# Step 6: Visualize
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Growth rates
genotypes = ['WT', 'Δbo3', 'Δbo3 + bd↑']
growth_rates = [wt_growth, ko_growth, comp_growth]
colors = ['#2E86AB', '#F24236', '#F6AE2D']

ax1.bar(genotypes, growth_rates, color=colors)
ax1.set_ylabel('Growth Rate (h⁻¹)', fontsize=12)
ax1.set_title('Growth Rate Comparison', fontsize=14)
ax1.set_ylim(0, max(growth_rates) * 1.2)

for i, rate in enumerate(growth_rates):
    ax1.text(i, rate + 0.05, f'{rate:.3f}', ha='center')

# Cytochrome fluxes
bo3_fluxes = [wt_bo3_flux, 0, 0]
bd_fluxes = [wt_bd_flux, ko_bd_flux, comp_bd_flux]

width = 0.35
x = range(3)
ax2.bar([i - width/2 for i in x], bo3_fluxes, width,
        label='Cytochrome bo3', color='#2E86AB')
ax2.bar([i + width/2 for i in x], bd_fluxes, width,
        label='Cytochrome bd', color='#F24236')

ax2.set_ylabel('Flux (mmol/gDW/h)', fontsize=12)
ax2.set_title('Terminal Oxidase Fluxes', fontsize=14)
ax2.set_xticks(x)
ax2.set_xticklabels(genotypes)
ax2.legend()

plt.tight_layout()
plt.savefig('results/compensation_analysis.png', dpi=300)
print("\n✓ Figure saved: results/compensation_analysis.png")

# Step 7: Export results
results_df = pd.DataFrame({
    'Genotype': genotypes,
    'Growth_Rate': growth_rates,
    'Bo3_Flux': bo3_fluxes,
    'Bd_Flux': bd_fluxes
})
results_df.to_csv('results/compensation_results.csv', index=False)
print("✓ Results saved: results/compensation_results.csv")
```

## References

**Model Files:**
- `GSMM/iLC858_v1.1.sbml` - V. natriegens GSMM (recommended)
- `GSMM/iLC858.sbml` - Original version (legacy)
- `RBA/model/` - RBA model XML files

**Key Documentation:**
- `README.md` - General project documentation
- `GSMM/updates/update_v1.1.py` - Model v1.1 changes
- `RBA/generate_model.py` - RBA model generation script

**External Resources:**
- COBRApy skill documentation (for general COBRA methods)
- COBRApy docs: https://cobrapy.readthedocs.io/
- BiGG Models: http://bigg.ucsd.edu/
- ModelSEED: https://modelseed.org/

**Important Gene-Reaction Mappings:**
```
Electron Transport Chain:
- PN96_21415-21430 → rxn10113_c (Cytochrome bo3)
- PN96_08165-08170 → rxn10806_c (Cytochrome bd)
- PN96_05720-05740 → rxn14426_c (Cytochrome cbb3)
- PN96_22625-22640 → rxn19357_c (Cytochrome aa3)
- PN96_11285-11295 → rxn35348_c (bc1 complex)

Central Metabolism:
- PN96_09285 → rxn00288_c (PGI)
- PN96_01345/09275 → rxn08094_c (Citrate synthase)

PHB Metabolism:
- PN96_18045 → rxn01453_c (PHB synthesis)

ATP Synthase:
- Multiple genes → rxn10042_c
```
