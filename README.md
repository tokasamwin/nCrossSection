nCrossSection
Neutron XS data interrogation tool

Generates material data for complex mixtures of compounds.
There are two major aims:
1. Provide the full feature-set required to represent engineered designs
2. Use an intuitive approach to defining materials

Materials are made using a class hierarchy:
1. Mixtures, made of compounds in a volume-fraction (e.g. mixed gas streams, pebble beds, engine block with cylinders of gas)
2. Compounds, made of elements in ratio (e.g. chemical compounds such as H2O, metal alloys, ceramics)
3. Elements, made of isotopes in ratios/enrichments (e.g. uranium at natural enrichment or 2% enrichment)
4. Isotopes

All classes can be interrogated for cross-section data, while compounds and mixtures (as materials with densities) will have macroscopic cross sections and thus different units.
