# Composing named observation channels into a Disease Profile without illegal merge into Θ

**Thesis #22.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #3 (Disease Profile research object), Thesis #8 (named observation channel), Thesis #15 (multi-observation profile; merged θ refused).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

When do separately named observation channels compose into a Disease Profile without illegally merging into a single therapeutic parameter θ?

The legal operator is disjoint union of named maps. The composition record names the factors, leaves θ empty, and names neither a quotient nor a combiner. Identification, a product, an arithmetic mean, promotion of one coefficient into θ, a repeated channel identifier, and a join across two disease classes are refused. The JSON Schema of the legal contract still accepts the repeated identifier and the cross-class join. The rule list does not.

On a five-preparation linear toy the joint Fisher matrix is diagonal, with entries 3819.4444, 9155.5556, and 4000, and rank 3. Each channel alone has rank 1. A binary join has rank 2. Identifying the coefficients produces θ̂ = 1.3252986, against generating values α = 0.62, β = 1.15, and γ = 2.40. The profile of that symbol is closed and excludes each generating value. The product αγ = 1.488 has rank 1 and a flat profile along the hyperbola. Seed 20260922. The loadings are synthetic.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered. It does not re-tabulate the 2022 papaya assay, and it does not reuse the loadings of Thesis #15.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `schema/composed_profile.schema.json` | JSON Schema for a legal disjoint union (`1.2.0-compose`) |
| `profiles/legal_bundle.profile.yaml` | Legal join of a T08-style factor with a two-channel factor |
| `profiles/legal_singletons.profile.yaml` | The same three channels as singleton factors |
| `profiles/legal_binary.profile.yaml` | Legal binary join |
| `profiles/illegal_*.profile.yaml` | Files the rule list refuses |
| `sim/compose_profile.py` | Seeded ranks, profiles, and refusal checks (seed 20260922) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Spectra, profiles, and ranks |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/compose_profile.py
python3 build_pdf.py
```

NumPy, Matplotlib, PyYAML, and jsonschema are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json`, the profile files, and `sim/figures/`.

## Cite

Ogbonna KE. Composing named observation channels into a Disease Profile without illegal merge into Θ [Internet]. Thesis #22 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-22-observation-channel-profile-composition

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

Priors, cited as manuscripts and not re-derived here: [Thesis #3](https://github.com/cloudynirvana/thesis-03-disease-profile), [Thesis #8](https://github.com/cloudynirvana/thesis-08-papaya-agnp-observation-channel), [Thesis #15](https://github.com/cloudynirvana/thesis-15-nanobiocomposite-multiobservation-profile).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
