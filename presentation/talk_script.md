# Talk script — Spatial Architecture of Pediatric Brain Tumors (~20 min)

Target ~20 min, 20 slides. Times are cumulative-ish guides, not hard limits.
If running long, the atlas (slides 6–9) is the most compressible — you can pair the
bar + forest for each concept and move faster. Framing is correlational throughout:
say "tracks / relates to / is associated with," not "drives/causes."

---

## [1] Title — ~0:30
"Thanks everyone. Today I want to talk about the *spatial* architecture of pediatric
brain tumors — how the cells in these tumors are arranged in tissue, and what that
arrangement can tell us. "

---

## [2] Hook — why look spatially — ~1:15
"Here's the motivation, and it starts with a microscope slide. This is a pilocytic
astrocytoma on H&E. Notice it isn't uniform — a pathologist immediately sees two
textures: dense, compact-fibrillary tissue on one side, and loose, myxoid tissue on
the other. We call that *biphasic*, and it's a diagnostic hallmark.

The point is: a tumor is not a bag of cells. Its structure is spatially organized —
you can literally see it. But the tools we usually use — bulk sequencing, or
single-cell after dissociation — throw that spatial context away the moment you
grind up the tissue.

Single-cell *spatial* transcriptomics lets us keep it: we measure gene expression
while every cell stays in place. So the question this whole talk is built around is —
what recurrent tissue organization exists across these tumors, and what does it
reflect?"

---

## [3] The cohort — ~1:00
"Briefly, the cohort. Eleven tumors, four histologies, profiled on the 10x Xenium
brain panel — 266 genes. Four pilocytic astrocytomas, two gangliogliomas, two central
neurocytomas, and three cerebellar liponeurocytomas.

One honest caveat up front, and I'll come back to it: eleven samples is small, and
it's confounded — histology, driver alteration, and anatomical location all move
together. So please read everything today as hypothesis-generating spatial
biology, not a powered cohort study."

---

## [4] Four histologies — ~0:45
"Here are the four, just so you have the visual vocabulary: pilocytic astrocytoma —
biphasic, BRAF-driven; ganglioglioma — mixed neurons and glia, often BRAF-V600E;
central neurocytoma — intraventricular; and cerebellar liponeurocytoma — the lipidized
one. Four quite different tissues. Let's dive into them."

---

## [5] BANKSY — the composition lens — ~1:00
"The first way we describe the tissue is *composition* — what cell types are present,
and how they're arranged. We do that with a method called BANKSY.

The idea is simple: first we give every cell a type from its expression. Then BANKSY
clusters each cell using not just its own profile, but its neighborhood too — so cells
that sit in similar surroundings group together. What comes out are recurrent
multicellular neighborhoods — we call them *spatial niches*. So a niche is a
reproducible local mixture of cell types. That's our composition readout: which cell
types co-occur, and where."

---

## [6] Cell types across the cohort — ~0:45
"First the cell-type composition, across all eleven tumors. The usual suspects are
here — glial lineages, neurons, myeloid and lymphoid immune cells, endothelium,
stroma. What I want you to take away is that the *mix* differs by histology. Let me
make that quantitative."

---

## [7] Cell-type enrichment by histology — ~1:00
"This is a forest plot — each histology versus the rest of the cohort. And each tumor
type has its own signature cell-type makeup. Central neurocytoma is enriched for
radial glia. Pilocytic astrocytoma for OPCs. Ganglioglioma for astroglial cells.
Liponeurocytoma leans neuronal. None of that is surprising — it matches the known
biology — which is exactly what you want to see as a sanity check that the method is
reading real signal."

---

## [8] Recurrent spatial niches — ~0:45
"Now the niches — the composition side. Across the cohort we recover a set of
recurrent spatial niches. Again, look at how they're distributed."

---

## [9] Niches by histology — ~1:00
"And again — each histology has its own niches. The spatial neighborhoods aren't
shared uniformly; they're histology-specific, in the same way the cell types were.
So both the cells and their spatial arrangements are characteristic of tumor type."

---

## [10] Niches are ecosystems — Act I payoff — ~1:15
"Here's the payoff of the descriptive part. If you ask what each signature niche is
actually *made of*, you find they're ecosystems built around a signature cell type.
The central-neurocytoma niche is a radial-glia / neuronal ecosystem. The
liponeurocytoma niche is dominated by neurons — 55%. The ganglioglioma niche is
astrocyte-rich. And the pilocytic astrocytoma niche is over half OPCs.

So histology-specific niches are built from histology-specific cells. That's the
natural transition into the second half of the talk — because it raises a sharper
question. It's not surprising a PA niche is full of PA-type cells. The real question
is: are those cells *also* high in the biologically relevant pathway — MAPK, in
pilocytic astrocytoma? more over, do they correspond also to the morphlogical features a pathologist would see?
To answer that I need more than composition."

---

## [11] The pivot — three views, three questions — ~1:15
"So here's the framework for the rest of the talk. We've been describing composition
with BANKSY. But you can describe the *same cells* two other ways. ENVI gives us cell
*state* — I'll explain it in a second. And STalign gives us *morphology* — it puts the
pathologist's H&E annotations back onto our cells.

Three descriptions of one tissue — and every pairing is a question.
Q1: does a niche relate to cell state — are a niche's cells high in the relevant
pathway? Q2: does a niche match the morphology a pathologist sees? And Q3: does cell
state relate to morphology? I'll take them one at a time, starting with Q1, and Q1 is
where most of the story is. From here on I'm focusing on the pilocytic astrocytomas."

---

## [12] Q1 — niches organize glial state — ~1:30
"First, what is ENVI. We have a 266-gene panel in space, and matched single-nucleus
whole-transcriptome data from the same tumors. ENVI co-embeds them and imputes
genome-wide expression for every spatial cell — so we go from 266 measured genes to a
full-transcriptome estimate, *in place*. That lets us score pathway activity — like
MAPK — cell by cell, spatially.

When we do that: MAPK activity is concentrated in the tumor glial lineage — the OPCs
and astrocytes — and it's low in immune and stromal cells. Now, within that glial
lineage — and this is the Q1 result — MAPK is higher in specific niches. In these four
PA patients, glia sitting in the two PA-enriched niches, niches 3 and 5, are the
MAPK-high glia. So the niche a glial cell occupies tracks its state. Composition and
state are linked."

---

## [13] Q1 — MAPK-high glia and their neighbors — ~1:30
"If glial state varies by niche, does it reshape the neighborhood around those glia?
This is per patient — each panel is one of the four PA tumors, each with a different
MAPK driver. As glial MAPK rises, the local neighborhood gains myeloid cells — it's
enriched in the same direction across patients. That's a recruitment signal. And
endothelial cells move the opposite way — depleted as MAPK rises. So the
neighborhood shifts toward immune and away from vascular.

One honest point: this model adjusts for patient, but not yet for niche. Because
MAPK-high glia sit in niches 3 and 5, part of 'more myeloid neighbors' could just be
restating that those niches are myeloid-rich. A niche-adjusted model is the planned
next step to separate those. So I'll call this a strong, consistent association, with
that control still owed."

---

## [14] Q1 — an immunosuppressive-TAM shift — ~1:30
"It's not just how many myeloid cells — it's *which kind*. This is the focused result.
For myeloid cells sitting near the MAPK-high, PA-niche glia — the tumor core — versus
myeloid near other glia, four programs summarize the shift. Complement is up.
Activated microglia is up. And classical inflammatory response and interferon-gamma
response are both down. Up complement and activated microglia, down IFN-gamma and
inflammatory — that's the profile of an immunosuppressive, alternatively-activated
tumor-associated macrophage.

And look across the columns: same direction in all four patients, each a different
driver. The inferential unit here is the sample — four patients — not the individual
cells; per-cell statistics would be pseudoreplicated. On the honest sample-level test,
complement and activated microglia survive correction — those carry the stars. And the
effect persists when we adjust for the glia's own MAPK dose, so it tracks PA-niche
proximity itself. It's suggestive at n=4, but it's coherent biology: MAPK-high glia
sit in an immunosuppressive myeloid microenvironment."

---

## [15] Q2 setup — do niches match the pathologist? — ~1:00
"That's Q1 — composition and state are linked. Q2 asks something a pathologist would
ask: do our niches line up with what they annotate on the H&E — compact versus myxoid?

To even ask that, we have to put the pathologist's annotations onto our cells, and
that's STalign — it registers the H&E and its region annotations onto the spatial data,
cell by cell. Does the registration work? Here's the control: our endothelial cells
land squarely on the blood vessels a pathologist can see on the H&E. So the alignment
is real. Now we can ask the question directly."

---

## [16] Q2 result — the null — ~0:45
"And the answer is no. Compact-fibrillary versus loose-myxoid regions show no enriched
or depleted niches — everything is non-significant. So the transcriptional niches we
find and the morphological regions a pathologist draws are *not* the same partition of
the tissue. That's a genuine, if slightly deflating, result — and it sets up a puzzle.
If composition doesn't explain the morphology, what does?"

---

## [17] Q3 setup — does cell state differ by region? — ~0:45
"So we go to the third edge. Composition didn't separate compact from myxoid — so
let's look at cell *state* directly. Same H&E: compact on the green side, loose myxoid
on the blue side. The question is simply — do the cells in these two regions differ in
their state, holding cell type and niche fixed so we're comparing like with like?"

---

## [18] Q3 result — myxoid carries a glycan-ECM state — ~1:15
"They do. Across OPCs, astrocytes, and myeloid cells, the myxoid regions turn *up*
glycosaminoglycan and proteoglycan synthesis — the machinery that builds gelatinous,
glycan-rich matrix — and turn *down* the pathways that break it down. It's a
coordinated matrix-accumulation state. And that's essentially what myxoid tissue *is*
histologically — loose, mucoid, glycan-rich matrix. So the morphology the pathologist
sees looks like a readout of a cell-state program.

Caveat, same as before: this is three PA samples, and it's partly definitional — but
it's a clean, consistent signal and it closes the loop back to that first H&E."

---

## [19] Synthesis — cell state is the common thread — ~1:15
"So let's put the triangle back together, now with answers. Q1: niche, state, and
immune infiltration are all associated — composition and state travel together. Q2:
niches do *not* match the pathologist's morphology — that edge is null. And Q3: cell
state *does* track morphology.

The way I read this: cell state is the common thread. Composition and morphology each
reflect the underlying state, but they read *different axes* of it — composition reads
an immune axis, morphology reads a matrix, secretory axis. That's exactly why Q2 can be
null even though Q1 and Q3 hold: two honest descriptions of the same tissue don't have
to agree, because they're reading different things. And it brings us full circle — the
biphasic appearance we started with, on that very first H&E, is what the tissue looks
like when the cell state deposits glycan-rich myxoid matrix."

---

## [20] Acknowledgements — ~0:30
"With that — thank you. This is a collaboration with the Sahm and Sievers groups at
DKFZ and Heidelberg, with Dilan Savran, and with Jared Collins, Pratiti Bandopadhayay,
and my advisor Rameen Beroukhim, plus Natalie and Sabrina. Happy to take questions."

---

## Anticipated Q&A (keep in your back pocket)
- **"Aren't you underpowered?"** — Yes, n=11 / n=4 PA, confounded. Every claim is
  direction + cross-sample consistency, not a p-value. The PA per-cell models use the
  sample as the inferential unit.
- **"Is the myeloid recruitment just niche composition?"** — Possibly in part; the
  current model is patient-adjusted, niche-adjustment is the planned control. The
  reprogramming (TAM axis) does survive adjusting for glial MAPK dose.
- **"How reliable is ENVI imputation on 266 genes?"** — It's imputed, so treat as
  hypothesis-generating; MAPK is concentrated in the expected lineage, which is a
  sanity check. Would validate key programs orthogonally.
- **"Why is Q2 null — bad registration?"** — Registration is validated (endothelial
  cells on vessels). The null is real: niches and morphology partition the tissue
  differently. Q3 explains what morphology *does* track (state).
- **"Causality?"** — We only show associations. The triangle edges are correlational;
  "state is the hub" is an interpretation, not a demonstrated mechanism.
