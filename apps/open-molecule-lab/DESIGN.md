# Open Molecule Lab Prompt Workbench Design System

## 1. Atmosphere & Identity

Open Molecule Lab is a prompt-first scientific computing workbench. The first screen is the product itself: a researcher describes a scientific screening request, confirms the target and run scale, then reviews an auditable RunSpec before any computation is started. The atmosphere is calm, technical, and trustworthy rather than dashboard-heavy. The signature is a single research-prompt console paired with an evidence-led execution plan.

The product must not feel like an internal engineering cockpit. Avoid side navigation, dense agent taxonomies, decorative model metrics, and exposed implementation jargon in the primary flow. Engineering details remain available only as audit evidence: task id, artifact count, commands, gate status, and blocked post-validation items.

## 2. Color

### Palette

| Role | Token | Value | Usage |
| --- | --- | --- | --- |
| Surface/canvas | `--surface-canvas` | `#f4efea` | Page background |
| Surface/primary | `--surface-primary` | `#fffdf9` | Form and report surfaces |
| Surface/secondary | `--surface-secondary` | `#f8f4ed` | Subtle bands and disabled areas |
| Surface/ink | `--surface-ink` | `#18201d` | Result header and dark actions |
| Text/primary | `--text-primary` | `#202420` | Main text |
| Text/secondary | `--text-secondary` | `#66706a` | Supporting copy |
| Text/muted | `--text-muted` | `#8a928c` | Metadata |
| Text/inverse | `--text-inverse` | `#fffdf9` | Text on dark surfaces |
| Border/default | `--border-default` | `rgba(32, 36, 32, 0.14)` | Inputs, rows, panels |
| Border/strong | `--border-strong` | `rgba(32, 36, 32, 0.24)` | Active controls |
| Accent/primary | `--accent-primary` | `#1d7a68` | Primary CTA and focus |
| Accent/primary-hover | `--accent-primary-hover` | `#146252` | CTA hover |
| Accent/yellow | `--accent-yellow` | `#ffde00` | Small product signal and highlights |
| Accent/yellow-hover | `--accent-yellow-hover` | `#f0cf00` | Primary yellow CTA hover |
| Accent/blue | `--accent-blue` | `#2f88cc` | Links and informational status |
| Status/success | `--status-success` | `#187108` | Completed/pass |
| Status/warning | `--status-warning` | `#9a6515` | Needs review/skipped |
| Status/error | `--status-error` | `#b42318` | Failed/fail-closed |

### Rules

- The canvas is warm cream, not dark.
- Green is reserved for the primary action and completed computational steps.
- Yellow is a small signature accent, never a dominant background.
- Red appears only for fail-closed or unavailable states.
- Raw color values may appear only in this document and `src/styles.css`.

## 3. Typography

| Level | Size | Weight | Line Height | Usage |
| --- | --- | --- | --- | --- |
| Display | 44px | 760 | 1.05 | Product promise |
| H1 | 32px | 720 | 1.12 | Page and report titles |
| H2 | 22px | 700 | 1.22 | Section headings |
| H3 | 16px | 700 | 1.35 | Candidate and step titles |
| Body | 15px | 440 | 1.6 | Main copy |
| Small | 13px | 520 | 1.45 | Labels, hints |
| Data | 12px | 520 | 1.5 | SMILES, IDs, command snippets |

Font stack:
- UI: `Inter`, `Geist`, `-apple-system`, `BlinkMacSystemFont`, `Segoe UI`, `sans-serif`
- Data: `SFMono-Regular`, `ui-monospace`, `Menlo`, `Consolas`, `monospace`

## 4. Spacing & Layout

Base spacing unit is 4px.

| Token | Value | Usage |
| --- | --- | --- |
| `--space-1` | 4px | Tight icon gap |
| `--space-2` | 8px | Inline spacing |
| `--space-3` | 12px | Compact controls |
| `--space-4` | 16px | Standard padding |
| `--space-5` | 20px | Field groups |
| `--space-6` | 24px | Panel padding |
| `--space-8` | 32px | Section gap |
| `--space-10` | 40px | Large composition gap |
| `--space-14` | 56px | First-screen vertical rhythm |

Layout:
- Max content width: 1180px.
- Header is a slim brand row with service state.
- First screen is a two-column product console: request form on the left, live result/status on the right.
- Result sections stack below the first screen.
- Mobile collapses to a single column with form before result.

## 5. Components

### ProductShell
- Structure: header, main request console, result report, evidence details.
- States: runner online, runner unavailable, task running, completed with review blockers.
- Accessibility: `<main>` landmark, one H1, visible service state.

### RequestConsole
- Structure: research prompt textarea, candidate-pool number input, desired candidate count input, target id, route selector, CSV upload, plan button and run button.
- Variants: ready, submitting, disabled, error.
- States: hover, focus, loading, unavailable.
- Backend contract: submit the prompt to `/api/prompt-plan`; parse and content-address the CSV through `/api/molecule-sets`, then attach it through `/api/runs`.

### ModelSelector
- Structure: selectable module rows with model name, backend availability, and route.
- Required baseline: RDKit is locked on.
- Optional modules: Vina, GNINA, PoseBusters, Boltz-2, Chai-1, PLIP, ProLIF.
- States: selected, unselected, unavailable, planned.
- Boundary: unavailable selected modules remain fail-closed and cannot synthesize fake evidence.

### LiveRunCard
- Structure: plan id/run id, status, preflight checks, persisted prepare/score/dock/report rows, attempt count, timestamps, terminal error code, and latest real metrics.
- States: no run, queued, running, complete, failed, blocked, cancelled.
- Resume: show the `RotateCcw` command only when the stage API reports `resumable`; never infer resumability from elapsed time or log text.
- Boundary: status is workflow state only, not drug efficacy.

### CandidateReport
- Structure: headline metrics, ranked rows, SMILES, four layer statuses, final score, gate status/reason and preserved failed rows.
- Variants: empty, populated, needs review.
- Boundary: candidate rows are computational priorities only.

### EvidenceLedger
- Structure: quality gate rows, post-validation rows, blocker summary, next computational actions.
- States: pass, warn, fail, missing.
- Boundary: missing structural evidence must stay visible and fail closed.

### AuditDrawer
- Structure: task request, exact options, commands, artifact list.
- States: collapsed, expanded.
- Boundary: audit material is for traceability, not efficacy claims.

## 6. Motion & Interaction

- Motion uses transform and opacity only.
- Buttons shift 1px on press and lift 1px on hover.
- Running state uses a subtle progress sweep.
- Respect `prefers-reduced-motion`.

## 7. Depth & Surface

Strategy: warm paper surfaces with subtle borders and small shadows.

| Level | Treatment | Use |
| --- | --- | --- |
| Flat | Cream canvas | Page background |
| Contained | `1px solid var(--border-default)` | Inputs, rows |
| Elevated | subtle shadow + border | Request and report surfaces |
| Inverted | dark ink surface | Live run emphasis |

No nested card stacks. Use framed surfaces only for the request console, result report, candidate rows, and audit details.

## 8. Verified Stage Resume

- The UI renders `GET /api/runs/:id/stages` after launch and on every run poll; planned rows are used only before launch.
- Library displays dock as `skipped`; cascade requires a complete dock attempt with at least one real structure-docking success.
- Resume uses the same run ID and immutable RunSpec. It appends an attempt only after manifest, input, asset, code, upstream checkpoint and output hashes match.
- `checkpoint_mismatch` is terminal `blocked` evidence. The UI displays the persisted stage/error code and never estimates a percentage.
- Recovery is stage-boundary recovery. Mid-model and mid-docking process state are intentionally outside the product contract.
