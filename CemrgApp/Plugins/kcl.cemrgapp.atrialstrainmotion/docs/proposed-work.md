# Atrial Strain Motion — work to finish the integration

Handover document, 2026-09-04. It proposes what is left before
`kcl.cemrgapp.atrialstrainmotion` can merge from `plugin/atrailstrainmotion` into `development`.

Read [`pipeline-state.md`](pipeline-state.md) first. It describes what the pipeline does now, what
the external tools are, and what one full run produced. This document assumes it.

**Every claim below was checked against the code on `plugin/atrialstrainmotion-fix-step-13-crop` on
2026-09-04.** Where a claim comes from reading the `afmotion` container rather than this
repository, it names the file in `discovery_report/recovered/code/` inside the `charlie_pipeline`
hand-over bundle.

---

## The five actions, in priority order

The order is: first what makes the results quietly wrong, then what makes a failure invisible, then
what makes the code reviewable, then robustness on a second dataset, then the external artefact.

Actions 1 and 4 share one prerequisite — a written input contract — and are best done by the same
person. Action 3 is the cheapest and touches nothing the others touch, so it can run in parallel.
Action 5 is the slowest and the riskiest; start it early and let it run alongside the rest.

| # | Action | Cost | Risk |
|---|---|---|---|
| [1](#1--derive-the-frame-count-from-the-data) | Derive the frame count from the data | medium | low |
| [2](#2--make-a-failed-afmotion-step-fail-visibly) | Make a failed `afmotion` step fail visibly | small | low |
| [3](#3--delete-the-dead-code-and-settle-the-mitral-valve-question) | Delete the dead code, settle the mitral valve | medium | low |
| [4](#4--choose-step-1s-reference-frame-explicitly) | Choose step 1's reference frame explicitly | small | low |
| [5](#5--rebuild-the-image-as-cemrgafmotion11) | Rebuild the image as `cemrg/afmotion:1.1` | large | **high** |

---

## 1 — Derive the frame count from the data

**Why this is first.** The frame count is hardcoded, inconsistently, in four places. Nothing counts
the files in `nifti/`. Both demo datasets hold 21 frames, and the pipeline registers ten of them.
Strain is computed over part of the cardiac cycle and **the result looks plausible**. 

| Where | Loop as written | Frames covered |
|---|---|---|
| `CropImage` | `for (int frame = 0; frame <= 20; frame++)` | 21, and it skips a frame that is absent |
| `DeformMesh` | `for (int frame = 0; frame < 10; frame++)` | 10, frames 0–9 |
| `GenerateCellAreaStrains` | container-side | 10, frames 0–9 |
| `JacobianThreshold` | container-side | 9, frames 1–9 |

The container adds two more: `uac_pipeline.py` sets `numTimes=10`, and
`GenerateLST_20_TDownsampled.sh` emits ten hardcoded phases naming frames 0, 2, 4 … 18. So the
invariant is **twenty acquired, every other one used**. Preserve it, or change it deliberately in
both places at once.

`DeformMesh` also passes `frame * 10` as the MIRTK `-St` time stamp. That stride matches the second
column of `imgTimes.lst` but is not read from it.

### What to do

**First, write the input contract down** — in this folder, and in the plugin's user manual. How
many frames, the `dcm-<N>.nii` convention, which frame is the reference, and whether every frame or
every other frame is registered.

**Then derive the count once and use it everywhere.** A helper on the view, computed after step 1
and stored as a member:

```cpp
int NumberOfFrames() {
    // dcm-crop-<N>.nii also matches dcm-*.nii, so exclude the frames step 13 writes.
    QStringList frames = QDir(directory + "/nifti/").entryList({"dcm-*.nii"}, QDir::Files);
    return frames.filter(QRegularExpression("^dcm-\\d+\\.nii$")).size();
}
```

The exclusion is needed only until [action 5](#5--rebuild-the-image-as-cemrgafmotion11) moves the
cropped frames out of `nifti/`. Write the filter anyway — it is one line, and without it the count
doubles after step 13 runs.

Derive the `-St` stride from the second column of `tracking/imgTimes.lst` rather than assuming 10,
and add a comment on `JacobianThreshold`'s 1-based start saying that frame 0 is the reference and
its Jacobian is 1 by definition.

### Done when

- One symbol holds the frame count, and all four loops read it.
- A dataset with a frame count other than 21 either runs correctly or stops with a message that
  names the count it found and the count it expected.
- The contract is written down where a user will find it.

---

## 2 — Make a failed `afmotion` step fail visibly

**Why.** `CemrgCommandLine::IsOutputSuccessful` tests `exists() && size() > 0`. Three real failures
pass that test. The 2026-09-03 run passed every check legitimately, which means the checks have
still **never been tested against a real failure**.

| # | Problem | Steps |
|---|---|---|
| 1 | A failed run still writes a small file | 16, 19 |
| 2 | The subcommand saves onto its own input | 12 |
| 3 | No check separates fresh output from stale | all |

**Problem 1 is the concrete one, and the numbers now exist.** Run against an empty project
directory, `generateCellAreaStrains` writes ten 6-byte CSVs holding `,Area` — a pandas header with
no rows — and `calcFiberStrains` writes eighteen 52-byte files holding `Number of Cells RefMsh: 0`
and `Number of Cells TrkMsh: 0`. Both are non-empty. Both pass. Real outputs are 6.45 to 7.90 MB.

### What to do

**A size floor.** `AtrialStrainMotionView::IsOutputFileCorrect` already builds a per-file report
through the file-scope `DescribeOutputProblem`. Give that function a minimum size, and let each call
site pass one. A floor of 1 MB clears every real output by a factor of six and every failure
artefact by four orders of magnitude — the two populations are five orders of magnitude apart. The
reference sizes are in [`pipeline-state.md` §9](pipeline-state.md#9-what-one-full-run-produced).

**A freshness check** for problems 2 and 3. Capture a timestamp before the container call and
compare the output's modification time against it. This is the general fix, not a special case for
step 12's in-place pair — a re-run whose input has gone missing currently passes on the previous
run's file, because output is host-owned and survives.

**Do not rely on checking the last file alone.** Steps 16 and 19 loop in `bash` without `set -e`, so
the highest-numbered file is written even when every frame failed. Check every expected file, or
check the first as well as the last. Step 17 loops in Python and stops at the first error, so its
last-file check is sound.

### Two things this does not fix

`CemrgCommandLine::ExecuteTracking` returns `void` and only logs on failure, so step 14's real
signal — whether `tracking/tsffd.dof` arrived — never reaches the user. Widening it to `bool` is
source-compatible with its three callers (`mmeasurement`, `mmcwplugin`, `CemrgCommandLineTest`), but
it is a change to a shared module and belongs in its own PR.

Step 14's container check is weak by nature. `GenerateCFG.sh` is a `cat` of a fixed template and
`GenerateLST_20_TDownsampled.sh` is ten hardcoded `echo`s. Neither reads any input, so neither can
fail on bad data. A pass shows only that the container ran and could write to `tracking/`. The
useful check there is different: **verify that every frame `imgTimes.lst` names exists on disk**
before handing it to MIRTK. That also catches step 14 being run before step 13, which nothing
detects today.

### Done when

- A deliberately broken run of step 16 or 19 — for example, delete `tracking/cLr-aligned-5.vtp`
  first — reports a failure instead of reporting success.
- Re-running a step whose input has been removed reports a failure rather than passing on the
  previous run's output.

---

## 3 — Delete the dead code, and settle the mitral valve question

**Why.** `AtrialStrainMotionView.cpp` is about 1,970 lines, of which roughly 150 are unreachable and
several members are indeterminate at the moment they are read. A reviewer coming to this plugin for
the first time cannot tell which half matters. This is the cheapest action and it makes every other
change easier to review.

### Verified dead

| Symbol | Evidence |
|---|---|
| `AutomaticAnalysis` | ~100 lines. Declared and defined. **No caller anywhere.** |
| `UserIncludeLgeAnalysis` | called only from `AutomaticAnalysis` |
| `ClipperMV` | declared as a slot, defined, **never connected and never called** |
| `automaticPipeline` | read in exactly one place, the `if (automaticPipeline)` block in `ClipperPV`. `SegmentExtract` forces `uiSelector_pipeline = 2` and then calls `SetAutomaticPipeline(false)` unconditionally, so the block is unreachable |
| `resurfaceMesh` | written by `CheckLoadedMeshQuality`, never read |
| `userHasSetMeshingParams` | written by `GetUserMeshingInputs` and `IdentifyPV`, never read |
| `cnnPath` | cleared by `SetDefaultPanelState`; every other use is inside the two dead methods |
| `m_UISelector`, `m_UIMeshing` | declared in the header with their generated UI includes, never `setupUi`'d |
| `uiUac_meshtype_labelled`, `uiUac_fibreFileIndex`, `uiUac_surftypeIndex` | collected and persisted to `prodUacMetadata.txt`, never read by any downstream call |
| `uiSelector_imgauto_skipCemrgNet`, `uiSelector_imgauto_skipLabel`, `uiSelector_man_useCemrgNet` | persisted to `prodMetadata.txt`; every use is on the dead path |
| `SetUseDockerContainers` | called **18 times** in this view. `_useDockerContainers` is read in exactly one place in the whole module, `CemrgCommandLine::DockerDicom2Nifti`, which this plugin never calls. Every `Execute*` method runs the local `MLib/` binaries unconditionally; every `Docker*` method shells out to `docker` unconditionally |
| `SetDockerImageOpenCarp()` in `UAC_Stage1` and `UAC_Stage2` | `OpenCarpDocker` sets the image itself on its first line |

**A real bug hides in that list.** `SetDefaultPanelState` initialises `tagName`, `cnnPath` and the
four `uiSelector_*` flags — and nothing else. `automaticPipeline`, `analysisOnLge`, `resurfaceMesh`
and `userHasSetMeshingParams` are `bool` members with no initialiser and no constructor. Clicking
step 7 as the first action after opening the view reads `automaticPipeline` uninitialised; clicking
step 4 first reads `analysisOnLge` uninitialised. Deleting the dead members removes most of this.
Give a default to whatever survives.

### The mitral valve decision

The `if (automaticPipeline)` block is not only dead code — it is the reason **the mitral valve is
never clipped**. That block once held the only call to `ClipperMV`, plus the
`CemrgAtrialTools::SetSurfaceLabels` correction that maps naive clipper labels onto the default
11/13/15/17/19 convention. The branch HEAD commit message *"without automated clipping"* records
this as the current intended state.

The mesh therefore retains the MV plane, which is a hole in the anatomy that UAC's Laplace solves
assume is closed. Whether that materially corrupts the coordinates is **untested either way**.
Decide, and record the decision:

- **If MV clipping belongs here:** give it its own button rather than hiding it behind a dead flag.
  Note that `ClipperMV` reads `UAC_CT/prodMVI.nii`, which nothing on the live path produces —
  `CleanAutomaticSegmentation` writes `prodMVI.nii` to the project root, and only from the dead
  `AutomaticAnalysis`. So re-enabling the branch as it stands would fail. The `atrialfibres`
  plugin's interactive route, `AtrialFibresView::SelectMvLandmarks` plus `AtrialFibresView::ClipMV`,
  is the alternative if the automatic version proves unreliable.
- **If it does not:** delete the block, `ClipperMV`, `AutomaticAnalysis`, `UserIncludeLgeAnalysis`,
  `automaticPipeline` and the whole `uiSelector_pipeline` machinery, and state in the user manual
  that MV clipping is handled elsewhere or is not needed.

Either way, remove `uiSelector_pipeline`. A three-way selector that is always forced to one value is
pure confusion, and `GetUserAnalysisSelectorInputs` never shows the dialog it exists for.

### The rest of the tidy-up, in the same pass

- **`files.cmake`.** `UI_FILES` names 16 entries, of which 15 are distinct —
  `AtrialFibresLandmarksViewRough.ui` is listed **twice**. Only eight are needed. Three come from
  this view's own members: `AtrialStrainMotionViewControls.ui`, `AtrialFibresViewUIUacSelector.ui`
  and `AtrialFibresViewUIEditLabels.ui`. Five more come from `AtrialFibresClipperView.h`, which
  `AtrialStrainMotionView.cpp` includes, and which itself includes the generated headers for
  `AtrialFibresClipperViewControls.ui`, `AtrialFibresClipperViewLabels.ui`,
  `AtrialFibresClipperViewUIRadius.ui`, `AtrialFibresClipperViewUICorridor.ui` — **and
  `AtrialFibresViewUIMeshing.ui`**. Note that last one: it stays needed even after `m_UIMeshing` is
  deleted from this view, because the clipper header pulls it in. `AtrialFibresViewUIAnalysisSelector.ui`
  can go with `m_UISelector`. The other six — `AtrialFibresViewControls.ui`,
  `AtrialFibresViewUIRemesh.ui`, `AtrialFibresViewUIConvert.ui`,
  `AtrialFibresLandmarksViewRough.ui`, `AtrialFibresVisualiseViewControls.ui` and the scar plugin's
  `AtrialScarViewUIcemrgnet.ui` — are unused today.
- **Commented-out bodies.** `OnSelectionChanged` is entirely commented out. `GetUserMeshingInputs`
  hardcodes its four parameters and never shows a dialog. `SetLgeAnalysis` carries a commented-out
  `m_Controls.button_z_scar->setEnabled(b)`.
- **A developer's path in a comment.** `UAC_Stage1` opens with a commented `docker run` line naming
  a home directory on the original author's workstation. Replace it with the argument list the code
  actually builds, or delete it.
- **Discarded return values.** `MITK_INFO(<call>) << ...` is MITK's conditional-log idiom, and it is
  used here to invoke functions whose result then decides only whether a line is logged. Two
  instances remain: `MITK_INFO(LoadSurfaceChecks())` in `UacCalculationVerifyLabels`, and
  `MITK_INFO(QFile::copy(...))` in `IdentifyPV` and `CreateModel`, where a failed copy is
  indistinguishable from a successful one. Each should be an `if` that reacts to the result. A third
  instance of the same idiom is why the missing `petsc_opts` folder went unnoticed for four years.
- **Naming.** Rename the slot `PlotFibersTrains` to `PlotFibreStrains`. The Docker argument string
  `plotFiberStrains` must not change. Widget names still carry meaningless `man`/`auto` prefixes
  from `atrialfibres`, and the panel has no button 6 — steps 2 and 6 read "(optional)" instead of a
  number.
- **Documentation.** `documentation/UserManual/Manual.dox` and `documentation/doxygen/modules.dox`
  are unedited MITK templates. The manual still claims the plugin "Generates PhD thesis" and "Brings
  world peace". `AtrialStrainMotionView.h` still carries the generated warning *"This class is not
  yet documented. Use 'git blame' and ask the author to provide basic documentation"*. Replace all
  three; [`pipeline-state.md`](pipeline-state.md) §1 and §5 can be adapted directly.
- **Tests.** Nothing under `Modules/CemrgAppModule/test/` covers this plugin or the two
  `CemrgCommandLine` methods it introduced. Add argument-construction tests for
  `DockerCctaMultilabelSegmentation` and `DockerAtrialStrainMotion`, following the `QtTest` pattern
  there and registering the new class in `test/files.cmake`. Argument building is pure string logic
  and testable without Docker — the `--user` pair and the trailing host path are both easy to break
  by accident.

### Done when

- No unreachable method remains, and every surviving `bool` member has a defined value at the moment
  it is first read.
- `files.cmake` names only the forms the plugin uses, once each.
- The two `.dox` files describe this plugin.
- `ctest -R Cemrg` passes with at least one new test covering `DockerAtrialStrainMotion`.

---

## 4 — Choose step 1's reference frame explicitly

**Why.** `SegmentExtract` picks its input like this:

```cpp
QDir nifti_dir = QDir(directory + "/nifti/");
QStringList nifti_entries = nifti_dir.entryList(QDir::AllEntries | QDir::NoDotAndDotDot);
// ... empty check ...
QString input_file_name = nifti_entries.at(0);
QString input_file_path = directory + "/" + input_file_name;
```

The empty-directory guard was added and is correct. Three problems remain:

1. **`QDir::AllEntries` includes directories.** A stray subfolder in `nifti/` would be picked as the
   input file.
2. **The sort is lexical, and there is no name filter.** With `dcm-1.nii` … `dcm-20.nii`,
   `dcm-10.nii` sorts before `dcm-2.nii`, so which frame becomes the reference is effectively
   arbitrary. On the demo data `dcm-0.nii` wins, and step 13's `dcm-crop-<N>.nii` files sort after
   it only because `'0'` is 0x30 and `'c'` is 0x63. **That is luck.** A series starting at 1, or a
   frame prefix sorting above `dcm-crop-`, silently changes the reference frame or takes a crop as
   the input.
3. **The copy target is the project root**, so `dcm-<N>.nii` and `dcm-<N>_label_maps.nii.gz` sit
   loose beside the pipeline's real outputs. Step 13 does the same with `dcm-10.nii`.

### What to do

```cpp
QStringList frames = QDir(directory + "/nifti/")
                        .entryList({"dcm-*.nii"}, QDir::Files, QDir::Name);
frames = frames.filter(QRegularExpression("^dcm-\\d+\\.nii$"));   // exclude dcm-crop-<N>.nii
```

Then **sort by the numeric suffix, not lexically**, and choose the reference frame explicitly:
either index 0 with a comment saying frame 0 is the reference, or a user choice. This is the same
decision as [action 1](#1--derive-the-frame-count-from-the-data) — do them together.

Give the working copies their own folder rather than the project root. The step-1 and step-13 copies
have nothing to do with each other and neither belongs beside `LA.nii` and `Labelled.nii`.

Keep the existing rule: **the pipeline does not modify the data the user supplies.** Step 13 already
reads each frame through a scratch copy at `<project>/crop-input.nii` for exactly that reason.

### Done when

- A `nifti/` folder holding a subfolder, or holding `dcm-1.nii` … `dcm-20.nii` with no `dcm-0.nii`,
  selects a frame the user can predict from the documented contract.
- Running step 13 and then step 1 selects the same frame step 1 selected before.

---

## 5 — Rebuild the image as `cemrg/afmotion:1.1`

**Read the risk before starting.** `cemrg/afmotion:1.0` is the September 2024 build that Bei Zhou
made on his own workstation, recovered through `docker save` and published unchanged. There is no
newer build, and its author does not recall its state. The reconstructed Dockerfile is complete —
22 layers map 1:1 onto 22 non-empty history entries — **but it needs a VTK 9.3 source tree, and the
exact tag used in 2024 is recorded nowhere in the image.** A rebuild has to pick a 9.3.x tag and
demonstrate the result behaves the same. That is why 1.0 was published rather than rebuilt, and it
is the whole cost of this action.

Everything a rebuild needs, apart from the VTK source, is in `discovery_report/recovered/` inside
the `charlie_pipeline` hand-over bundle: `Dockerfile.reconstructed`, `requirements.lock` (41 pinned
versions, captured from the running container), the whole `/code` payload, and the `strains.cxx` /
`AreaStrain.cxx` sources. **Commit all of it into this repository as part of this action.** Today
the image has no source anywhere in version control, and that is the single largest maintenance risk
in the pipeline.

Wire the lockfile in, because the original `pip install` line was entirely unpinned:

```dockerfile
COPY requirements.lock /tmp/
RUN pip install -r /tmp/requirements.lock
```

### The three changes

**5a. Move the cropped frames out of `nifti/`.** This is the one that pays off in this repository.
`GenerateLST_20_TDownsampled.sh` line 5 is the **only** place in the whole image that names the
cropped frames:

```bash
echo "$1/nifti/dcm-crop- .nii"		#cropped niftis
```

Confirmed by grep over the entire `/code` payload. Because the container hardcodes that path into
`tracking/imgTimes.lst`, and the host-side MIRTK `register` reads that file, the plugin has no
choice but to write into the user's own input folder. Change the path here — to `nifti-cropped/`,
say — and `CropImage` can write to a folder the pipeline owns. Then:

- `CropImage` creates `nifti-cropped/` and clears stale `dcm-crop-*.nii` first. `imgTimes.lst` names
  frames by number, so a second run on a shorter series would otherwise leave frames from the first
  that the list still names.
- Actions [1](#1--derive-the-frame-count-from-the-data) and
  [4](#4--choose-step-1s-reference-frame-explicitly) no longer need their crop filters. Leave the
  filters in; they cost one line and they are correct either way.

**5b. Let an internal failure propagate.** `uac_pipeline.py` dispatches every subcommand through
`os.system(cmd)` and discards the return value — 17 times, with no error handling in the file.
Replace each with `subprocess.run(cmd, shell=True, check=True)`.

> **On its own this gives CemrgApp zero extra signal.** `CemrgCommandLine::ExecuteCommand` never
> inspects `process->exitCode()`. The host-side check is the prerequisite, and it is about three
> lines — but `ExecuteCommand` is shared by every plugin in the repository, so making it react to a
> non-zero exit status changes behaviour well outside this pipeline. Give it **its own PR and its own
> review**, and check the other callers first. Until then, treat 5b as preparation.

Two details a naive conversion gets wrong. The `registration` branch has its MIRTK `register` block
commented out but leaves the following `print(cmd)` and `os.system(cmd)` uncommented, so it
**re-runs the previous command** — `GenerateLST_20_TDownsampled.sh` a second time. It is harmless
today because `>` truncates, so it is idempotent. `check=True` would faithfully preserve the
duplicate. Delete it. And while in the file, delete the dead `MLIB_PATH` constant and the
`deformMesh` branch, or comment that the plugin runs those natively — as written they invite someone
to "fix" the container by installing MIRTK in it, which would break the host-side `register`.

**5c. Parameterise `cropImages`.** `uac_pipeline.py` hardcodes the crop box:

```python
x_start = 35 ; y_start = 33 ; z_start = 70
x_end   = 366; y_end   = 320; z_end   = 244
cmd = f"bash {SCRIPTS_PATH}/CropNiftiLoop.sh {NIFFITY_PATH} {x_start} … {z_end}"
```

**Be clear about what this buys.** `cropImages` is **not on the live path** — step 13 crops on the
host through `CemrgCommonUtils::CropImage`, and has done since PR #97. Parameterising it is
preparation for a future where the container does the crop again, not a fix to anything that runs
today. There is a cheaper option: **delete the subcommand** and record in the container's README
that the host does the crop. Recommended, unless someone wants the container usable standalone.

If you do parameterise it, three faults must be fixed at the same time or it will be wrong in a new
way:

1. **The provenance of `35, 33, 70, 366, 320, 244` is unknown.** Bei does not recall whether they
   came from the upstream `LA-cropping.xlsx` or were measured by hand. Treat them as a constant for
   one dataset, never as a default.
2. **The last three arguments are sizes, not end indices.** `uac_pipeline.py` names them `x_end`,
   `y_end`, `z_end`; `CropNiftiLoop.sh` echoes them as `x_size`; `CropNifti_v3.py` slices
   `img_ar[index_x:index_x+size_x, …]`. The crop actually taken is `[35:401, 33:353, 70:314]`. Fix
   the names or the values.
3. **`CropNifti_v3.py` guesses the sign of the affine shift per axis** — `-index_x*|a00|`,
   `+index_y*|a11|`, `+index_z*|a22|`. That is right for this data's direction matrix and wrong in
   general. `mitk::BoundingObjectCutter`, which the host crop uses, sets the origin with
   `IndexToWorld` on the region start and is correct for any direction matrix.

If the plugin ever passes its own box, it must also **reorder**: `CropImage` measures the box in ITK
bounding-box order, axis-interleaved `[x_min, x_max, y_min, y_max, z_min, z_max]`, while
`cropImages` wants axis-grouped `(x_start, y_start, z_start, x_size, y_size, z_size)`. Passing the
box positionally would put `x_max` into `y_start`.

Also note `CropNiftiLoop.sh` loops `seq 0 1 19` — twenty frames — against the 21 the host loop
covers, and against the ten phases `imgTimes.lst` names. Reconcile that with
[action 1](#1--derive-the-frame-count-from-the-data).

### Done when

- The Dockerfile, the lockfile and the `/code` payload are committed to this repository.
- A build of 1.1 reproduces the 2024 regression fixture in the hand-over bundle: `autoLM` against
  `02_tests_from_gt/01_autoLM_tests/v-0005` writes `prodRefinedLandmarks.txt` with md5
  `14a5ed485da22b3e6a20009a26db8c54`.
- Steps 12, 14 and 16 to 20 run to completion on the `V-0004` demo dataset, and every per-frame
  output matches the reference sizes in
  [`pipeline-state.md` §9](pipeline-state.md#9-what-one-full-run-produced) and holds distinct
  content by `md5sum`.
- `cemrg/afmotion:1.1` is pushed, and `SetDockerImageAfmotion`'s default tag is updated. Keep 1.0
  published — it is the artefact the 2026-09-03 run was verified against.

---

## Lower-priority gaps

Kept short. None of these stopped the 2026-09-03 run. They change results quietly, waste work, or
wait to bite a second dataset.

### In this plugin

| Gap | Effect |
|---|---|
| `AtrialFibresClipperView` holds `directory` and `fileName` as private statics, and its destructor clears them. `IdentifyPV`, `MeshPreprocessing` and `UacCalculationVerifyLabels` all call `ResetPerspective()` **before** `SetDirectoryFile`. BlueBerry does not guarantee the old view is disposed inside `ResetPerspective` — teardown goes through `deleteLater`, and `ShowView` runs a nested event loop that drains the deferred-delete queue | The destructor can run *after* the statics are set. The clipper then falls back to a file-open dialogue and overwrites `automaticPipeline` from whichever button the user clicks, so the view can open in the wrong mode. **Two statics plus a clearing destructor cannot be made safe by ordering** — pass the values to the created view instead. `ShowView` returns the `berry::IViewPart`, so the caller can set them on the instance it receives |
| The `R` branch of `AtrialFibresClipperView::KeyCallBackFunc` passes `automaticPipeline` to **both** `SetAutomaticModeButtons` and `SetManualModeButtons`, where `CreateQtPartControl` passes it and its negation | In manual mode `R` hides both sets and the button row goes empty until the view is reopened. The same branch never re-enables `button_man1_ctrlines`, so point selection cannot be resumed. Pass `!automaticPipeline` to the second call |
| `AtrialFibresClipperView::SetFocus` focuses `button_man1_ctrlines`, which `CreateQtPartControl` has already hidden in automatic mode | The focus lands on whichever widget is next in the tab chain. Branch on the mode |
| The `D` branch of `KeyCallBackFunc` rebuilds the corridor point sets but does not call `Render()`, where the `Delete` branch does | The removed point stays on screen until another event redraws. The user presses `D` again and loses a second point |
| `CemrgAtrialTools::AssignOstiaLabelsToVeins` is called from `CreateLabelledMesh` with bare `directory` | `Labelled.nii` lands in the project root rather than `UAC_CT/`. Harmless today, but a misplaced intermediate is how a later step reads a stale file. Note the inconsistency *inside* `CreateLabelledMesh`: it passes `directory + "/UAC_CT"` to `ProjectTagsOnExistingSurface` and bare `directory` to `AssignOstiaLabelsToVeins` a few lines earlier. A `UacPath()` helper beside the existing inline `Path()` would fix the whole family |
| `GetUserEditLabelsInputs` has a full RA branch, but `UAC_Stage1`, `UAC_Stage2` and `CreateModel` all pass `"la"` literals and the LA tag list `1 19 11 13 15 17` | Choosing RA changes only which fields the dialog shows. An RA dataset produces silently wrong results. Given the pipeline is about *left* atrial strain, removing the RA branch is probably the honest choice |
| Phase B loads nothing into the Data Manager — no `AddToStorage` call appears from `CropImage` onwards — and there is no progress indication outside step 13, though the view includes `mitkProgressBar.h` | After potentially many minutes of computation the user has no feedback and must hunt through `tracking/` by hand. `MmcwViewPlot` and `powertransViewPlot` are the existing in-repo patterns for a plot view |
| Step 13's box is loose on the x axis — see [`pipeline-state.md` §9](pipeline-state.md#9-what-one-full-run-produced) | The crop keeps 396 of 512 voxels on x. A largest-connected-component filter on label 4 before the bounding box would tighten it and protect it from a single stray voxel |
| Step 13 hardcodes `dcm-10.nii` as the frame it measures | Nothing records why frame 10 was chosen. It stays until something does |

### Upstream, in `CemrgAppModule`

Not this plugin's code, but it bites this plugin. Listed so nobody re-diagnoses them from scratch.

| Symbol | Defect | Effect here |
|---|---|---|
| `CemrgAtrialTools::AssignOstiaLabelsToVeins` | applies labels in raster order from a zero-initialised table | vein voxels seen before their ostium are set to **0**, not the ostium label — silently erodes veins in step 4. **The one on this list that affects output quality.** Warrants a two-pass rewrite: collect ostium labels per component, then apply |
| `CemrgCommonUtils::PadImageWithConstant` | creates the output with `SetRegions` + `Allocate` only, so spacing, origin and direction are lost; the `constant` parameter is unimplemented and the fill loop always writes 0 | step 1's `UAC_CT/LA_msh.nii` no longer shares a coordinate frame with `LA.nii`, so it does not overlay in the Data Manager and manual correction at step 2 is harder than it should be. Partly masked, because the mask was rasterised with spacing (1,1,1) anyway |
| `CemrgAtrialTools::SetRIPV` | `inline void SetRIPV(int ri){lspv=ri;}` — writes `lspv` | latent; this plugin never calls the setters |
| `CemrgAtrialTools::SetNaiveSegmentationTags`, `CemrgCommandLine::DockerUacFibreMapping` | declared in the header, never defined | link error if ever called |
| `CemrgCommandLine::DockerUacMainMode` | leaves the output list empty for an unrecognised stage, then calls `outputs.at(0)` | throws on a typo'd stage string |
| `CemrgCommonUtils::OpenCarpParamFileGenerator` | hardcodes `gregion[0].ID[0] = 2`, and can return the literal string `"ERROR_IN_PROCESSING"` as a filename | the region matches no element on a mesh tagged 1, 11, 13, 15, 17, 19. Neither call site checks the error string |
| `CemrgAtrialTools::FixSingleLabelConnectivityInSurface` | `int maxObjectSize = -1` truncates a `double` area | regions differing by less than 1.0 in area compare equal; the first wins |
| `CemrgCommandLine::ExecuteRegistration` | appends `.dof` to the *moving image* path, not the output; pushes moving first and fixed second, the opposite of what the parameter names suggest | only on the dead LGE path |
| `CemrgCommonUtils::Binarise` | assigns to its local parameter copy, so the caller never receives the result | use `ReturnBinarised`, which this plugin correctly does |

`SetRIPV` and the two missing definitions are one-line changes with no behavioural risk. The rest
are worth an issue each rather than a speculative fix.

---

## Decisions already taken — do not re-open without new evidence

| Decision | Why |
|---|---|
| **The cropped frames stay in `nifti/`** until the image is rebuilt | Weighed and rejected on 2026-09-03. The container hardcodes `$1/nifti/dcm-crop- .nii`; moving the files means patching the container's own output before MIRTK reads it — a second step's worth of change plus a verification path, to buy tidiness. [Action 5a](#the-three-changes) is the real fix |
| **The plugin does not use `afmotion cropImages`** | Its box is hardcoded to another dataset, its last three arguments are consumed as sizes while named as end indices, and its affine shift guesses a sign per axis |
| **The originals under `nifti/` are never modified** | The pipeline repairs only files it created itself. This is why step 13 reads each frame through a scratch copy, and it is the rule any new step must follow |
| **The 21-frame loop in step 13 stays** until the input layout is documented | Agreed 2026-08-14. The loop skips a frame that is not on disk, and otherwise matches its neighbours |
| **`afmotion:1.0` was published rather than rebuilt** | The exact image is verified against the 2024 regression fixture, no newer build exists, and a rebuild cannot be shown equivalent without the VTK 9.3 tag |
| **A slim runtime for the image is deferred indefinitely** | Two-thirds of the 7.22 GB is a from-source VTK build kept in the same layer as its build tree. Shrinking it would invalidate the verification for a size win nobody is blocked on |
| **The `strains` upstream repository is closed as won't-fix** | No access. The image's own copy of `strains.cxx`, recovered in the hand-over bundle, is the canonical source |
| **Review every diff against `plugin/atrailstrainmotion`** | Not against `master`. The branch name carries a typo and that is the real name |

---

## One open question worth testing

`main_LA_landmarks.py`, which `autoLM` drives, calls `get_mv_clipper_vtk`. When
`UAC_CT/pvClipper_MV.vtk` is absent, that function **renames** `pvClipper_<id>.vtk` — the entry
labelled `21` — into place, rather than copying it. So a re-run after a partial failure may behave
differently the second time. The 2024 regression fixture ships `pvClipper_MV.vtk` directly and never
reaches this path, but a real `UAC_CT/` produced by steps 5 and 7 does. Worth one deliberate test:
run step 9 twice on the same project and compare the landmarks.
