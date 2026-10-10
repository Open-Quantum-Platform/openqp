# Art Rendering

Art produces a progressively ray-traced image from the currently selected
Analysis content. Open a project in Analysis and select the desired molecular
structure or 3D map before moving to Art.

## Content

The **Content** menu can include:

- the molecule alone;
- the selected SCF molecular orbital or Dyson orbital;
- NTO hole or particle orbitals;
- attachment, detachment, difference, transition, or state density; and
- an existing or combined cube surface.

Items appear only when the selected project contains the required data. If an
NTO or density cannot be derived from the exported results, the item is marked
unavailable with the reason rather than rendering a substitute surface.

## Rendering controls

**Quality** changes pixel ratio and progressive sample count. Draft responds
quickly while adjusting the view; Publication accumulates more samples and
takes longer.

**Background** offers Studio dark, neutral gray, white, and black. **Exposure**
controls overall brightness. **Depth of field** is off, subtle, or strong;
strong depth of field is an artistic effect and can hide scientific detail.

For volumetric content, choose matte, gloss, or translucent surface material,
opacity, and positive/negative phase visibility. Phase colors and the selected
isovalue originate in Analysis.

## Interaction and export

Drag to rotate, use the normal pointer controls to frame the molecule, and
select **Reset view** to restore the camera. Progressive sampling restarts after
the view or a rendering setting changes. **Pause** freezes or resumes sample
accumulation. Select **PNG** after the image has converged visually.

Art is intended for presentation graphics. Numeric interpretation, isovalues,
state labels, and orbital identities should be checked in Analysis before an
image is exported.
