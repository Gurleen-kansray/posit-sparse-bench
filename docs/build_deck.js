const pptxgen = require("pptxgenjs");

const NAVY = "1E2761";
const NAVY_DARK = "151D49";
const ICE = "CADCFC";
const ICE_SOFT = "EAF1FC";
const WHITE = "FFFFFF";
const AMBER = "D98E2E";
const GRAY = "5B6478";
const LIGHTBG = "F7F9FD";
const CARD = "FFFFFF";

const SERIF = "Cambria";
const SANS = "Calibri";

const W = 13.33, H = 7.5;
const MARGIN = 0.55;

let pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.author = "Gurleen Kaur";
pres.title = "Posit Arithmetic Meets Real-World Sparse Matrices";

function makeShadow() {
  return { type: "outer", color: "1E2761", blur: 8, offset: 3, angle: 90, opacity: 0.12 };
}

function lightSlide() {
  let s = pres.addSlide();
  s.background = { color: LIGHTBG };
  return s;
}

function pageNum(s, n) {
  s.addText(String(n), {
    x: W - 0.9, y: H - 0.45, w: 0.5, h: 0.3, fontFace: SANS, fontSize: 10, color: GRAY,
    align: "right", margin: 0,
  });
}

function header(s, kicker, title) {
  s.addText(kicker.toUpperCase(), {
    x: MARGIN, y: 0.4, w: W - 2 * MARGIN, h: 0.3, fontFace: SANS, fontSize: 11.5,
    color: AMBER, bold: true, charSpacing: 1.5, margin: 0,
  });
  s.addText(title, {
    x: MARGIN, y: 0.68, w: W - 2 * MARGIN, h: 0.65, fontFace: SERIF, fontSize: 28,
    color: NAVY, bold: true, margin: 0,
  });
}

function card(s, x, y, w, h, fill) {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, {
    x, y, w, h, rectRadius: 0.07, fill: { color: fill || CARD }, line: { type: "none" },
    shadow: makeShadow(),
  });
}

function statCallout(s, x, y, w, h, value, label, valueColor) {
  card(s, x, y, w, h, WHITE);
  s.addText(value, {
    x: x + 0.15, y: y + 0.1, w: w - 0.3, h: h * 0.55, fontFace: SERIF, fontSize: 26,
    bold: true, color: valueColor || NAVY, align: "center", valign: "bottom", margin: 0,
  });
  s.addText(label, {
    x: x + 0.15, y: y + h * 0.6, w: w - 0.3, h: h * 0.38, fontFace: SANS, fontSize: 11.5,
    color: GRAY, align: "center", valign: "top", margin: 0,
  });
}

// ---------- Slide 1: Title ----------
{
  let s = pres.addSlide();
  s.background = { color: NAVY };
  s.addShape(pres.shapes.OVAL, { x: 10.6, y: -1.3, w: 4.5, h: 4.5, fill: { color: ICE, transparency: 90 }, line: { type: "none" } });
  s.addShape(pres.shapes.OVAL, { x: 11.6, y: 4.5, w: 3.2, h: 3.2, fill: { color: AMBER, transparency: 88 }, line: { type: "none" } });
  s.addShape(pres.shapes.OVAL, { x: -1.2, y: 5.6, w: 3.0, h: 3.0, fill: { color: ICE, transparency: 92 }, line: { type: "none" } });

  s.addText("POSITS JOURNAL CLUB  ·  IEEE HPEC 2025 THEME  ·  LFX MENTORSHIP", {
    x: MARGIN, y: 1.5, w: 10, h: 0.35, fontFace: SANS, fontSize: 12.5, color: ICE,
    bold: true, charSpacing: 1.5, margin: 0,
  });
  s.addText("Posit Arithmetic Meets\nReal-World Sparse Matrices", {
    x: MARGIN, y: 2.0, w: 10.5, h: 1.9, fontFace: SERIF, fontSize: 44, bold: true,
    color: WHITE, margin: 0, lineSpacingMultiple: 1.05,
  });
  s.addText("Quire-exact CG accuracy on FEM structural matrices — posit16 / posit32 / posit64 vs double64", {
    x: MARGIN, y: 3.95, w: 9.8, h: 0.6, fontFace: SANS, fontSize: 16.5, color: ICE, italic: true, margin: 0,
  });

  s.addText([
    { text: "Gurleen Kaur", options: { bold: true, breakLine: true } },
    { text: "RISC-V High Precision", options: { breakLine: true } },
    { text: "Mentors: Kurt Keville (MIT R&D Labs) · Joshua Gyllinsky", options: {} },
  ], { x: MARGIN, y: 5.15, w: 9, h: 1.0, fontFace: SANS, fontSize: 13.5, color: WHITE, margin: 0, lineSpacing: 20 });

  s.addText("June 17, 2026", {
    x: MARGIN, y: 6.65, w: 4, h: 0.35, fontFace: SANS, fontSize: 11.5, color: ICE, margin: 0,
  });
}

// ---------- Slide 2: Background ----------
{
  let s = lightSlide();
  header(s, "Background", "Posit Arithmetic & the Quire");
  pageNum(s, 2);

  const colW = 5.55, colY = 1.65, colH = 3.7;
  card(s, MARGIN, colY, colW, colH, WHITE);
  card(s, MARGIN + colW + 0.5, colY, colW, colH, WHITE);

  s.addText("IEEE 754", { x: MARGIN + 0.3, y: colY + 0.22, w: colW - 0.6, h: 0.4, fontFace: SERIF, fontSize: 17, bold: true, color: GRAY, margin: 0 });
  s.addText([
    { text: "Rounding error accumulates across every multiply and add in a dot product", options: { bullet: true, breakLine: true } },
    { text: "~17 significant decimal digits at float64 — fixed, not adaptive", options: { bullet: true, breakLine: true } },
    { text: "Multiple exception encodings: qNaN, sNaN, +Inf, −Inf, −0", options: { bullet: true, breakLine: true } },
    { text: "No bitwise-reproducibility guarantee across platforms", options: { bullet: true } },
  ], { x: MARGIN + 0.3, y: colY + 0.75, w: colW - 0.6, h: colH - 1.0, fontFace: SANS, fontSize: 13.5, color: NAVY, margin: 0, paraSpaceAfter: 10 });

  const c2x = MARGIN + colW + 0.5;
  s.addText("Posit + Quire", { x: c2x + 0.3, y: colY + 0.22, w: colW - 0.6, h: 0.4, fontFace: SERIF, fontSize: 17, bold: true, color: NAVY, margin: 0 });
  s.addText([
    { text: "Tapered precision — more fraction bits clustered near ±1", options: { bullet: true, breakLine: true } },
    { text: "Quire: exact fixed-point accumulator (512 bits for posit32)", options: { bullet: true, breakLine: true } },
    { text: "Dot product rounding: 2n−1 roundings (IEEE) → 1 round, at the end", options: { bullet: true, breakLine: true } },
    { text: "Only two exceptions: 0 and NaR (Not a Real)", options: { bullet: true } },
  ], { x: c2x + 0.3, y: colY + 0.75, w: colW - 0.6, h: colH - 1.0, fontFace: SANS, fontSize: 13.5, color: NAVY, margin: 0, paraSpaceAfter: 10 });

  card(s, MARGIN, 5.65, W - 2 * MARGIN, 1.15, NAVY);
  s.addText([
    { text: "The claim we test: ", options: { bold: true, color: AMBER } },
    { text: "Gustafson, HPEC 2025 — exact dot products via the quire mean 16-bit posits can replace 64-bit doubles in stable computations. We test this directly on real FEM sparse matrices, not synthetic benchmarks.", options: { color: WHITE } },
  ], { x: MARGIN + 0.3, y: 5.8, w: W - 2 * MARGIN - 0.6, h: 0.9, fontFace: SANS, fontSize: 14, valign: "middle", margin: 0, lineSpacing: 19 });
}

// ---------- Slide 3: Methodology ----------
{
  let s = lightSlide();
  header(s, "Methodology", "Two FEM Matrices, One Instrumented CG Loop");
  pageNum(s, 3);

  const tx = MARGIN, ty = 1.7, tw = 7.1;
  s.addTable([
    [
      { text: "Matrix", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "Size", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "Nonzeros", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "Dynamic range", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "Condition", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
    ],
    ["bcsstk03 (FEM)", "112 × 112", "640", "~10^16.6", "~10^6"],
    ["bcsstk14 (FEM)", "1,806 × 1,806", "63,454", "~10^15.7", "~10^10"],
  ], {
    x: tx, y: ty, w: tw, colW: [1.7, 1.4, 1.2, 1.4, 1.4],
    fontFace: SANS, fontSize: 13, border: { type: "solid", color: "E2E8F0", pt: 1 },
    fill: { color: WHITE }, autoPage: false, valign: "middle", rowH: 0.55,
  });

  s.addText("SuiteSparse HB structural-stiffness matrices, both real and symmetric positive definite by construction.", {
    x: tx, y: ty + 1.9, w: tw, h: 0.6, fontFace: SANS, fontSize: 12.5, italic: true, color: GRAY, margin: 0,
  });

  const cx = MARGIN + tw + 0.45, cw = W - MARGIN - cx;
  card(s, cx, ty, cw, 3.7, NAVY);
  s.addText("What we measure", { x: cx + 0.3, y: ty + 0.2, w: cw - 0.6, h: 0.4, fontFace: SERIF, fontSize: 16, bold: true, color: AMBER, margin: 0 });
  s.addText([
    { text: "Jacobi-preconditioned CG; dotX (pAp, r·r) instrumented every iteration", options: { bullet: true, breakLine: true } },
    { text: "Logged per call: double64 reference vs posit+quire vs posit naive", options: { bullet: true, breakLine: true } },
    { text: "CG itself always runs in double64 — posit values are logged for comparison only, never drive the solver path here", options: { bullet: true, breakLine: true } },
    { text: "Universal C++ v3.80, WSL2/x86_64 today → RISC-V (Milk-V) next", options: { bullet: true } },
  ], { x: cx + 0.3, y: ty + 0.65, w: cw - 0.6, h: 2.9, fontFace: SANS, fontSize: 12.5, color: WHITE, margin: 0, paraSpaceAfter: 9 });

  s.addText("A third matrix (add32, circuit/SPICE) was tested and excluded — see slide 9.", {
    x: MARGIN, y: 5.75, w: W - 2 * MARGIN, h: 0.4, fontFace: SANS, fontSize: 12.5, color: GRAY, italic: true, margin: 0,
  });
}

// ---------- Ladder slide builder ----------
function ladderSlide(pageN, kicker, title, rows, callouts, footnote) {
  let s = lightSlide();
  header(s, kicker, title);
  pageNum(s, pageN);

  const tx = MARGIN, ty = 1.7, tw = W - 2 * MARGIN;
  const head = [
    { text: "Precision", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
    { text: "Max rel error (quire)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
    { text: "Max rel error (naive)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
    { text: "Verdict", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
  ];
  let body = rows.map(r => [
    { text: r[0], options: { bold: true } }, r[1], r[2],
    { text: r[3], options: { color: r[4] ? "B91C1C" : (r[5] ? "0D7A4F" : NAVY) } },
  ]);

  s.addTable([head, ...body], {
    x: tx, y: ty, w: tw, colW: [1.8, 3.3, 3.3, tw - 8.4],
    fontFace: SANS, fontSize: 13, border: { type: "solid", color: "E2E8F0", pt: 1 },
    fill: { color: WHITE }, autoPage: false, valign: "middle", rowH: 0.5,
  });

  const cw = (tw - 0.6) / 3;
  callouts.forEach((c, i) => {
    statCallout(s, tx + i * (cw + 0.3), 4.25, cw, 1.5, c.value, c.label, c.color);
  });

  if (footnote) {
    s.addText(footnote, { x: tx, y: 6.0, w: tw, h: 0.45, fontFace: SANS, fontSize: 12, italic: true, color: GRAY, margin: 0 });
  }
}

// ---------- Slide 4: bcsstk03 ladder ----------
ladderSlide(4, "Result 1 — FEM Structural", "bcsstk03 — Precision Ladder",
  [
    ["posit8", "6.71×10^8", "2.82×10^10", "overflow — unusable", true, false],
    ["posit16", "1.59×10^2", "6.39×10^2", "NOT viable", true, false],
    ["posit32", "5.44×10^-7", "7.37×10^-6", "14x quire wins — viable", false, true],
    ["posit64", "0.000", "0.000", "exact double", false, true],
  ],
  [
    { value: "159", label: "posit16 max rel error — fails", color: "B91C1C" },
    { value: "5.44×10⁻⁷", label: "posit32 max rel error — viable", color: NAVY },
    { value: "14x", label: "quire improvement over naive at posit32", color: AMBER },
  ]
);

// ---------- Slide 5: bcsstk14 ladder + posit64 ----------
ladderSlide(5, "Result 2 — Larger FEM Matrix", "bcsstk14 — Scaling Up, Verifying posit64",
  [
    ["posit8", "6.87×10^5", "1.36×10^5", "overflow — unusable", true, false],
    ["posit16", "5.26×10^-2", "3.68×10^-1", "7x quire wins — NOT viable", true, false],
    ["posit32", "3.65×10^-8", "7.83×10^-7", "21x quire wins — viable", false, true],
    ["posit64", "0.000", "0.000", "exact double", false, true],
  ],
  [
    { value: "2.55×10⁻⁵", label: "posit32 max rel error (300-iter run)", color: NAVY },
    { value: "4.14×10⁻¹⁵", label: "posit64 max rel error — verified", color: AMBER },
    { value: "10 orders", label: "of magnitude improvement, posit32 → posit64", color: NAVY },
  ],
  "posit64+quire matches double64 machine epsilon exactly — no mixed-precision infrastructure needed for near-double accuracy."
);

// ---------- Slide 6: SLFFEA ----------
{
  let s = lightSlide();
  header(s, "Real-Codebase Validation", "SLFFEA — Beyond Standalone Benchmarks");
  pageNum(s, 6);

  const tx = MARGIN, ty = 1.75, tw = 7.3;
  card(s, tx, ty, tw, 3.55, WHITE);
  s.addText([
    { text: "The precision-ladder results above use standalone harnesses on extracted SuiteSparse matrices. As a follow-up, posit32+quire dotX was patched directly into the conjugate-gradient solver of SLFFEA v1.5 (beam module) — a real, otherwise-unmodified open-source FEM package, not a synthetic test case.", options: { breakLine: true, breakLine: true } },
    { text: "On a cantilever test model, the patched solver converges in 12 CG iterations, with posit32+quire vs double absolute differences in the 1e-9 to 1e-10 range — consistent with the precision-ladder findings above.", options: {} },
  ], { x: tx + 0.3, y: ty + 0.25, w: tw - 0.6, h: 3.0, fontFace: SANS, fontSize: 14, color: NAVY, margin: 0, lineSpacing: 21 });

  const cx = tx + tw + 0.4, cw = W - MARGIN - cx;
  statCallout(s, cx, ty, cw, 1.1, "12", "CG iterations to converge", NAVY);
  statCallout(s, cx, ty + 1.25, cw, 1.1, "864 dof", "cantilever beam, 6 elements", AMBER);
  statCallout(s, cx, ty + 2.5, cw, 1.05, "1e-9–1e-10", "abs diff vs double64", NAVY);
}

// ---------- Slide 7: Convergence boundary ----------
{
  let s = lightSlide();
  header(s, "Result 3 — Solver-Level Correctness", "Quire Keeps the Solver Converging, Not Just Accurate");
  pageNum(s, 7);

  const tx = MARGIN, ty = 1.65, tw = W - 2 * MARGIN;
  s.addText("Experiment: scale bcsstk03's largest diagonal entries to push dynamic range from ~10^6 to ~10^18, then run preconditioned CG under double, posit32+quire, and posit32-naive — both with a double-precision matvec and with a full posit32 matvec (a realistic RISC-V deployment).", {
    x: tx, y: ty, w: tw, h: 0.65, fontFace: SANS, fontSize: 12.5, color: GRAY, margin: 0, lineSpacing: 16,
  });

  const t2y = ty + 0.85;
  s.addTable([
    [
      { text: "Scale", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "Dyn. range", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "quire (dbl-mv)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "naive (dbl-mv)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "quire (p32-mv)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "naive (p32-mv)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
    ],
    ["1e+08", "1.52e+14", "265", "343", "299", { text: "DIV", options: { color: "B91C1C", bold: true } }],
    ["1e+09", "1.52e+15", "257", "280", "330", { text: "313*", options: { color: "B45309", bold: true } }],
    ["1e+10", "1.52e+16", "292", { text: "DIV", options: { color: "B91C1C", bold: true } }, "367", { text: "DIV", options: { color: "B91C1C", bold: true } }],
    ["1e+12", "1.52e+18", "430", { text: "DIV", options: { color: "B91C1C", bold: true } }, "439", { text: "DIV", options: { color: "B91C1C", bold: true } }],
  ], {
    x: tx, y: t2y, w: tw, colW: [1.4, 1.7, 1.95, 1.95, 1.95, 1.95],
    fontFace: SANS, fontSize: 12.5, border: { type: "solid", color: "E2E8F0", pt: 1 },
    fill: { color: WHITE }, autoPage: false, valign: "middle", rowH: 0.45,
  });

  card(s, tx, t2y + 2.4, tw, 1.5, NAVY);
  s.addText([
    { text: "Finding: ", options: { bold: true, color: AMBER, breakLine: true } },
    { text: "naive posit32 diverges (or enters chaotic-residual behavior) by dynamic range ~10^14–10^16. Quire's exact accumulation preserves the CG descent direction all the way to ~10^18 — the first empirical demonstration of this boundary on a real FEM matrix.", options: { color: WHITE } },
  ], { x: tx + 0.3, y: t2y + 2.53, w: tw - 0.6, h: 1.25, fontFace: SANS, fontSize: 13, margin: 0, lineSpacing: 18 });

  s.addText("* 313-iteration \"convergence\" at scale=1e9 was checked and found to be a false positive — chaotic residual oscillation crossing the threshold once, not real convergence. Corrected before reporting.", {
    x: tx, y: t2y + 4.1, w: tw, h: 0.4, fontFace: SANS, fontSize: 11, italic: true, color: GRAY, margin: 0,
  });
}

// ---------- Slide 8: Performance overhead ----------
{
  let s = lightSlide();
  header(s, "Performance", "The RISC-V Motivation");
  pageNum(s, 8);

  const tx = MARGIN, ty = 1.7, tw = 6.5;
  s.addTable([
    [
      { text: "matrix (n)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "double (ns)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "posit32+quire (ns)", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "ratio", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
    ],
    ["bcsstk03 (112)", "80.9", "577,264", "7,132x"],
    ["bcsstk14 (1,806)", "1,042.6", "8,244,849", "7,908x"],
    ["add32 (4,960)*", "2,473.1", "26,422,731", "10,684x"],
  ], {
    x: tx, y: ty, w: tw, colW: [1.9, 1.55, 2.05, 1.0],
    fontFace: SANS, fontSize: 12.5, border: { type: "solid", color: "E2E8F0", pt: 1 },
    fill: { color: WHITE }, autoPage: false, valign: "middle", rowH: 0.5,
  });

  s.addText("* add32 here is a vector-length timing data point only — not a CG accuracy result. See slide 9.", {
    x: tx, y: ty + 2.15, w: tw, h: 0.4, fontFace: SANS, fontSize: 11, italic: true, color: GRAY, margin: 0,
  });

  s.addText([
    { text: "quire/naive ratio ≈ 1.0 across all cases — the overhead is posit conversion (software emulation), not quire accumulation. On hardware with a native posit FPU, this overhead disappears.", options: {} },
  ], { x: tx, y: ty + 2.65, w: tw, h: 1.0, fontFace: SANS, fontSize: 12.5, color: NAVY, margin: 0, lineSpacing: 17 });

  const cx = tx + tw + 0.45, cw = W - MARGIN - cx;
  s.addChart(pres.charts.BAR, [{
    name: "quire / double", labels: ["bcsstk03\n(n=112)", "bcsstk14\n(n=1806)", "add32\n(n=4960)"],
    values: [7132, 7908, 10684],
  }], {
    x: cx, y: ty, w: cw, h: 3.55, barDir: "col",
    chartColors: [NAVY], showTitle: true, title: "Software-emulation overhead (×)",
    titleFontSize: 13, titleColor: NAVY, titleFontFace: SANS,
    chartArea: { fill: { color: WHITE } }, plotArea: { fill: { color: WHITE } },
    catAxisLabelColor: GRAY, valAxisLabelColor: GRAY, catAxisLabelFontSize: 10.5,
    valGridLine: { color: "E2E8F0", size: 0.5 }, catGridLine: { style: "none" },
    showValue: true, dataLabelPosition: "outEnd", dataLabelColor: NAVY, dataLabelFontSize: 11,
    showLegend: false,
  });
}

// ---------- Slide 9: Why add32 excluded ----------
{
  let s = lightSlide();
  header(s, "Scientific Discipline", "Why add32 Was Dropped");
  pageNum(s, 9);

  const tx = MARGIN, ty = 1.7, tw = 7.3;
  card(s, tx, ty, tw, 4.0, WHITE);
  s.addText([
    { text: "add32 (a circuit/SPICE matrix) was initially planned as a third, lower-condition-number domain alongside the two FEM matrices.", options: { breakLine: true, breakLine: true } },
    { text: "Direct inspection of the matrix found it is ", options: {} },
    { text: "not symmetric", options: { bold: true } },
    { text: " — 5,172 off-diagonal pairs where A[i,j] ≠ A[j,i] — and not positive definite.", options: { breakLine: true, breakLine: true } },
    { text: "Conjugate Gradient assumes a symmetric positive-definite matrix; its step-size formula and convergence guarantee don't hold otherwise. Any \"CG accuracy\" number computed on add32 isn't measuring a well-defined quantity — even though the numbers looked clean.", options: { breakLine: true, breakLine: true } },
    { text: "Decision: dropped add32 from every CG-based result rather than keep numbers that looked good but rested on an invalid solver assumption.", options: { bold: true } },
  ], { x: tx + 0.3, y: ty + 0.25, w: tw - 0.6, h: 3.5, fontFace: SANS, fontSize: 13.5, color: NAVY, margin: 0, lineSpacing: 19 });

  const cx = tx + tw + 0.4, cw = W - MARGIN - cx;
  statCallout(s, cx, ty, cw, 1.85, "5,172", "asymmetric off-diagonal\npairs, verified directly", "B91C1C");
  statCallout(s, cx, ty + 2.1, cw, 1.85, "0", "CG accuracy claims kept\nfor a non-SPD matrix", NAVY);
}

// ---------- Slide 10: Gustafson claim + bottom line ----------
{
  let s = lightSlide();
  header(s, "Connection to HPEC 2025", "Testing Gustafson's Claim, Component by Component");
  pageNum(s, 10);

  const tx = MARGIN, ty = 1.65, tw = W - 2 * MARGIN;
  s.addTable([
    [
      { text: "Claim component", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
      { text: "Our evidence", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
    ],
    ["Quire gives exact dot-product accumulation", { text: "Confirmed — accumulation error = 0 by construction on both FEM matrices", options: { color: "0D7A4F" } }],
    ["16-bit posits sufficient for stable computation", { text: "Not supported — posit16 fails by 2–4 orders of magnitude on both matrices tested", options: { color: "B91C1C" } }],
    ["posit32 sufficient for stable computation", { text: "Confirmed — viable on both matrices, 14x–21x better than naive", options: { color: "0D7A4F" } }],
    ["posit64+quire matches double64", { text: "Confirmed — 4.14×10⁻¹⁵, machine epsilon", options: { color: "0D7A4F" } }],
  ], {
    x: tx, y: ty, w: tw, colW: [5.4, tw - 5.4],
    fontFace: SANS, fontSize: 13.5, border: { type: "solid", color: "E2E8F0", pt: 1 },
    fill: { color: WHITE }, autoPage: false, valign: "middle", rowH: 0.62,
  });

  card(s, tx, ty + 3.35, tw, 1.35, NAVY);
  s.addText([
    { text: "Bottom line: ", options: { bold: true, color: AMBER, breakLine: true } },
    { text: "the quire's exact-accumulation mechanism holds exactly as claimed. The 16-bit precision claim does not — for FEM structural matrices with dynamic range ~10^16, posit32 is the practical floor.", options: { color: WHITE } },
  ], { x: tx + 0.3, y: ty + 3.5, w: tw - 0.6, h: 1.1, fontFace: SANS, fontSize: 14.5, margin: 0, lineSpacing: 19, valign: "middle" });
}

// ---------- Slide 11: Roadmap ----------
{
  let s = pres.addSlide();
  s.background = { color: NAVY };
  s.addShape(pres.shapes.OVAL, { x: -1.4, y: -1.6, w: 4.2, h: 4.2, fill: { color: ICE, transparency: 91 }, line: { type: "none" } });
  s.addShape(pres.shapes.OVAL, { x: 11.8, y: 5.4, w: 3.4, h: 3.4, fill: { color: AMBER, transparency: 89 }, line: { type: "none" } });

  s.addText("NEXT STEPS", { x: MARGIN, y: 0.5, w: 10, h: 0.35, fontFace: SANS, fontSize: 12, color: AMBER, bold: true, charSpacing: 1.5, margin: 0 });
  s.addText("Roadmap", { x: MARGIN, y: 0.82, w: 10, h: 0.7, fontFace: SERIF, fontSize: 32, bold: true, color: WHITE, margin: 0 });

  const items = [
    ["1", "Milk-V RISC-V execution", "Run posit32+quire dotX on real RISC-V hardware. Bitwise-identical to x86 per §4.3? Key HPEC paper contribution."],
    ["2", "Ozaki / MTL5 comparison", "Side-by-side: quire vs Ozaki's emulation method on bcsstk03 / bcsstk14, in collaboration with Dr. Omtzigt."],
    ["3", "A second domain, done right", "Replace add32 with a real SPD circuit or sparse matrix, verified symmetric and positive definite before any CG run."],
    ["4", "Scale up", "Validate on a larger matrix (e.g. a bigger circuit or FEM case from SuiteSparse) once SPD is confirmed."],
  ];
  let iy = 1.85;
  items.forEach(([n, t, d]) => {
    s.addShape(pres.shapes.OVAL, { x: MARGIN, y: iy + 0.04, w: 0.5, h: 0.5, fill: { color: AMBER }, line: { type: "none" } });
    s.addText(n, { x: MARGIN, y: iy + 0.04, w: 0.5, h: 0.5, fontFace: SERIF, fontSize: 18, bold: true, color: NAVY_DARK, align: "center", valign: "middle", margin: 0 });
    s.addText(t, { x: MARGIN + 0.75, y: iy, w: 10.8, h: 0.4, fontFace: SANS, fontSize: 15.5, bold: true, color: WHITE, margin: 0 });
    s.addText(d, { x: MARGIN + 0.75, y: iy + 0.4, w: 10.8, h: 0.45, fontFace: SANS, fontSize: 12.5, color: ICE, margin: 0 });
    iy += 1.12;
  });

  card(s, MARGIN, 6.45, W - 2 * MARGIN, 0.7, "16224F");
  s.addText([
    { text: "Target: ", options: { bold: true, color: AMBER } },
    { text: "HPEC 2026 — \"Precision requirements for posit arithmetic in sparse iterative solvers\"   ·   Mentors: Kurt Keville (MIT R&D Labs) · Joshua Gyllinsky", options: { color: WHITE } },
  ], { x: MARGIN + 0.25, y: 6.45, w: W - 2 * MARGIN - 0.5, h: 0.7, fontFace: SANS, fontSize: 12, valign: "middle", margin: 0 });
}

pres.writeFile({ fileName: "Posit_Sparse_Bench_Deck.pptx" }).then(() => {
  console.log("DONE");
});
