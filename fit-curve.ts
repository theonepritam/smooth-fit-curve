// Attempt at implementing Schneider's curve fitting algorithm from Graphics Gems.
// Takes a bunch of points and turns them into smooth cubic beziers.

type Point = [number, number];
export type BezierCurve = [Point, Point, Point, Point]; // [start, cp1, cp2, end]

// --- vector math ---

function subtract(a: Point, b: Point): Point {
  return [a[0] - b[0], a[1] - b[1]];
}

function add(a: Point, b: Point): Point {
  return [a[0] + b[0], a[1] + b[1]];
}

function scale(v: Point, s: number): Point {
  return [v[0] * s, v[1] * s];
}

function dot(a: Point, b: Point): number {
  return a[0] * b[0] + a[1] * b[1];
}

function norm(v: Point): number {
  return Math.hypot(v[0], v[1]);
}

function normalize(v: Point): Point {
  const len = norm(v);
  if (len === 0) return [0, 0];
  return [v[0] / len, v[1] / len];
}

// --- bezier math ---

/** evaluate cubic bezier at parameter t */
function bezierQ(ctrlPoly: BezierCurve, t: number): Point {
  const s = 1 - t;
  const a = scale(ctrlPoly[0], s * s * s);
  const b = scale(ctrlPoly[1], 3 * s * s * t);
  const c = scale(ctrlPoly[2], 3 * s * t * t);
  const d = scale(ctrlPoly[3], t * t * t);
  return add(add(a, b), add(c, d));
}

/** first derivative at t */
function bezierQPrime(ctrlPoly: BezierCurve, t: number): Point {
  const s = 1 - t;
  const a = scale(subtract(ctrlPoly[1], ctrlPoly[0]), 3 * s * s);
  const b = scale(subtract(ctrlPoly[2], ctrlPoly[1]), 6 * s * t);
  const c = scale(subtract(ctrlPoly[3], ctrlPoly[2]), 3 * t * t);
  return add(add(a, b), c);
}

/** second derivative at t */
function bezierQPrimePrime(ctrlPoly: BezierCurve, t: number): Point {
  return add(
    scale(add(subtract(ctrlPoly[2], scale(ctrlPoly[1], 2)), ctrlPoly[0]), 6 * (1 - t)),
    scale(add(subtract(ctrlPoly[3], scale(ctrlPoly[2], 2)), ctrlPoly[1]), 6 * t),
  );
}

// --- core algorithm ---

function createTangent(a: Point, b: Point): Point {
  return normalize(subtract(a, b));
}

/** chord-length parameterization — space out t values based on distance between points */
function chordLengthParameterize(points: Point[]): number[] {
  const u: number[] = [0];
  for (let i = 1; i < points.length; i++) {
    u.push(u[i - 1] + norm(subtract(points[i], points[i - 1])));
  }
  const totalLen = u[u.length - 1];
  if (totalLen === 0) return u.map(() => 0);
  return u.map((v) => v / totalLen);
}

/** newton-raphson: nudge each t value closer to its best fit */
function newtonRaphsonRootFind(bez: BezierCurve, point: Point, u: number): number {
  const d = subtract(bezierQ(bez, u), point);
  const qp = bezierQPrime(bez, u);
  const numerator = dot(d, qp);
  const denominator = dot(qp, qp) + 2 * dot(d, bezierQPrimePrime(bez, u));
  if (denominator === 0) return u;
  return u - numerator / denominator;
}

function reparameterize(bez: BezierCurve, points: Point[], parameters: number[]): number[] {
  return parameters.map((p, i) => newtonRaphsonRootFind(bez, points[i], p));
}

/** map t values to distances along the curve (for better error measurement) */
function mapTtoRelativeDistances(bez: BezierCurve, parts: number): number[] {
  const dist = [0];
  let prev = bez[0];
  let sumLen = 0;
  for (let i = 1; i <= parts; i++) {
    const curr = bezierQ(bez, i / parts);
    sumLen += norm(subtract(curr, prev));
    dist.push(sumLen);
    prev = curr;
  }
  if (sumLen === 0) return dist.map(() => 0);
  return dist.map((d) => d / sumLen);
}

function findT(param: number, tDistMap: number[], parts: number): number {
  if (param <= 0) return 0;
  if (param >= 1) return 1;
  for (let i = 1; i <= parts; i++) {
    if (param <= tDistMap[i]) {
      const tMin = (i - 1) / parts;
      const tMax = i / parts;
      const lenMin = tDistMap[i - 1];
      const lenMax = tDistMap[i];
      return ((param - lenMin) / (lenMax - lenMin)) * (tMax - tMin) + tMin;
    }
  }
  return 1;
}

/** least-squares: find the best control points for a cubic bezier through these points */
function generateBezier(
  points: Point[],
  parameters: number[],
  leftTangent: Point,
  rightTangent: Point,
): BezierCurve {
  const firstPoint = points[0];
  const lastPoint = points[points.length - 1];
  const bezCurve: BezierCurve = [firstPoint, [0, 0], [0, 0], lastPoint];

  // Compute the A matrices
  const A: [Point, Point][] = parameters.map((u) => {
    const ux = 1 - u;
    return [
      scale(leftTangent, 3 * u * ux * ux),
      scale(rightTangent, 3 * ux * u * u),
    ];
  });

  // Create C and X matrices
  const C: [[number, number], [number, number]] = [[0, 0], [0, 0]];
  const X: [number, number] = [0, 0];

  for (let i = 0; i < points.length; i++) {
    const u = parameters[i];
    const a = A[i];
    C[0][0] += dot(a[0], a[0]);
    C[0][1] += dot(a[0], a[1]);
    C[1][0] += dot(a[0], a[1]);
    C[1][1] += dot(a[1], a[1]);

    const tmp = subtract(
      points[i],
      bezierQ([firstPoint, firstPoint, lastPoint, lastPoint], u),
    );
    X[0] += dot(a[0], tmp);
    X[1] += dot(a[1], tmp);
  }

  // Compute determinants
  const detC0C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
  const detC0X = C[0][0] * X[1] - C[1][0] * X[0];
  const detXC1 = X[0] * C[1][1] - X[1] * C[0][1];

  // Derive alpha values
  const alphaL = detC0C1 === 0 ? 0 : detXC1 / detC0C1;
  const alphaR = detC0C1 === 0 ? 0 : detC0X / detC0C1;

  // If alpha negative, use the Wu/Barsky heuristic
  const segLength = norm(subtract(firstPoint, lastPoint));
  const epsilon = 1.0e-6 * segLength;

  if (alphaL < epsilon || alphaR < epsilon) {
    bezCurve[1] = add(firstPoint, scale(leftTangent, segLength / 3));
    bezCurve[2] = add(lastPoint, scale(rightTangent, segLength / 3));
  } else {
    bezCurve[1] = add(firstPoint, scale(leftTangent, alphaL));
    bezCurve[2] = add(lastPoint, scale(rightTangent, alphaR));
  }

  return bezCurve;
}

/** how far off is our curve from the original points? returns max distance + where */
function computeMaxError(
  points: Point[],
  bez: BezierCurve,
  parameters: number[],
): [number, number] {
  const B_PARTS = 10;
  const tDistMap = mapTtoRelativeDistances(bez, B_PARTS);
  let maxDist = 0;
  let splitPoint = Math.floor(points.length / 2);

  for (let i = 0; i < points.length; i++) {
    const t = findT(parameters[i], tDistMap, B_PARTS);
    const v = subtract(bezierQ(bez, t), points[i]);
    const dist = v[0] * v[0] + v[1] * v[1];
    if (dist > maxDist) {
      maxDist = dist;
      splitPoint = i;
    }
  }

  return [maxDist, splitPoint];
}

/** try to fit one bezier. if it doesn't work, split and try both halves. */
function fitCubic(
  points: Point[],
  leftTangent: Point,
  rightTangent: Point,
  error: number,
): BezierCurve[] {
  const MAX_ITERATIONS = 20;

  // Two-point heuristic
  if (points.length === 2) {
    const dist = norm(subtract(points[0], points[1])) / 3;
    return [[
      points[0],
      add(points[0], scale(leftTangent, dist)),
      add(points[1], scale(rightTangent, dist)),
      points[1],
    ]];
  }

  // Parameterize and attempt to fit
  let u = chordLengthParameterize(points);
  let bezCurve = generateBezier(points, u, leftTangent, rightTangent);
  let [maxError, splitPoint] = computeMaxError(points, bezCurve, u);

  if (maxError === 0 || maxError < error) {
    return [bezCurve];
  }

  // If error not too large, try reparameterization
  if (maxError < error * error) {
    let uPrime = u;
    let prevErr = maxError;
    let prevSplit = splitPoint;

    for (let i = 0; i < MAX_ITERATIONS; i++) {
      uPrime = reparameterize(bezCurve, points, uPrime);
      bezCurve = generateBezier(points, uPrime, leftTangent, rightTangent);
      [maxError, splitPoint] = computeMaxError(points, bezCurve, u);

      if (maxError < error) {
        return [bezCurve];
      }

      if (splitPoint === prevSplit) {
        const errChange = maxError / prevErr;
        if (errChange > 0.9999 && errChange < 1.0001) break;
      }

      prevErr = maxError;
      prevSplit = splitPoint;
    }
  }

  // Fitting failed — split at max error point and fit recursively
  const centerVector = subtract(points[splitPoint - 1], points[splitPoint + 1]);
  if (centerVector[0] === 0 && centerVector[1] === 0) {
    // Points are the same; rotate tangent 90°
    const diff = subtract(points[splitPoint - 1], points[splitPoint]);
    centerVector[0] = -diff[1];
    centerVector[1] = diff[0];
  }
  const toCenterTangent = normalize(centerVector);
  const fromCenterTangent = scale(toCenterTangent, -1);

  return [
    ...fitCubic(points.slice(0, splitPoint + 1), leftTangent, toCenterTangent, error),
    ...fitCubic(points.slice(splitPoint), fromCenterTangent, rightTangent, error),
  ];
}

// --- public api ---

/**
 * Fit bezier curves to freehand points.
 * maxError controls the tradeoff — small = accurate, large = smoother.
 */
export function fitCurve(points: Point[], maxError: number): BezierCurve[] {
  // Remove duplicate consecutive points
  const cleaned = points.filter(
    (p, i) => i === 0 || p[0] !== points[i - 1][0] || p[1] !== points[i - 1][1],
  );

  if (cleaned.length < 2) return [];

  const leftTangent = createTangent(cleaned[1], cleaned[0]);
  const rightTangent = createTangent(cleaned[cleaned.length - 2], cleaned[cleaned.length - 1]);

  return fitCubic(cleaned, leftTangent, rightTangent, maxError);
}

/** turn fitted curves into an SVG path d string */
export function bezierCurvesToSvgPath(curves: BezierCurve[]): string {
  if (curves.length === 0) return '';

  let d = `M ${curves[0][0][0]} ${curves[0][0][1]}`;
  for (const curve of curves) {
    d += ` C ${curve[1][0]} ${curve[1][1]}, ${curve[2][0]} ${curve[2][1]}, ${curve[3][0]} ${curve[3][1]}`;
  }
  return d;
}

/** bounding box of the fitted curves (sampled + control points) */
export function bezierCurvesBounds(curves: BezierCurve[]): {
  minX: number; minY: number; maxX: number; maxY: number;
  width: number; height: number;
} {
  let minX = Infinity, minY = Infinity;
  let maxX = -Infinity, maxY = -Infinity;

  const SAMPLES_PER_CURVE = 20;
  for (const curve of curves) {
    for (let i = 0; i <= SAMPLES_PER_CURVE; i++) {
      const t = i / SAMPLES_PER_CURVE;
      const p = bezierQ(curve, t);
      if (p[0] < minX) minX = p[0];
      if (p[1] < minY) minY = p[1];
      if (p[0] > maxX) maxX = p[0];
      if (p[1] > maxY) maxY = p[1];
    }
  }

  // Also check exact control points for safety
  for (const curve of curves) {
    for (const p of curve) {
      if (p[0] < minX) minX = p[0];
      if (p[1] < minY) minY = p[1];
      if (p[0] > maxX) maxX = p[0];
      if (p[1] > maxY) maxY = p[1];
    }
  }

  return {
    minX, minY, maxX, maxY,
    width: maxX - minX,
    height: maxY - minY,
  };
}

/** shift curves so the path starts near (0,0). gives you offset + dimensions. */
export function normalizeFittedCurves(curves: BezierCurve[]): {
  curves: BezierCurve[];
  path: string;
  offsetX: number;
  offsetY: number;
  width: number;
  height: number;
} {
  if (curves.length === 0) {
    return { curves: [], path: '', offsetX: 0, offsetY: 0, width: 0, height: 0 };
  }

  const bounds = bezierCurvesBounds(curves);
  const PAD = 1; // 1px padding to avoid clipping strokes
  const ox = bounds.minX - PAD;
  const oy = bounds.minY - PAD;

  const normalized: BezierCurve[] = curves.map((curve) => [
    [curve[0][0] - ox, curve[0][1] - oy],
    [curve[1][0] - ox, curve[1][1] - oy],
    [curve[2][0] - ox, curve[2][1] - oy],
    [curve[3][0] - ox, curve[3][1] - oy],
  ]);

  return {
    curves: normalized,
    path: bezierCurvesToSvgPath(normalized),
    offsetX: ox,
    offsetY: oy,
    width: bounds.width + PAD * 2,
    height: bounds.height + PAD * 2,
  };
}

/** map a 1-100 slider to maxError. 1 = wobbles preserved, 100 = super smooth. */
export function smoothingToError(smoothing: number): number {
  // Clamp to 1–100
  const s = Math.max(1, Math.min(100, smoothing));
  // Exponential mapping: 1 → ~0.5, 50 → ~8, 100 → ~200
  return 0.5 * Math.pow(1.06, s);
}
