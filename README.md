# @bunny-editor/smooth-fit-curve

Turn freehand drawn points into smooth cubic Bezier curves.

The kind of thing you need when building a pen/pencil tool. Based on Schneider's curve fitting from Graphics Gems (1990). No dependencies.

## Install

```bash
npm install @bunny-editor/smooth-fit-curve
```

## Basic usage

```ts
import { fitCurve } from '@bunny-editor/smooth-fit-curve';

// raw points from mouse/pen input
const points: [number, number][] = [
  [0, 0], [10, 5], [20, 20], [30, 15], [40, 30],
  [50, 25], [60, 40], [70, 35], [80, 50],
];

// lower maxError = more accurate, higher = smoother
const curves = fitCurve(points, 4);

// each curve: [startPoint, controlPoint1, controlPoint2, endPoint]
```

## Get an SVG path

```ts
import { fitCurve, bezierCurvesToSvgPath } from '@bunny-editor/smooth-fit-curve';

const d = bezierCurvesToSvgPath(fitCurve(points, 4));
// "M 0 0 C 5.2 1.3, 9.8 4.1, 10 5 C ..."
```

## Normalize for element positioning

Shifts the path to start near (0,0) and gives you the offset + bounding box — handy when you're storing it as a canvas element:

```ts
import { fitCurve, normalizeFittedCurves } from '@bunny-editor/smooth-fit-curve';

const { path, offsetX, offsetY, width, height } = normalizeFittedCurves(fitCurve(points, 4));
```

## Smoothing slider

If you want to expose a 1–100 slider to your users:

```ts
import { smoothingToError, fitCurve } from '@bunny-editor/smooth-fit-curve';

const maxError = smoothingToError(50); // ~8
const curves = fitCurve(points, maxError);
```

1 = keeps all the wobbles, 50 = balanced, 100 = very smooth with few curves.

## API

| Function | What it does |
|----------|-------------|
| `fitCurve(points, maxError)` | Fit bezier curves to `[x,y][]` points. Returns `BezierCurve[]`. |
| `bezierCurvesToSvgPath(curves)` | Convert to SVG path `d` string. |
| `bezierCurvesBounds(curves)` | Bounding box: `{ minX, minY, maxX, maxY, width, height }`. |
| `normalizeFittedCurves(curves)` | Shift to origin + get offset/dimensions. |
| `smoothingToError(1–100)` | Map slider value to `maxError`. |

`BezierCurve` is `[[x,y], [x,y], [x,y], [x,y]]` — start, cp1, cp2, end.

## How it works

Takes your raw points, assigns parameter values based on distance along the polyline (chord-length parameterization), then finds the best-fit cubic Bezier using least-squares. If the fit isn't good enough, it refines the parameters with Newton-Raphson iteration. If it's still off, it splits at the worst point and fits both halves recursively.

The result is the minimum number of cubic Bezier segments needed to approximate your freehand stroke within the given error tolerance.

## License

MIT
