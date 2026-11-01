# Algorithm overview

`algorithm-overview.svg` is the editable vector source for the README splash.
It follows the teal/slate visual vocabulary of the vdjtools analysis atlas,
with original schematic artwork. The neighbourhood and histogram do not show
measured observations. Receptor annotation and repertoire inference retain
their distinct scoring, null and calibration contracts.

Render the GitHub/PyPI-compatible raster with:

```sh
rsvg-convert assets/algorithm-overview.svg -o assets/algorithm-overview.png
```

The manuscript uses its existing Graphviz source and build targets separately;
this README image is sized for screen reading.
