setOption("BlackBackground", false);
run("Convert to Mask");
run("Despeckle");
run("Watershed");
run("Analyze Particles...", "size=0.00004-Infinity display exclude summarize");
run("Revert");