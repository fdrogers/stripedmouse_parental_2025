run("Measure");
run("Crop");
setBackgroundColor(0, 0, 0);
run("Subtract Background...", "rolling=7");
run("Clear Outside");
run("Select None");
run("16-bit");
run("Invert");
setAutoThreshold("Default");
//run("Threshold...");
