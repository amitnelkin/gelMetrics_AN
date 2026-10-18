// Scratch copy of the repo's BatchParticleTracker.ijm (revision d89bf68).
// Only change: inputFolder is set to the scratch nd2 folder, plus a final exit
// so the headless run returns. Tracker parameters are unchanged.

inputFolder  = "C:/Users/roeeamit-lab/Desktop/AN/GelMetrics/zigzag_scratch/4PP7_nd2/";    // Folder containing nd2 images

// Particle Tracker parameters
// These parameters worked for our data - user may need to tweak them

radius       = 5;
cutoff       = 0.001;
percentile   = 0.1;
linkrange    = 2;
displacement = 10;

list = getFileList(inputFolder);
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".nd2")) {
        print("Processing: " + list[i]);
	run("Bio-Formats Importer", "open=" + inputFolder + list[i] + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename(list[i]);
	getDimensions(w, h, c, z, t);
	print("dims: c=" + c + " z=" + z + " t=" + t + " nSlices=" + nSlices);

        // Run Mosaic Particle Tracker
	print("Running tracker on " + list[i]);
	run("Particle Tracker 2D/3D", "radius=" + radius + " cutoff=" + cutoff + " per/abs=" + percentile + " link=" + 	linkrange + " displacement=" + displacement + " dynamics=Brownian saveMss");

	selectImage(list[i]);
	close;
    }
}
print("Batch processing complete!");
eval("script", "System.exit(0);");
