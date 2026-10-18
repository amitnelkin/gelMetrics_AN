// ===== Batch Particle Tracker Macro =====

inputFolder  = "";    // Folder containing nd2 images
outputFolder = "";    // Folder for the tracked fluorescence-channel TIFFs and Mosaic CSVs

// Particle Tracker parameters
// These parameters worked for our data - user may need to tweak them

radius       = 5;
cutoff       = 0.001;
percentile   = 0.1;
linkrange    = 2;
displacement = 10;

// Multi-channel .nd2 movies (red fluorescence + DIC, or + an empty dark
// channel) are reduced to the fluorescence channel before tracking: the
// channel with the highest max/mean contrast on its first frame (red ~15,
// DIC ~3.5, dark ~2; the lowest-minimum rule picked the dark channel of
// 4PP7 012-recovered-after-corruption). Tracking the opened
// 2-channel hyperstack made Mosaic walk the raw interleaved planes, so its
// "frames" 0-360 were red/DIC planes of only the first 180 time points and it
// detected particles on DIC (analysis_outputs/zigzag_diagnosis). The channel
// is saved as <name>.tif in outputFolder and tracked from there, so Mosaic
// writes Traj_<name>.tif.csv next to it; frames are then time points 0-360.

File.makeDirectory(outputFolder);
list = getFileList(inputFolder);
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".nd2")) {
        print("Processing: " + list[i]);
	setBatchMode(true);
	run("Bio-Formats Importer", "open=" + inputFolder + list[i] + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	getDimensions(w, h, c, z, t);
	red = 1;
	if (c > 1) {
	    best = -1;
	    for (k = 1; k <= c; k++) {
	        Stack.setChannel(k);
	        getStatistics(area, mean, min, max);
	        if (max / mean > best) { best = max / mean; red = k; }
	    }
	    run("Duplicate...", "title=red duplicate channels=" + red);
	}
	print("  channels=" + c + " frames=" + t + " tracked channel=" + red);
	base = replace(list[i], ".nd2", "");
	saveAs("Tiff", outputFolder + base + ".tif");
	run("Close All");
	setBatchMode(false);
	open(outputFolder + base + ".tif");

        // Run Mosaic Particle Tracker
	print("Running tracker on " + base + ".tif");
	run("Particle Tracker 2D/3D", "radius=" + radius + " cutoff=" + cutoff + " per/abs=" + percentile + " link=" + 	linkrange + " displacement=" + displacement + " dynamics=Brownian saveMss");

	run("Close All");
    }
}
print("Batch processing complete!");
