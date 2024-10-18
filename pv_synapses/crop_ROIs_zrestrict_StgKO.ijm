// ************* FIJI SCRIPT FOR AUTOMATICALLY CROPPING MANUALLY-IDENTIFIED ROIS ************* \\
// Last modified by MES on 6/8/24
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

//CHANGE ALL OF THE BELOW AS NEEDED
parent = "A:/Margaret/mExR/2023.10_StgKO_TRN/registered/round1_noSMIGFAP/"; //where whole field of view images are stored
target = "A:/Margaret/mExR/2023.10_StgKO_TRN/cropped_synaptic_rois/"; //where to save cropped images
roiparent = "A:/Margaret/mExR/2023.10_StgKO_TRN/ROIs/"; //where the .zip files for Fiji ROIs are saved

//roitypes = newArray("PV","PN");

ROIs = newArray( //fields of view to process
	"mouseA-ROI1",
	"mouseA-ROI1-2",
	"mouseA-ROI2",
	"mouseA-ROI3",
	"mouseA-ROI4",
	"mouseA-ROI5",
	
	"mouseB-ROI1",
	"mouseB-ROI2",
	"mouseB-ROI3",
	"mouseB-ROI4",
	"mouseB-ROI5",
	
	"mouseD-ROI1",
	"mouseD-ROI2",
	"mouseD-ROI4",
	"mouseD-ROI5",
	
	"mouseE-ROI2",
	"mouseE-ROI3",
	"mouseE-ROI5",
	
	"mouseG-ROI1",
	"mouseG-ROI2",
	"mouseG-ROI3",
	"mouseG-ROI4",
	"mouseG-ROI5"
	
	//"mouseA-ROI1",
	//"mouseA-ROI1-2",
	//"mouseA-ROI1-3",
	//"mouseA-ROI2",
	//"mouseA-ROI3-1",
	//"mouseA-ROI3-2",
	
	//"mouseB-ROI1",
	//"mouseB-ROI2",
	//"mouseB-ROI3",
	//"mouseB-ROI4",
	
	//"mouseD-ROI4",
	//"mouseD-ROI5",
	
	//"mouseE-ROI1",
	//"mouseE-ROI3",
	//"mouseE-ROI5",
	
	//"mouseF-ROI2",
	//"mouseF-ROI4",
	
	//"mouseG-ROI1",
	//"mouseG-ROI2",
	//"mouseG-ROI3",
	//"mouseG-ROI4",
	//"mouseG-ROI5"	
	
	);
				
nslices_range = 15; //number of z-slices in each direction from ROI midpoint to crop. total z-stack length for each ROI will be double this

// NO NEED TO MODIFY THE BELOW
image_folder = parent;
list = getFileList(image_folder);									 //find all files in current data folder					
extension = "tif";	
run("Close All");
for  (i=0; i<ROIs.length; i++){
	ROI = ROIs[i];

	targetdir = target;// + roitype + "/";
	
	for (j=0; j<list.length; j++) {		//loop through all images in the image folder
		if (endsWith(toLowerCase(list[j]), "." + extension)) {
			if (indexOf(list[j], ROI) >= 0){
				name =  list[j];
		    	aIndex = indexOf(name, ".tif");										//record position of dot that separates extension from filename	
		    	title = substring(name, 0, aIndex); //record the abbreviated file name (excluding file extension; eg:'stack1')
				
				roidir = roiparent;
			    roiManager("Open", roidir + ROI + ".zip"); //open the manually-identified ROIs. NOTE: they must be named the same as the prefix for the .tif stacks
			    open(parent + name);
			    getDimensions(width, height, channels, slices, frames);
			    
				// Iterate all ROIs in ROI Manager
			    for (n=0; n<roiManager("count"); ++n) {
			    	roiManager("Select", n);
			    	zidx = getSliceNumber();
			    	zrange1 = zidx-nslices_range;
			    	zrange2 = zidx+nslices_range;

			    	if (zrange1<1){
			    		zrange1=1;
			    	}
			    	if (zrange2>slices){
			    		zrange=slices;
			    	}

			    	getSelectionBounds(x,y,width,height);
			    	
			    	run("Duplicate...", "duplicate range=" + toString(zrange1) + "-" + toString(zrange2));
			        //saveAs("Tiff", target + title + "_ROI_"+toString(x)+"_"+toString(y)+".tif");
			        if (n+1 < 10){
			        	roino = "0" + toString(n+1);
			        } else {
			        	roino = toString(n+1);
			        }
			        title =  replace(title, "round002", "round001");//NOTE: replaced round001 with GFAP/SMI with round001 WITHOUT GFAP/SMI (has Cav2.1 alone)
			        saveAs("Tiff", targetdir + title + "_ROI"+roino+".tif");
			        close();
			    }
			    selectWindow(name);
			    run("Close");
			    roiManager("Delete");
			    selectWindow("ROI Manager");
			    run("Close");

			}
		}
	}
}
run("Close All");
run("Collect Garbage");
