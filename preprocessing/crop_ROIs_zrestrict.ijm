// ************* FIJI SCRIPT FOR AUTOMATICALLY CROPPING MANUALLY-IDENTIFIED ROIS ************* \\
// Last modified by MES on 12/15/22
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

//CHANGE ALL OF THE BELOW AS NEEDED
parent = "E:/Margaret/mExR/2022.10_5xFAD/cropped_z/"; //where whole field of view images are stored
target = "E:/Margaret/mExR/2022.10_5xFAD/cropped_abeta/"; //where to save cropped images
roidir = "E:/Margaret/mExR/2022.10_5xFAD/abeta_ROIs/"; //where the .zip files for Fiji ROIs are saved

ROIs = newArray( //fields of view to process
	"5xFAD-ctx-ROI1",
	"5xFAD-ctx-ROI2",
	"5xFAD-ctx-ROI3",
	"5xFAD-ctx-ROI4",
	"5xFAD-ctx-ROI5"	
				);
				
nslices_range = 20; //number of z-slices in each direction from ROI midpoint to crop. total z-stack length for each ROI will be double this

// NO NEED TO MODIFY THE BELOW
image_folder = parent;
list = getFileList(image_folder);									 //find all files in current data folder					
extension = "tif";	
run("Close All");
for  (i=0; i<ROIs.length; i++){
	ROI = ROIs[i];
	
	for (j=0; j<list.length; j++) {		//loop through all images in the image folder
		if (endsWith(toLowerCase(list[j]), "." + extension)) {
			if (indexOf(list[j], ROI) >= 0){
				name =  list[j];
		    	aIndex = indexOf(name, ".tif");										//record position of dot that separates extension from filename	
		    	title = substring(name, 0, aIndex); //record the abbreviated file name (excluding file extension; eg:'stack1')
		
			    roiManager("Open", roidir + ROI + "_rois.zip"); //open the manually-identified ROIs. NOTE: they must be named the same as the prefix for the .tif stacks
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
			        saveAs("Tiff", target + title + "_ROI_"+toString(x)+"_"+toString(y)+".tif");
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
